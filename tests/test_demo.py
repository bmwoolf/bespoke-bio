"""Sequence handoff, cardinality, retry behavior and isolated artifact access."""

import json
from datetime import UTC, datetime, timedelta
from hashlib import sha256
from io import BytesIO
from urllib.error import HTTPError

import pytest
from fastapi.testclient import TestClient

import bespoke_bio.demo.service as demo_service
from bespoke_bio.demo.catalog import Catalog, build_index
from bespoke_bio.demo.client import DesignService, ServiceError
from bespoke_bio.demo.published import PublishedSnapshot, artifact_key, export_snapshot
from bespoke_bio.demo.service import create_app
from bespoke_bio.demo.structures import StructureLibrary, publish_structures
from bespoke_bio.demo.vault import materialize
from bespoke_bio.stage1.io import artifact

HEADERS = {"X-Bespoke-Request": "1"}
REMOTE = "a" * 32


def test_runtime_credential_is_private_and_expired_refresh_preserves_valid_file(tmp_path):
    path = tmp_path / "credential.json"
    config = {
        "base_url": "https://design.example",
        "api_token": "synthetic-test-token",
        "expires_at": (datetime.now(UTC) + timedelta(hours=1)).isoformat(),
    }
    materialize(json.dumps(config), path)
    assert path.stat().st_mode & 0o777 == 0o600
    assert DesignService.from_config(path).available
    original = path.read_bytes()
    config["expires_at"] = "2000-01-01T00:00:00+00:00"
    with pytest.raises(ServiceError, match="expired"):
        materialize(json.dumps(config), path)
    assert path.read_bytes() == original
    assert list(tmp_path.iterdir()) == [path]
    config["expires_at"] = "2030-01-01T00:00:00"
    with pytest.raises(ValueError, match="timezone"):
        materialize(json.dumps(config), path)
    assert path.read_bytes() == original


def test_expired_credential_cannot_submit_or_mark_a_prepared_job_dispatched(index, tmp_path):
    runs = tmp_path / "runs"
    with TestClient(create_app(index, runs, FakeService())) as client:
        job = prepare(client)
    expired = DesignService("https://design.example", "fixture", expires_at="2000-01-01T00:00:00Z")
    with TestClient(create_app(index, runs, expired, allow_run=True)) as client:
        metadata = client.get("/api/catalog").json()
        assert metadata["service_configured"] is False
        assert metadata["gpu_execution_enabled"] is False
        assert metadata["service_expires_at"]
        assert "fixture" not in json.dumps(metadata)
        for path, payload in [
            ("/api/prepare", {"target_id": "fixture-0"}),
            (f"/api/jobs/{job['job_id']}/run", {}),
        ]:
            assert client.post(path, json=payload, headers=HEADERS).status_code == 503
    saved = json.loads((runs / (job["job_id"] + ".json")).read_text())
    assert saved["status"] == "prepared_pending_gpu_execution"
    assert not saved.get("dispatch_attempted")


def test_public_operator_link_does_not_enable_mutations(index, published):
    with TestClient(create_app(index, published=published[0], operator_access=True)) as client:
        meta = client.get("/api/catalog").json()
        assert meta["operator_access"] is True
        assert meta["read_only"] is True
        assert meta["service_configured"] is False
        assert meta["service_expires_at"] is None
        assert client.post("/api/prepare", json={}, headers=HEADERS).status_code == 403


@pytest.fixture
def structure_source(index, tmp_path):
    from io import StringIO

    from Bio.PDB import MMCIFIO, PDBParser

    folder = tmp_path / "structure-source"
    folder.mkdir()
    lines = []
    for residue, name in enumerate(("MET", "GLN"), 1):
        for offset, atom in enumerate(("N", "CA", "C", "O")):
            serial = (residue - 1) * 4 + offset + 1
            lines.append(
                f"ATOM  {serial:5d} {atom:^4s} {name} A{residue:4d}    "
                f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{80.0:6.2f}          "
                f"{atom[0]:>2s}  \n"
            )
    writer = MMCIFIO()
    writer.set_structure(PDBParser(QUIET=True).get_structure("fixture", StringIO("".join(lines))))
    writer.save(str(folder / "target.cif"))
    (folder / "confidence.json").write_text('{"synthetic":true}')
    source = {
        "version": 1,
        "catalog_sha256": Catalog(index).metadata()["manifest_sha256"],
        "structures": [
            {
                "amino_acid_sequence": "MQ",
                "sequence_sha256": sha256(b"MQ").hexdigest(),
                "file": "target.cif",
                "file_sha256": sha256((folder / "target.cif").read_bytes()).hexdigest(),
                "confidence": {
                    "file": "confidence.json",
                    "sha256": sha256((folder / "confidence.json").read_bytes()).hexdigest(),
                },
                "source": "Synthetic test",
                "mean_plddt": 80.0,
                "provenance": {"synthetic": True, "peptide_geometry_screen_flags": 1},
            }
        ],
    }
    (folder / "manifest.json").write_text(json.dumps(source))
    return folder / "manifest.json"


def test_structures_preserve_all_records_and_filter_exact_ready_sequences(
    index, structure_source, tmp_path, package
):
    root = tmp_path / "structures"
    assert publish_structures(index, structure_source, root) == {
        "structured_records": 2,
        "structured_sequences": 1,
    }
    with TestClient(create_app(index, structures=root)) as client:
        assert client.get("/api/catalog").json()["target_records"] == 3
        assert client.get("/api/targets").json()["total"] == 3
        pages = [client.get(f"/api/targets?ready=true&limit=1&offset={i}").json() for i in range(3)]
        assert all(p["total"] == 2 for p in pages)
        assert [r["target_id"] for p in pages for r in p["records"]] == ["fixture-1", "fixture-2"]
        assert client.get("/api/targets?ready=true&q=missing").json()["total"] == 0
        assert client.get("/api/targets?ready=true&changed=true").json()["total"] == 2
        assert client.get("/api/targets/fixture-1").json() == package[1][1]
        a = client.get("/api/targets/fixture-1/structure").json()
        b = client.get("/api/targets/fixture-2/structure").json()
        assert a == b
        assert a["geometry_screen_flags"] == 1
        assert client.get("/api/targets/fixture-0/structure").json() is None
        response = client.get(a["structure_url"])
        assert response.content == (structure_source.parent / "target.cif").read_bytes()
        assert "immutable" in response.headers["cache-control"]
        assert client.get("/api/structures/secret.json").status_code == 404
    with pytest.raises(FileExistsError):
        publish_structures(index, structure_source, root)


def test_structure_publication_rejects_mislabeled_coordinates_and_corruption(
    index, structure_source, tmp_path
):
    source = json.loads(structure_source.read_text())
    source["structures"][0].update(
        amino_acid_sequence="ME", sequence_sha256=sha256(b"ME").hexdigest()
    )
    structure_source.write_text(json.dumps(source))
    with pytest.raises(ValueError, match="exact complete sequence"):
        publish_structures(index, structure_source, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()
    source["structures"][0].update(
        amino_acid_sequence="MQ", sequence_sha256=sha256(b"MQ").hexdigest()
    )
    structure_source.write_text(json.dumps(source))
    root = tmp_path / "good"
    publish_structures(index, structure_source, root)
    path = next((root / "files").glob("*.cif"))
    path.write_bytes(b"x" * path.stat().st_size)
    with pytest.raises(ValueError, match="checksum mismatch"):
        StructureLibrary(root, Catalog(index).metadata()["manifest_sha256"])


def test_structure_library_rejects_a_different_catalog(index, structure_source, tmp_path):
    root = tmp_path / "structures"
    publish_structures(index, structure_source, root)
    with pytest.raises(ValueError, match="does not match"):
        StructureLibrary(root, "0" * 64)


@pytest.fixture
def published(gateway, tmp_path):
    client, _, _ = gateway
    job = client.post("/api/demo", json={}, headers=HEADERS).json()
    destination = tmp_path / "published"
    export_snapshot(
        destination,
        lambda path: client.get(path).json(),
        lambda path: client.get(path).content,
        [job["job_id"]],
    )
    return destination, job["job_id"]


def test_public_snapshot_serves_exact_data_without_upstream(index, package, published, gateway):
    root, job_id = published
    _, upstream, _ = gateway
    upstream.calls.clear()
    upstream.failure = ServiceError(503, "Upstream is now offline")
    with TestClient(create_app(index, published=root, hosts=("demo.example",))) as client:
        metadata = client.get("/api/catalog", headers={"Origin": "https://other.example"})
        assert metadata.json()["read_only"] is True
        assert metadata.json()["service_configured"] is False
        assert metadata.headers["access-control-allow-origin"] == "*"
        assert "max-age=60" in metadata.headers["cache-control"]
        assert client.get("/api/targets").json()["total"] == 3
        assert client.get("/api/targets/fixture-1").json() == package[1][1]
        assert [j["job_id"] for j in client.get("/api/jobs").json()] == [job_id]
        result = client.get(f"/api/jobs/{job_id}/results").json()
        candidate = result["candidates"][0]
        assert candidate["native_filter_passed"] is False
        response = client.get(candidate["complex_structure_url"])
        assert response.content == b"artifact bytes"
        assert response.headers["x-content-type-options"] == "nosniff"
        assert client.get(result["download_url"]).content == b"artifact bytes"
        assert client.get("/api/jobs", headers={"Host": "demo.example"}).status_code == 200
        assert client.get("/api/jobs", headers={"Host": "evil.example"}).status_code == 400
        assert client.get(f"/api/jobs/{'c' * 32}").status_code == 404
        assert client.get(f"/api/jobs/{job_id}/artifacts/secret.txt").status_code == 404
        for path in (
            "/api/prepare",
            "/api/demo",
            f"/api/jobs/{job_id}/run",
            f"/api/jobs/{job_id}/cancel",
        ):
            assert client.post(path, json={}, headers=HEADERS).status_code == 403
        schema = client.get("/openapi.json").json()
        assert all(set(methods) == {"get"} for methods in schema["paths"].values())
    assert not upstream.calls


def test_public_snapshot_cannot_use_live_credentials_or_private_receipts(
    index, published, tmp_path
):
    root, _ = published
    for kwargs in ({"runs": tmp_path}, {"service": FakeService()}, {"allow_run": True}):
        with pytest.raises(ValueError, match="Public snapshots"):
            create_app(index, published=root, **kwargs)
    with pytest.raises(ValueError, match="require a verified catalog"):
        create_app(published=root)
    with pytest.raises(ValueError, match="exact public hostnames"):
        create_app(index, published=root, hosts=("*",))


def test_public_snapshot_rejects_corruption_and_catalog_mismatch(index, published):
    root, _ = published
    with pytest.raises(ValueError, match="does not match"):
        PublishedSnapshot(root, "0" * 64)
    artifact_path = next((root / "files").iterdir())
    artifact_path.write_bytes(b"x" * artifact_path.stat().st_size)
    with pytest.raises(ValueError, match="checksum mismatch"):
        PublishedSnapshot(root, Catalog(index).metadata()["manifest_sha256"])


def test_publication_is_explicit_complete_and_never_overwrites(gateway, published, tmp_path):
    client, _, _ = gateway
    root, job_id = published
    with pytest.raises(FileExistsError):
        export_snapshot(root, None, None, [job_id])
    prepared = prepare(client)
    output = tmp_path / "incomplete"
    with pytest.raises(ValueError, match="Only completed"):
        export_snapshot(output, lambda path: client.get(path).json(), None, [prepared["job_id"]])
    assert not output.exists()
    assert not list(tmp_path.glob(".published-*"))


@pytest.mark.parametrize("suffix", ["../secret", "%2e%2e/secret", "/secret", "a\\b", "a?x", "a#x"])
def test_published_artifact_paths_cannot_escape_job(suffix):
    with pytest.raises(ValueError, match="does not belong"):
        artifact_key(f"/api/jobs/{REMOTE}/artifacts/{suffix}", REMOTE)


def test_frontend_build_is_optional_for_api_and_assets_are_compressed(tmp_path, monkeypatch):
    monkeypatch.setattr(demo_service, "STATIC", tmp_path)
    with TestClient(create_app()) as client:
        assert client.get("/").status_code == 503
        assert client.get("/api/jobs").json() == []
        (tmp_path / "assets").mkdir()
        content = "/* compiled asset */" * 200
        (tmp_path / "assets/bundle.js").write_text(content)
        (tmp_path / "index.html").write_text("<main>React application</main>")
        response = client.get("/assets/bundle.js", headers={"Accept-Encoding": "gzip"})
        assert response.text == content
        assert response.headers["content-encoding"] == "gzip"
        assert "immutable" in response.headers["cache-control"]
        assert client.get("/").headers["cache-control"] == "no-store"
        assert client.get("/api/jobs").headers["cache-control"] == "no-store"


@pytest.fixture
def package(tmp_path):
    root = tmp_path / "package"
    (root / "evidence").mkdir(parents=True)
    records = []
    for i, sequence in enumerate(["ME", "MQ", "MQ"]):
        records.append(
            {
                "target_id": f"fixture-{i}",
                "donor": "synthetic",
                "gene_id": "ENSG_FIXTURE",
                "gene_symbol": "DEMO",
                "transcript_id": f"ENST_FIXTURE.{i}",
                "phase_scenario": f"synthetic:{i}",
                "amino_acid_sequence": sequence,
                "sequence_sha256": sha256(sequence.encode()).hexdigest(),
                "reference_sequence_sha256": sha256(b"ME").hexdigest(),
                "expression": None,
            }
        )
    files = {
        "protein-targets.jsonl": "".join(json.dumps(r) + "\n" for r in records),
        "proteins.fasta": "".join(
            f">sha256:{sha256(s.encode()).hexdigest()}\n{s}\n" for s in ["ME", "MQ"]
        ),
        "evidence/alleles.jsonl": "".join(
            json.dumps(
                {
                    "status": "resolved",
                    "phase_scenario": r["phase_scenario"],
                    "protein_sha256": r["sequence_sha256"],
                }
            )
            + "\n"
            for r in records
        ),
        "evidence/dispositions.jsonl": "".join(
            json.dumps(
                {"status": "selected", "transcript_id": r["transcript_id"], "released_alleles": 1}
            )
            + "\n"
            for r in records
        ),
    }
    artifacts = []
    for name, contents in files.items():
        path = root / name
        path.write_text(contents)
        artifacts.append(
            artifact(path, name, "FASTA" if name.endswith("fasta") else "JSONL").model_dump(
                mode="json"
            )
        )
    (root / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "stage1-proteins-1",
                "status": "sequence_export_complete",
                "counts": {"target_records": 3},
                "unique_sequences": 2,
                "reference_transcripts": 3,
                "artifacts": artifacts,
            }
        )
    )
    return root, records


@pytest.fixture
def index(package, tmp_path):
    destination = tmp_path / "catalog.sqlite3"
    build_index(package[0], destination)
    return destination


class FakeService:
    def __init__(self):
        self.calls = []
        self.failure = None
        self.foreign = False
        self.wrong_id = False
        self.cached = False

    def request(self, path, payload=None, *, binary=False):
        self.calls.append((path, payload))
        if self.failure:
            raise self.failure
        if binary:
            return b"artifact bytes"
        if path.endswith("/results"):
            prefix = f"/v1/designs/{'b' * 32 if self.foreign else REMOTE}/artifacts/"
            return {
                "run_id": REMOTE,
                "candidates": [
                    {
                        "id": "c1",
                        "sequence": "ME",
                        "native_filter_passed": False,
                        "binder_structure_url": prefix + "binder.cif",
                        "complex_structure_url": prefix + "complex.cif",
                    }
                ],
                "download_url": f"/v1/designs/{REMOTE}/download",
            }
        if payload is not None and path in {"/v1/designs", "/v1/designs/demo"}:
            self.cached = path.endswith("/demo")
        cached = self.cached
        return {
            "run_id": "b" * 32 if self.wrong_id else REMOTE,
            "status": "cached_demo_completed" if cached else "prepared_pending_gpu_execution",
            "execution_performed": False,
            "mode": "actual_cached_replay" if cached else "prepared",
        }


@pytest.fixture
def gateway(index, tmp_path):
    upstream = FakeService()
    app = create_app(index, tmp_path / "jobs", upstream)
    with TestClient(app) as client:
        yield client, upstream, tmp_path / "jobs"


def prepare(client):
    result = client.post("/api/prepare", json={"target_id": "fixture-1"}, headers=HEADERS)
    assert result.status_code == 201, result.text
    return result.json()


def test_catalog_preserves_biological_records_with_duplicate_sequences(index, package):
    catalog = Catalog(index)
    assert catalog.metadata()["target_records"] == 3
    assert catalog.metadata()["unique_sequences"] == 2
    assert catalog.search(changed=True)["total"] == 2
    assert catalog.search("demo", limit=1)["total"] == 3
    assert catalog.search("%")["total"] == 0
    assert catalog.search("' OR 1=1 --")["total"] == 0
    for query in ("", "d", "demo", "ENST_FIXTURE", "fixture-1", "%"):
        for changed in (False, True):
            expected = [
                r
                for r in package[1]
                if (not changed or r["sequence_sha256"] != r["reference_sequence_sha256"])
                and (
                    r["gene_symbol"].lower().startswith(query.lower())
                    or r["transcript_id"].lower().startswith(query.lower())
                    or r["target_id"] == query
                )
            ]
            pages = [catalog.search(query, changed, limit=1, offset=i) for i in range(4)]
            assert all(page["total"] == len(expected) for page in pages)
            assert [r["target_id"] for page in pages for r in page["records"]] == [
                r["target_id"] for r in expected
            ]
    assert catalog.get("fixture-1") == package[1][1]
    with pytest.raises(FileExistsError):
        build_index(package[0], index)


def test_corrupt_package_never_publishes_index(package, tmp_path):
    (package[0] / "protein-targets.jsonl").write_text("corrupt\n")
    destination = tmp_path / "bad.sqlite3"
    with pytest.raises(ValueError, match="artifact integrity"):
        build_index(package[0], destination)
    assert not destination.exists()


def test_prepare_keeps_exact_record_and_is_idempotent(gateway, package):
    client, upstream, root = gateway
    first = prepare(client)
    assert prepare(client) == first
    assert upstream.calls == [
        (
            "/v1/designs",
            {
                **package[1][1],
                "binder_length": 75,
                "num_designs": 2,
            },
        )
    ]
    receipt = json.loads((root / (first["job_id"] + ".json")).read_text())
    assert receipt["manifest_sha256"]
    assert (
        client.post(f"/api/jobs/{first['job_id']}/run", json={}, headers=HEADERS).status_code == 403
    )
    assert len(upstream.calls) == 1


def test_timeout_remains_unknown_and_never_reposts(gateway):
    client, upstream, _ = gateway
    upstream.failure = ServiceError(502, "Timed out")
    request = {"target_id": "fixture-1"}
    assert client.post("/api/prepare", json=request, headers=HEADERS).status_code == 502
    retry = client.post("/api/prepare", json=request, headers=HEADERS)
    assert retry.json()["status"] == "submission_unknown"
    assert len(upstream.calls) == 1


def test_definite_auth_failure_can_be_retried(gateway):
    client, upstream, _ = gateway
    upstream.failure = ServiceError(401, "Missing token")
    assert (
        client.post("/api/prepare", json={"target_id": "fixture-1"}, headers=HEADERS).status_code
        == 401
    )
    upstream.failure = None
    assert prepare(client)["status"] == "prepared_pending_gpu_execution"
    assert len(upstream.calls) == 2


def test_cached_replay_is_separate_and_keeps_failed_candidates(gateway):
    client, upstream, _ = gateway
    result = client.post("/api/demo", json={}, headers=HEADERS).json()
    assert result["target_id"] is None
    assert result["mode"] == "actual_cached_replay"
    assert result["execution_performed"] is False
    base = f"/api/jobs/{result['job_id']}"
    candidate = client.get(base + "/results").json()["candidates"][0]
    assert candidate["native_filter_passed"] is False
    assert candidate["binder_structure_url"] == base + "/artifacts/binder.cif"
    assert client.get(candidate["binder_structure_url"]).content == b"artifact bytes"
    assert client.get(base + "/download").headers["content-type"] == "application/zip"
    upstream.foreign = True
    assert client.get(base + "/results").status_code == 502


def test_unknown_targets_and_jobs_never_reach_upstream(gateway):
    client, upstream, _ = gateway
    assert (
        client.post("/api/prepare", json={"target_id": "missing"}, headers=HEADERS).status_code
        == 404
    )
    assert client.get(f"/api/jobs/{'c' * 32}/download").status_code == 404
    assert client.get("/api/targets?limit=1000").status_code == 422
    assert client.get(f"/api/targets?offset={2**63}").status_code == 422
    assert not upstream.calls


def test_status_from_a_different_job_is_rejected(gateway):
    client, upstream, _ = gateway
    job = prepare(client)
    upstream.wrong_id = True
    assert client.get(f"/api/jobs/{job['job_id']}").status_code == 502


def test_browser_mutations_require_same_origin_and_custom_header(gateway):
    client, upstream, _ = gateway
    assert client.post("/api/demo", json={}).status_code == 403
    assert (
        client.post(
            "/api/demo", json={}, headers={**HEADERS, "Origin": "https://untrusted.example"}
        ).status_code
        == 403
    )
    assert client.get("/api/catalog", headers={"Host": "untrusted.example"}).status_code == 400
    assert not upstream.calls


def test_unknown_dispatch_is_not_repeated_after_status_refresh(index, tmp_path):
    upstream = FakeService()
    with TestClient(create_app(index, tmp_path / "jobs", upstream, allow_run=True)) as client:
        job = prepare(client)
        base = f"/api/jobs/{job['job_id']}"
        upstream.failure = ServiceError(502, "Timeout after sending")
        assert client.post(base + "/run", json={}, headers=HEADERS).status_code == 502
        upstream.failure = None
        assert client.get(base).json()["dispatch_attempted"] is True
        assert client.post(base + "/run", json={}, headers=HEADERS).status_code == 409
        assert sum(path.endswith("/run") for path, _ in upstream.calls) == 1


def test_budget_rejection_can_be_explicitly_retried_after_capacity_is_restored(index, tmp_path):
    upstream = FakeService()
    with TestClient(create_app(index, tmp_path / "jobs", upstream, allow_run=True)) as client:
        job = prepare(client)
        base = f"/api/jobs/{job['job_id']}"
        upstream.failure = ServiceError(402, "Design service budget is exhausted")
        assert client.post(base + "/run", json={}, headers=HEADERS).status_code == 402
        saved = client.get("/api/jobs").json()[0]
        assert saved["status"] == "prepared_pending_gpu_execution"
        assert saved["dispatch_attempted"] is False
        upstream.failure = None
        response = client.post(base + "/run", json={}, headers=HEADERS)
        assert response.status_code == 200
        assert response.json()["error"] is None
        assert response.json()["dispatch_attempted"] is True


@pytest.mark.parametrize("url", ["http://remote.example", "https://a:b@host", "https://host/path"])
def test_client_rejects_unsafe_origins(url):
    with pytest.raises(ValueError):
        DesignService(url)


def test_client_sends_token_only_to_configured_origin_and_redacts_errors():
    client = DesignService("https://service.example", "private-token")

    class Opener:
        def open(self, request, timeout):
            assert request.full_url == "https://service.example/v1/designs"
            assert request.get_header("Authorization") == "Bearer private-token"
            assert request.get_header("User-agent") == "bespoke-bio/0.1"
            assert json.loads(request.data) == {"sequence": "ME"}
            raise HTTPError(request.full_url, 401, "private-token", {}, BytesIO(b"private-token"))

    client.opener = Opener()
    with pytest.raises(ServiceError) as error:
        client.request("/v1/designs", {"sequence": "ME"})
    assert error.value.status == 401
    assert "private-token" not in str(error.value)
    with pytest.raises(ValueError):
        client.request("https://elsewhere.example/v1/designs")


def test_design_catalog_keeps_distinct_jobs_and_native_failures(index, published):
    root, job_id = published
    manifest = json.loads((root / "manifest.json").read_text())
    original = manifest["jobs"][job_id]
    # Distinct jobs with the same biological label must remain distinct rows.
    for number in range(27):
        jid = f"{number:032x}"
        entry = json.loads(json.dumps(original))
        entry["job"].update(job_id=jid, target_id=f"external:{number}")
        entry["results"].update(run_id=jid, target_id=f"external:{number}")
        entry["design"] = {
            "target_id": f"external:{number}",
            "name": "FGF2",
            "amino_acid_sequence": "ME",
            "sequence_sha256": sha256(b"ME").hexdigest(),
        }
        for url in list(manifest["files"]):
            if job_id in url:
                manifest["files"][url.replace(job_id, jid)] = manifest["files"][url]
        entry["results"] = json.loads(json.dumps(entry["results"]).replace(job_id, jid))
        manifest["jobs"][jid] = entry
    (root / "manifest.json").write_text(json.dumps(manifest))
    with TestClient(create_app(index, published=root)) as client:
        a = client.get("/api/designs?q=fgf2&limit=25").json()
        b = client.get("/api/designs?q=fgf2&limit=25&offset=25").json()
        assert a["total"] == b["total"] == 27
        assert a["candidate_total"] == b["candidate_total"] == 27
        assert a["structured_candidates"] == b["structured_candidates"] == 27
        assert len(a["records"]) == 25 and len(b["records"]) == 2
        assert len({r["job_id"] for r in a["records"] + b["records"]}) == 27
        assert all(r["candidate_count"] == 1 and r["passed_count"] == 0 for r in a["records"])
        detail = client.get("/api/designs/" + a["records"][0]["job_id"]).json()
        assert detail["name"] == "FGF2" and detail["protein_structure"] is None
        assert client.get("/api/designs/" + "f" * 32).status_code == 404
    first = manifest["jobs"]["0" * 32]
    first["design"]["sequence_sha256"] = "0" * 64
    (root / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="sequence checksum"):
        PublishedSnapshot(root, manifest["catalog_sha256"])


def test_campaign_import_rejects_sequence_and_chain_role_mismatches(structure_source):
    from io import StringIO

    from Bio.PDB import MMCIFIO

    from bespoke_bio.demo.import_campaigns import complex_roles
    from bespoke_bio.demo.import_science import read_structure

    data = (structure_source.parent / "target.cif").read_bytes()
    model, _ = read_structure(data)
    key = model[0]["A"].copy()
    key.id = "Z"
    key[(" ", 2, " ")].resname = "GLU"
    model[0].add(key)
    stream = StringIO()
    writer = MMCIFIO()
    writer.set_structure(model)
    writer.save(stream)
    data = stream.getvalue().encode()
    _, roles = complex_roles(data, "MQ", "ME")
    assert roles == {"target": ["A"], "binder": ["Z"]}
    with pytest.raises(ValueError, match="sequences or chain roles"):
        complex_roles(data, "MM", "ME")
    with pytest.raises(ValueError, match="sequences or chain roles"):
        complex_roles(data, "MQ", "MQ")


def test_campaign_bundle_requires_hashes_for_every_structure(tmp_path):
    from bespoke_bio.demo.import_campaigns import verified_bundle, within

    (tmp_path / "binder.cif").write_bytes(b"synthetic coordinates")
    (tmp_path / "complex.cif").write_bytes(b"synthetic complex")
    hashes = {
        name: sha256((tmp_path / name).read_bytes()).hexdigest()
        for name in ("binder.cif", "complex.cif")
    }
    (tmp_path / "manifest.json").write_text(json.dumps(hashes))
    bundle = {
        "target_states": [{"chains": {"A": "MQ"}}],
        "components": [
            {
                "id": "binder_abc",
                "sequence": "ME",
                "sequence_sha256": sha256(b"ME").hexdigest(),
                "structure": "binder.cif",
                "evidence": {"complex_structure": "complex.cif"},
            }
        ],
    }
    (tmp_path / "bundle.json").write_text(json.dumps(bundle))
    assert verified_bundle(tmp_path)[1] == "MQ"
    bundle["components"][0]["sequence"] = "MM"
    (tmp_path / "bundle.json").write_text(json.dumps(bundle))
    with pytest.raises(ValueError, match="sequence checksum"):
        verified_bundle(tmp_path)
    (tmp_path / "complex.cif").write_bytes(b"corrupted coordinates")
    with pytest.raises(ValueError, match="checksum mismatch"):
        verified_bundle(tmp_path)
    with pytest.raises(ValueError, match="escapes"):
        within(tmp_path, "../escape")
