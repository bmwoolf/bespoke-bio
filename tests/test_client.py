import io
import json
import zipfile
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlsplit

import pytest
from fastapi.testclient import TestClient

from bespoke_bio.api import TranslateRequest, build_bundle
from biocad.client import DesignAPIError, DesignClient, bundle_payload
from biocad.service import create_app

ROOT = Path(__file__).resolve().parents[1]


def upstream_package(tmp_path):
    package = tmp_path / "upstream"
    build_bundle(
        TranslateRequest.model_validate_json((ROOT / "examples/case.json").read_text()), package
    )
    return package


def connected_client(tmp_path):
    api = TestClient(create_app(tmp_path / "runs"))

    def opener(request, timeout):
        response = api.request(
            request.get_method(),
            urlsplit(request.full_url).path,
            content=request.data,
            headers=dict(request.header_items()),
        )
        if response.status_code >= 400:
            raise HTTPError(
                request.full_url,
                response.status_code,
                response.text,
                response.headers,
                io.BytesIO(response.content),
            )
        return io.BytesIO(response.content)

    return DesignClient("http://testserver", opener=opener), api


def test_actual_upstream_package_to_independent_design_jobs(tmp_path):
    package = upstream_package(tmp_path)
    client, _ = connected_client(tmp_path)
    submitted = client.submit_bundle(package, binder_length=80, num_designs=3)
    assert submitted["prepared_targets"] == 1
    record = submitted["records"][0]
    child = client.status(record["child_run_id"])
    assert child["plan"]["binder_residues"] == 80
    assert child["plan"]["num_designs"] == 3
    assert child["plan"]["target_input_mode"] == "sequence_only"
    assert client.wait(child["run_id"], timeout=0)["execution_performed"] is False
    with zipfile.ZipFile(client.download(child["run_id"], tmp_path / "prepared.zip")) as archive:
        state = json.loads(archive.read("prepared/target_state.json"))
    assert state["context"]["upstream_record"]["source_variant_ids"] == ["demo:4:G:C"]
    assert state["measurements"][0]["value"] == 42.7
    with pytest.raises(DesignAPIError, match="503"):
        client.start(child["run_id"])
    with pytest.raises(DesignAPIError, match="409"):
        client.results(child["run_id"])


def test_results_client_reads_real_cached_candidates_and_downloads(tmp_path):
    client, api = connected_client(tmp_path)
    demo = api.post("/v1/designs/demo").json()
    result = client.results(demo["run_id"])
    assert result["candidate_count"] == 2
    assert result["candidates"][0]["native_filter_metrics"]["filter_rmsd"] > 0
    assert result["candidates"][0]["native_filter_metric_units"]["filter_rmsd"] == "angstrom"
    assert {c["native_filter_passed"] for c in result["candidates"]} == {True, False}
    with zipfile.ZipFile(client.download(demo["run_id"], tmp_path / "actual.zip")) as archive:
        assert archive.testzip() is None
        assert "bundle/bundle.json" in archive.namelist()


def test_declared_structure_embedding_and_path_containment(tmp_path):
    package = upstream_package(tmp_path)
    record = json.loads((package / "targets.jsonl").read_text().splitlines()[0])
    structure = package / "targets/one.cif"
    structure.write_text("test structure bytes")
    record["structure"] = {"status": "complete", "file": "targets/one.cif"}
    (package / "targets.jsonl").write_text(json.dumps(record) + "\n")
    assert bundle_payload(package)["artifact_texts"] == {"targets/one.cif": "test structure bytes"}
    outside = tmp_path / "outside.cif"
    outside.write_text("outside")
    structure.unlink()
    structure.symlink_to(outside)
    with pytest.raises(ValueError, match="inside"):
        bundle_payload(package)
    record["structure"]["file"] = "../outside.cif"
    (package / "targets.jsonl").write_text(json.dumps(record) + "\n")
    with pytest.raises(ValueError, match="inside"):
        bundle_payload(package)


def test_hash_mismatch_rejected_before_network(tmp_path):
    package = upstream_package(tmp_path)
    record = json.loads((package / "targets.jsonl").read_text().splitlines()[0])
    record["sequence_hash"] = "sha256:wrong"
    (package / "targets.jsonl").write_text(json.dumps(record) + "\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        bundle_payload(package)


def test_client_rejects_cross_origin_returned_paths_and_wait_timeout(tmp_path):
    client, _ = connected_client(tmp_path)
    with pytest.raises(ValueError, match="relative"):
        client._request("https://other.example/collect")
    client.status = lambda run_id: {"status": "running"}
    with pytest.raises(TimeoutError, match="not canceled"):
        client.wait("a" * 32, timeout=0)


def test_service_token_is_sent_only_to_configured_origin(monkeypatch):
    monkeypatch.setenv("BIOCAD_API_TOKEN", "local-test-token")
    calls = []

    def opener(request, timeout):
        calls.append(request)
        return io.BytesIO(b'{"status":"cancel_requested"}')

    client = DesignClient("https://service.example", opener=opener)
    assert client.cancel("a" * 32)["status"] == "cancel_requested"
    assert calls[0].get_header("Authorization") == "Bearer local-test-token"
    assert calls[0].get_method() == "POST"
    with pytest.raises(ValueError, match="relative"):
        client._request("//other.example/steal")
    assert len(calls) == 1


def test_direct_target_json_client(tmp_path):
    client, _ = connected_client(tmp_path)
    response = client.submit_target({"sequence": "ACDEFGHIKLMNPQRSTVWY", "num_designs": 1})
    assert response["status"] == "prepared_pending_gpu_execution"
    assert response["plan"]["num_designs"] == 1
    with pytest.raises(TypeError, match="object"):
        client.submit_target(["ACDE"])
