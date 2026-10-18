from pathlib import Path

from fastapi.testclient import TestClient

from biocad.service import create_app
from biocad.structures import sha

ROOT = Path(__file__).resolve().parents[1]


def client(tmp_path):
    return TestClient(create_app(tmp_path))


def test_sequence_preparation_and_safe_download(tmp_path):
    api = client(tmp_path)
    response = api.post("/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY"})
    assert response.status_code == 201
    body = response.json()
    assert body["status"] == "prepared_pending_gpu_execution"
    assert body["cloud_launched"] is False
    assert body["structural_validation"] == "pending_model_execution"
    status = api.get(body["status_url"]).json()
    assert "prepared/design.yaml" in status["artifacts"]
    artifact = api.get(body["artifact_base_url"] + "prepared/design.yaml")
    assert artifact.status_code == 200
    assert artifact.json()["entities"][1]["protein"]["sequence"] == "ACDEFGHIKLMNPQRSTVWY"
    (tmp_path / body["run_id"] / "escape").symlink_to("/etc/passwd")
    assert api.get(body["artifact_base_url"] + "escape").status_code == 404
    assert api.get("/v1/designs/bad").status_code == 404
    assert api.post("/v1/designs", json={"sequence": "BAD!"}).status_code == 422
    assert (
        api.post("/v1/designs", json={"sequence": "ACD", "structure_text": "foo"}).status_code
        == 422
    )


def test_structure_demo_and_upstream_routes_coexist(tmp_path):
    api = client(tmp_path)
    text = (ROOT / "examples/5JDS.cif").read_text()
    prepared = api.post("/v1/designs", json={"structure_text": text, "chain": "A"})
    assert prepared.status_code == 201
    assert prepared.json()["plan"]["target_residues"] == 115
    demo = api.post("/v1/designs/demo").json()
    assert demo["status"] == "cached_demo_completed"
    assert demo["execution_performed"] is False
    assert demo["candidates"] == 2 and demo["native_filter_passed"] == 1
    bundle = api.get(demo["artifact_base_url"] + demo["bundle"]).json()
    assert len(bundle["components"]) == 2
    assert "/v1/translate" in api.get("/openapi.json").json()["paths"]


def upstream_payload():
    seq = "ACDEFGHIKLMNPQRSTVWY"
    return {
        "case_id": "teammate",
        "reference": {"genome_build": "synthetic", "annotation": "test"},
        "targets": [
            {
                "target_id": "protein1",
                "gene_id": "g1",
                "gene_symbol": "G1",
                "transcript_id": "t1",
                "haplotype": "reference",
                "source_variant_ids": [],
                "expression": {"tpm": 5.0},
                "reference_amino_acid_sequence": seq,
                "amino_acid_sequence": seq,
                "amino_acid_length": len(seq),
                "sequence_hash": "sha256:" + sha(seq.encode()),
                "first_changed_residue": None,
                "translation_status": "complete",
                "warnings": [],
                "structure": {"status": "unresolved"},
            }
        ],
        "target_count": 1,
        "unique_sequence_count": 1,
    }


def test_upstream_unresolved_sequence_fallback_and_bad_hash(tmp_path):
    api = client(tmp_path)
    payload = upstream_payload()
    response = api.post("/v1/designs/upstream", json=payload)
    assert response.status_code == 201
    body = response.json()
    assert body["prepared_targets"] == 1
    assert body["records"][0]["target_input_mode"] == "sequence_only"
    assert body["execution_performed"] is False
    payload["targets"][0]["sequence_hash"] = "sha256:bad"
    body = api.post("/v1/designs/upstream", json=payload).json()
    assert body["prepared_targets"] == 0
    assert body["records"][0]["status"] == "failed"
    payload = upstream_payload()
    payload["artifact_texts"] = {"../outside.cif": "unsafe"}
    assert api.post("/v1/designs/upstream", json=payload).status_code == 422


def test_upstream_embedded_structure_preserves_import_provenance(tmp_path):
    import gemmi

    payload = upstream_payload()
    text = (ROOT / "examples/5JDS.cif").read_text()
    structure = gemmi.read_structure(str(ROOT / "examples/5JDS.cif"))
    structure.setup_entities()
    seq = gemmi.one_letter_code(next(e for e in structure.entities if e.name == "2").full_sequence)
    protein = payload["targets"][0]
    protein.update(
        amino_acid_sequence=seq,
        reference_amino_acid_sequence=seq,
        amino_acid_length=len(seq),
        sequence_hash="sha256:" + sha(seq.encode()),
        structure={"status": "complete", "source": "RCSB PDB", "file": "targets/target.cif"},
    )
    payload["artifact_texts"] = {"targets/target.cif": text}
    payload["chain_selections"] = {"protein1": ["A"]}
    api = client(tmp_path)
    response = api.post("/v1/designs/upstream", json=payload)
    assert response.status_code == 201
    body = response.json()
    assert body["prepared_targets"] == 1
    assert body["records"][0]["target_input_mode"] == "observed_structure"
    imported = api.get(body["artifact_base_url"] + body["imported_targets"]).json()
    assert imported["target_states"][0]["measurements"][0]["value"] == 5.0
    assert imported["target_states"][0]["context"]["exact_protein_sequence"] == seq

    prepared = api.get(body["artifact_base_url"] + body["records"][0]["target_state_path"]).json()
    assert prepared["context"]["upstream_target_state_id"] == imported["target_states"][0]["id"]
    assert prepared["measurements"][0]["value"] == 5.0
    assert (
        prepared["context"]["upstream_imported_state"]["context"]["exact_protein_sequence"] == seq
    )


def configure_worker(monkeypatch, tmp_path):
    key = tmp_path / "worker_key"
    key.write_text("mock key; never used for a remote call")
    monkeypatch.setenv("BIOCAD_WORKER_HOST", "mock.invalid")
    monkeypatch.setenv("BIOCAD_WORKER_PORT", "22022")
    monkeypatch.setenv("BIOCAD_WORKER_KEY", str(key))


def test_live_run_disabled_without_operator_configuration(tmp_path, monkeypatch):
    for name in ("BIOCAD_WORKER_HOST", "BIOCAD_WORKER_PORT", "BIOCAD_WORKER_KEY"):
        monkeypatch.delenv(name, raising=False)
    api = client(tmp_path)
    job = api.post("/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY"}).json()
    assert api.post(job["status_url"] + "/run").status_code == 503
    assert api.get(job["status_url"]).json()["status"] == "prepared_pending_gpu_execution"


def test_live_run_success_status_artifacts_and_duplicate_guard(tmp_path, monkeypatch):
    configure_worker(monkeypatch, tmp_path)
    calls = []

    def fake_runner(prepared, output, **kwargs):
        calls.append((prepared, kwargs))
        status = __import__("json").loads((output.parent / "status.json").read_text())
        assert status["status"] == "running"
        assert status["execution_performed"] is False
        assert status["execution_attempted"] is True
        (output / "bundle").mkdir(parents=True)
        (output / "bundle/bundle.json").write_text('{"test_fixture": true}')
        return {"status": "completed", "bundle": "bundle/bundle.json"}

    api = TestClient(create_app(tmp_path / "runs", runner=fake_runner))
    job = api.post("/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY"}).json()
    response = api.post(job["status_url"] + "/run", json={})
    assert response.status_code == 202 and response.json()["status"] == "queued"
    done = api.get(job["status_url"]).json()
    assert done["status"] == "completed"
    assert done["mode"] == "live_model_outputs" and done["cloud_launched"] is False
    assert done["bundle"] == "live/bundle/bundle.json"
    assert api.get(job["artifact_base_url"] + done["bundle"]).json()["test_fixture"] is True
    assert len(calls) == 1 and calls[0][1]["timeout"] == 3600
    assert api.post(job["status_url"] + "/run").status_code == 409


def test_live_run_failure_and_bounded_design_count(tmp_path, monkeypatch):
    configure_worker(monkeypatch, tmp_path)

    def fake_runner(*args, **kwargs):
        return {
            "status": "failed_or_interrupted",
            "error": "mock GPU failure",
            "execution_performed": True,
            "remote_cleanup_required": True,
        }

    api = TestClient(create_app(tmp_path / "runs", runner=fake_runner))
    job = api.post("/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY"}).json()
    assert api.post(job["status_url"] + "/run").status_code == 202
    failed = api.get(job["status_url"]).json()
    assert failed["status"] == "failed" and failed["execution_performed"] is True
    assert failed["worker_receipt"]["error"] == "mock GPU failure"
    large = api.post(
        "/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY", "num_designs": 33}
    ).json()
    assert api.post(large["status_url"] + "/run").status_code == 422


def test_live_upstream_target_selection_and_runner_exception(tmp_path, monkeypatch):
    configure_worker(monkeypatch, tmp_path)

    def fake_runner(*args, **kwargs):
        raise RuntimeError("mock exporter error")

    api = TestClient(create_app(tmp_path / "runs", runner=fake_runner))
    job = api.post("/v1/designs/upstream", json=upstream_payload()).json()
    assert api.post(job["status_url"] + "/run", json={"target_index": 7}).status_code == 422
    child_url = job["records"][0]["child_status_url"]
    assert api.post(job["records"][0]["child_run_url"], json={}).status_code == 202
    failed = api.get(child_url).json()
    assert failed["status"] == "failed" and failed["error"] == "mock exporter error"


def test_download_all_contains_only_run_artifacts_and_no_escaping_symlink(tmp_path):
    import io
    import zipfile

    api = client(tmp_path)
    demo = api.post("/v1/designs/demo").json()
    (tmp_path / "unrelated.txt").write_text("outside this run")
    (tmp_path / demo["run_id"] / "escape").symlink_to(tmp_path / "unrelated.txt")
    response = api.get(demo["status_url"] + "/download")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/zip"
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        assert "bundle/bundle.json" in archive.namelist()
        assert "prepared/target_state.json" in archive.namelist()
        assert "escape" not in archive.namelist()
        assert "unrelated.txt" not in archive.namelist()
        assert not any(name.endswith(".zip") for name in archive.namelist())
        assert len(__import__("json").loads(archive.read("bundle/bundle.json"))["components"]) == 2
    assert not list((tmp_path / demo["run_id"]).glob("*.zip"))


def test_upstream_batch_children_run_independently(tmp_path, monkeypatch):
    configure_worker(monkeypatch, tmp_path)
    calls = []

    def fake_runner(prepared, output, **kwargs):
        calls.append(prepared)
        (output / "bundle").mkdir(parents=True)
        (output / "bundle/bundle.json").write_text('{"mock_output":true}')
        return {"status": "completed", "bundle": "bundle/bundle.json"}

    api = TestClient(create_app(tmp_path / "runs", runner=fake_runner))
    payload = upstream_payload()
    other = {**payload["targets"][0], "target_id": "protein2", "transcript_id": "t2"}
    payload["targets"].append(other)
    payload["target_count"] = 2
    job = api.post("/v1/designs/upstream", json=payload).json()
    first, second = job["records"]
    assert first["child_run_id"] != second["child_run_id"]
    assert api.post(job["status_url"] + "/run").status_code == 422
    assert api.get(first["child_status_url"]).json()["status"] == "prepared_pending_gpu_execution"
    assert api.get(second["child_status_url"]).json()["status"] == "prepared_pending_gpu_execution"
    assert api.post(first["child_run_url"]).json()["status"] == "queued"
    assert api.get(first["child_status_url"]).json()["status"] == "completed"
    assert api.get(second["child_status_url"]).json()["status"] == "prepared_pending_gpu_execution"
    assert api.post(second["child_run_url"]).json()["status"] == "queued"
    assert len(calls) == 2 and calls[0] != calls[1]
    done = api.get(job["status_url"]).json()
    assert done["status"] == "completed"
    assert all(r["status"] == "completed" for r in done["records"])
    assert api.get(first["child_download_url"]).status_code == 200
    assert api.get(second["child_download_url"]).status_code == 200


def test_openapi_describes_job_envelopes_and_target_inputs(tmp_path):
    api = client(tmp_path)
    schema = api.get("/openapi.json").json()
    for route, method, code in [
        ("/v1/designs", "post", "201"),
        ("/v1/designs/demo", "post", "201"),
        ("/v1/designs/upstream", "post", "201"),
        ("/v1/designs/{run_id}/run", "post", "202"),
        ("/v1/designs/{run_id}", "get", "200"),
    ]:
        response = schema["paths"][route][method]["responses"][code]
        assert response["content"]["application/json"]["schema"]["$ref"].endswith("/JobResponse")
    models = schema["components"]["schemas"]
    assert models["JobResponse"]["additionalProperties"] is True
    assert models["JobResponse"]["properties"]["bundle"]["description"]
    assert models["JobResponse"]["properties"]["records"]["description"]
    assert models["JobResponse"]["examples"][0]["status"] == "prepared_pending_gpu_execution"
    for field in ("sequence", "structure_text", "chain", "binder_length", "num_designs"):
        assert models["DesignInput"]["properties"][field]["description"]
    job = api.post("/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY"}).json()
    assert job["structural_validation"] == "pending_model_execution"
    assert job["target_state_path"] == "prepared/target_state.json"


def test_connection_failure_is_an_attempt_not_model_execution(tmp_path, monkeypatch):
    configure_worker(monkeypatch, tmp_path)

    def failed_before_model(*args, **kwargs):
        return {
            "status": "failed_or_interrupted",
            "execution_performed": False,
            "error": "mock SSH connection failure",
            "remote_cleanup_required": False,
        }

    api = TestClient(create_app(tmp_path / "runs", runner=failed_before_model))
    job = api.post("/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY"}).json()
    assert api.post(job["status_url"] + "/run").status_code == 202
    status = api.get(job["status_url"]).json()
    assert status["status"] == "failed"
    assert status["execution_attempted"] is True
    assert status["execution_performed"] is False


def test_status_reports_last_native_log_stage_without_implying_completion(tmp_path):
    api = client(tmp_path)
    job = api.post("/v1/designs", json={"sequence": "ACDEFGHIKLMNPQRSTVWY"}).json()
    live = tmp_path / job["run_id"] / "live"
    live.mkdir()
    (live / "worker.log").write_text(
        "Pipeline step 1 of 7: design\nnoise\nPipeline step 2 of 7: inverse_folding\n"
    )
    status = api.get(job["status_url"]).json()
    assert status["progress"] == {"stage_index": 2, "total_stages": 7, "stage": "inverse_folding"}
    assert status["status"] == "prepared_pending_gpu_execution"
    assert status["execution_performed"] is False
    (live / "worker.log").write_text("Pipeline step 7 of 7: filtering\n")
    status = api.get(job["status_url"]).json()
    assert status["progress"]["stage_index"] == 7
    assert status["status"] == "prepared_pending_gpu_execution"
    # An old stage outside the bounded tail must not be reported as current progress.
    (live / "worker.log").write_text("Pipeline step 1 of 7: design\n" + "x" * (1024 * 1024 + 1))
    assert "progress" not in api.get(job["status_url"]).json()
