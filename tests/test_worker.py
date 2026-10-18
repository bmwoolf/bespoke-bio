"""Exercise orchestration without SSH/network; cached native outputs are test fixtures."""

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from biocad import worker
from biocad.boltzgen_adapter import prepare_target

ROOT = Path(__file__).resolve().parents[1]


def prepared_job(tmp_path):
    prepared = tmp_path / "prepared"
    prepare_target(ROOT / "examples/5JDS.cif", prepared, chain="A")
    return prepared


def test_worker_transfers_collects_and_exports_native_candidates(tmp_path, monkeypatch):
    prepared = prepared_job(tmp_path)
    calls = []

    def remote(command, **kwargs):
        calls.append((command, kwargs))
        if command[0] == "scp" and command[-1] == str(tmp_path / "out/downloaded"):
            results = Path(command[-1]) / "results"
            ranked = results / "final_ranked_designs"
            ranked.mkdir(parents=True)
            shutil.copy2(ROOT / "examples/boltzgen_results/all_designs_metrics.csv", ranked)
            for candidate in (ROOT / "examples/boltzgen_results/candidates").glob("*.cif"):
                shutil.copy2(candidate, ranked)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(worker.subprocess, "run", remote)
    result = worker.run_job(prepared, tmp_path / "out", "worker.example", 2222, tmp_path / "key")
    assert result["status"] == "completed"
    assert result["candidates"] == 2
    assert result["execution_performed"] is True
    assert result["pod_shutdown_performed"] is False
    assert len(calls) == 4
    assert calls[1][0][0] == calls[3][0][0] == "scp"
    script = calls[2][1]["input"]
    assert "--protocol protein-anything" in script
    assert "--cache /workspace/cache/boltzgen" in script
    assert "--num_designs 2" in script
    assert "--esmfold2_python" in script
    bundle = json.loads((tmp_path / "out/bundle/bundle.json").read_text())
    assert len(bundle["components"]) == 2
    assert len({c["id"] for c in bundle["components"]}) == 2
    assert {e["status"] for e in bundle["evaluations"] if e["id"].endswith("_filters")} == {
        "pass",
        "fail",
    }
    assert json.loads((tmp_path / "out/worker_status.json").read_text())["status"] == "completed"


@pytest.mark.parametrize("failure", ["process", "timeout"])
def test_failed_inference_still_downloads_partial_outputs(tmp_path, monkeypatch, failure):
    prepared = prepared_job(tmp_path)
    calls = []

    def remote(command, **kwargs):
        calls.append(command)
        if command[-1] == "bash -s":
            if failure == "timeout":
                raise subprocess.TimeoutExpired(command, 1)
            raise subprocess.CalledProcessError(1, command)
        if command[0] == "scp" and command[-1] == str(tmp_path / "out/downloaded"):
            destination = Path(command[-1])
            destination.mkdir()
            (destination / "partial.log").write_text("partial native output")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(worker.subprocess, "run", remote)
    monkeypatch.setattr(
        worker, "export_bundle", lambda *a, **k: pytest.fail("Must not export failed inference")
    )
    result = worker.run_job(
        prepared, tmp_path / "out", "worker.example", 2222, tmp_path / "key", timeout=1
    )
    assert result["status"] == "failed_or_interrupted"
    assert result["execution_performed"] is True
    assert result["remote_cleanup_required"] is True
    assert "error" in result
    assert len(calls) == 4
    assert (tmp_path / "out/downloaded/partial.log").read_text() == "partial native output"
    assert not (tmp_path / "out/bundle").exists()


def test_plan_path_traversal_rejected_before_remote_execution(tmp_path, monkeypatch):
    prepared = prepared_job(tmp_path)
    plan_path = prepared / "boltzgen_plan.json"
    plan = json.loads(plan_path.read_text())
    plan["command"][2] = "../outside.yaml"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setattr(
        worker.subprocess, "run", lambda *a, **k: pytest.fail("Must not execute invalid plan")
    )
    with pytest.raises(ValueError, match="local prepared basename"):
        worker.run_job(prepared, tmp_path / "out", "worker.example", 2222, tmp_path / "key")
