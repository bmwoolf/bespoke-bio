"""Managed GPU safety and lifecycle tests, using a fake Runpod control plane."""

import json
import threading

import pytest

from biocad import cloud


@pytest.fixture
def harness(tmp_path, monkeypatch):
    key = tmp_path / "key"
    key.write_text("fake test key")
    monkeypatch.setenv("BIOCAD_SSH_KEY", str(key))
    monkeypatch.setenv("BIOCAD_CLOUD_STATE_DIR", str(tmp_path / "leases"))
    calls = []

    def cli(*args, **kwargs):
        calls.append(args)
        if args[:2] == ("gpu", "list"):
            return [
                {"gpuId": g, "available": True, "securePricePerHr": price}
                for g, price in zip(cloud.GPU_PRIORITY, [1.79, 3.99])
            ]
        if args[:2] == ("pod", "create"):
            return {"id": "owned-pod", "costPerHr": 1.79, "env": {"secret": "never save"}}
        raise AssertionError(args)

    monkeypatch.setattr(cloud, "_cli", cli)
    monkeypatch.setattr(cloud, "_start_watchdog", lambda path: None)
    monkeypatch.setattr(cloud, "_connect", lambda *args: ("127.0.0.1", 22))
    monkeypatch.setattr(cloud, "_bootstrap", lambda *args: None)
    deleted = []
    monkeypatch.setattr(cloud, "cleanup_owned_pod", lambda pod: deleted.append(pod) or True)
    monkeypatch.setattr(
        cloud,
        "run_job",
        lambda *args, **kwargs: {
            "status": "completed",
            "execution_performed": True,
            "candidates": 1,
        },
    )
    events = []
    return tmp_path, calls, deleted, lambda pod, **kw: events.append((pod, kw["state"])), events


def test_success_download_then_delete_and_sanitized_lease(harness):
    root, calls, deleted, callback, events = harness
    receipt = cloud.run_managed_job(
        root / "input",
        root / "out",
        job_id="job-1",
        cancel_event=threading.Event(),
        lease_callback=callback,
    )
    assert receipt["status"] == "completed"
    assert deleted == ["owned-pod"]
    assert events == [("owned-pod", "active"), ("owned-pod", "deleted")]
    assert receipt["pod_shutdown_performed"] is True
    assert 0 < receipt["actual_cost_usd"] < 1
    create = next(c for c in calls if c[:2] == ("pod", "create"))
    assert create[create.index("--min-cuda-version") + 1] == "13.0"
    assert create[create.index("--network-volume-id") + 1] == cloud.VOLUME_ID
    assert all("secret" not in p.read_text() for p in root.rglob("*.json"))


def test_bootstrap_failure_still_deletes(harness, monkeypatch):
    root, _, deleted, callback, _ = harness

    def fail(*args):
        raise cloud.CloudError("bootstrap failed")

    monkeypatch.setattr(cloud, "_bootstrap", fail)
    result = cloud.run_managed_job(
        root / "input",
        root / "out",
        job_id="job",
        cancel_event=threading.Event(),
        lease_callback=callback,
    )
    assert result["status"] == "failed_or_interrupted"
    assert deleted == ["owned-pod"]


def test_capacity_falls_back_and_rate_above_cap_never_provisions(harness, monkeypatch):
    root, calls, deleted, callback, _ = harness
    original = cloud._cli

    def cli(*args, **kwargs):
        if args[:2] == ("pod", "create") and cloud.GPU_PRIORITY[0] in args:
            raise cloud.CloudError("capacity", code="bad_request")
        return original(*args, **kwargs)

    monkeypatch.setattr(cloud, "_cli", cli)
    receipt = cloud.run_managed_job(
        root / "input",
        root / "out",
        job_id="job",
        cancel_event=threading.Event(),
        lease_callback=callback,
    )
    assert receipt["status"] == "completed"
    create = next(c for c in calls if c[:2] == ("pod", "create"))
    assert cloud.GPU_PRIORITY[1] in create
    assert deleted


def test_pre_cancelled_never_creates(harness):
    root, calls, deleted, callback, _ = harness
    cancel = threading.Event()
    cancel.set()
    result = cloud.run_managed_job(
        root / "input", root / "out", job_id="job", cancel_event=cancel, lease_callback=callback
    )
    assert result["status"] == "cancelled"
    assert not calls and not deleted


def test_deleted_callback_only_after_verified_cleanup(harness, monkeypatch):
    root, _, _, callback, events = harness
    monkeypatch.setattr(cloud, "cleanup_owned_pod", lambda pod: False)
    result = cloud.run_managed_job(
        root / "input",
        root / "out",
        job_id="job",
        cancel_event=threading.Event(),
        lease_callback=callback,
    )
    assert result["cleanup_required"] is True
    assert result["status"] == "failed_or_interrupted"
    assert events == [("owned-pod", "active")]
    assert json.loads(next((root / "leases").glob("*.json")).read_text())["state"] == "active"


def test_cleanup_requires_not_found(monkeypatch):
    def cli(*args):
        if args[:2] == ("pod", "get"):
            raise cloud.CloudError("not found", code="not_found")
        return {"deleted": True}

    monkeypatch.setattr(cloud, "_cli", cli)
    assert cloud.cleanup_owned_pod("explicitly-owned") is True


def test_watchdog_deletes_expired_owned_lease(tmp_path, monkeypatch):
    lease = tmp_path / "lease.json"
    lease.write_text(json.dumps({"pod_id": "owned", "state": "active", "deadline": 0}))
    called = []
    monkeypatch.setattr(cloud, "cleanup_owned_pod", lambda pod: called.append(pod) or True)
    cloud._watchdog(lease)
    assert called == ["owned"]
    assert json.loads(lease.read_text())["watchdog_deleted"] is True


def test_timeout_rejected_before_provision(harness):
    root, calls, _, callback, _ = harness
    with pytest.raises(ValueError):
        cloud.run_managed_job(
            root / "input",
            root / "out",
            job_id="job",
            cancel_event=threading.Event(),
            lease_callback=callback,
            timeout=7201,
        )
    assert not calls


def test_catalog_prices_above_cap_never_create(harness, monkeypatch):
    root, calls, deleted, callback, _ = harness
    monkeypatch.setattr(
        cloud,
        "_cli",
        lambda *args, **kw: [
            {"gpuId": g, "available": True, "securePricePerHr": 4.01} for g in cloud.GPU_PRIORITY
        ],
    )
    result = cloud.run_managed_job(
        root / "input",
        root / "out",
        job_id="job",
        cancel_event=threading.Event(),
        lease_callback=callback,
    )
    assert result["status"] == "failed_or_interrupted"
    assert not calls and not deleted


def test_mid_execution_cancel_deletes_gpu(harness, monkeypatch):
    root, _, deleted, callback, _ = harness
    cancel = threading.Event()

    def worker(*args, **kwargs):
        cancel.set()
        return {"status": "failed_or_interrupted", "execution_performed": True}

    monkeypatch.setattr(cloud, "run_job", worker)
    result = cloud.run_managed_job(
        root / "input", root / "out", job_id="job", cancel_event=cancel, lease_callback=callback
    )
    assert result["status"] == "cancelled"
    assert deleted == ["owned-pod"]


def test_restart_reconciles_interrupted_create_by_exact_owned_name(tmp_path, monkeypatch):
    lease = tmp_path / "pending.json"
    lease.write_text(json.dumps({"name": "biocad-job-unique", "state": "provisioning"}))
    monkeypatch.setattr(
        cloud,
        "_cli",
        lambda *args: [
            {"id": "unrelated", "name": "other"},
            {"id": "created", "name": "biocad-job-unique"},
        ],
    )
    deleted = []
    monkeypatch.setattr(cloud, "cleanup_owned_pod", lambda pod: deleted.append(pod) or True)
    assert cloud.reconcile_owned_leases(tmp_path) == []
    assert deleted == ["created"]
    assert json.loads(lease.read_text())["state"] == "deleted"


def test_cli_never_exposes_raw_template_secrets(monkeypatch):
    import subprocess

    monkeypatch.setattr(
        cloud.subprocess,
        "run",
        lambda *args, **kw: subprocess.CompletedProcess(
            args,
            1,
            stdout=json.dumps({"env": {"JUPYTER_PASSWORD": "very-secret"}, "code": "bad_request"}),
            stderr="",
        ),
    )
    with pytest.raises(cloud.CloudError) as caught:
        cloud._cli("pod", "create")
    assert "very-secret" not in str(caught.value)
    assert caught.value.code == "bad_request"


def test_create_returns_id_before_wait_and_names_gateway_owner(harness, monkeypatch):
    root, calls, _, callback, _ = harness
    monkeypatch.setenv("BIOCAD_GATEWAY_POD_ID", "gateway123")
    cloud.run_managed_job(
        root / "in",
        root / "out",
        job_id="job",
        cancel_event=threading.Event(),
        lease_callback=callback,
    )
    create = next(c for c in calls if c[:2] == ("pod", "create"))
    assert "--wait" not in create
    assert "--wait-timeout" not in create
    assert create[create.index("--name") + 1].startswith("biocad-job-gateway123-")
    lease = json.loads(next((root / "leases").glob("*.json")).read_text())
    assert lease["name"].startswith("biocad-job-gateway123-")


def test_atomic_receipt_writers_use_independent_temporary_files(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    path = tmp_path / "lease.json"
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(
            pool.map(
                lambda index: cloud._write(path, {"writer": index, "text": "x" * 10000}), range(32)
            )
        )
    result = json.loads(path.read_text())
    assert 0 <= result["writer"] < 32
    assert result["text"] == "x" * 10000
    assert not list(tmp_path.glob("*.partial"))
