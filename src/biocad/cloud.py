"""Disposable Runpod GPUs backed by the already verified persistent model cache.

Only this module creates resources. The service owns budget reservation and serial dispatch.
Credentials stay in the gateway environment; raw Runpod responses are never persisted.
"""

import json
import math
import os
import re
import shlex
import socket
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from pathlib import Path

from .worker import run_job

VOLUME_ID = "sjltm1vrwh"
DATA_CENTER = "US-MO-1"
GPU_PRIORITY = ("NVIDIA A100-SXM4-80GB", "NVIDIA H100 80GB HBM3")
MAX_RATE = 4.0


class CloudError(RuntimeError):
    def __init__(self, message, *, code=None, pod_id=None):
        super().__init__(message)
        self.code, self.pod_id = code, pod_id


def _cli(*args, timeout=90):
    binary = os.environ.get("RUNPODCTL", "runpodctl")
    result = subprocess.run(
        [binary, *args], capture_output=True, text=True, timeout=timeout, check=False
    )
    try:
        payload = json.loads(result.stdout or result.stderr)
    except (ValueError, TypeError):
        payload = {}
    if result.returncode:
        # Do not propagate raw stdout/stderr: create/get may contain template passwords.
        raise CloudError(
            "Runpod operation failed: " + " ".join(args[:2]),
            code=payload.get("code"),
            pod_id=payload.get("id"),
        )
    return payload


def _write(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", dir=path.parent, prefix=path.name + ".", suffix=".partial", delete=False
        ) as stream:
            temporary = Path(stream.name)
            json.dump(payload, stream, indent=2)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def cleanup_owned_pod(pod_id):
    """Delete an explicitly owned lease and confirm 404; never stop unrelated resources."""
    for _ in range(3):
        try:
            _cli("pod", "delete", pod_id)
        except (CloudError, OSError, subprocess.SubprocessError):
            pass
        try:
            _cli("pod", "get", pod_id)
        except CloudError as exc:
            if exc.code == "not_found":
                return True
        except (OSError, subprocess.SubprocessError):
            pass
        time.sleep(1)
    return False


def _find_created_pod(name):
    pods = _cli("pod", "list")
    return next((p["id"] for p in pods if p.get("name") == name), None)


def reconcile_owned_leases(state_root=None):
    """Gateway restart recovery, including create requests interrupted before an ID returned."""
    root = Path(state_root or os.environ.get("BIOCAD_CLOUD_STATE_DIR", "cloud_leases"))
    remaining = []
    for path in root.glob("*.json"):
        lease = json.loads(path.read_text())
        if lease.get("state") == "deleted":
            continue
        try:
            pod_id = lease.get("pod_id") or _find_created_pod(lease["name"])
            if pod_id and not cleanup_owned_pod(pod_id):
                remaining.append(pod_id)
                continue
            lease.update(state="deleted", deleted_at=time.time())
            _write(path, lease)
        except (CloudError, OSError, subprocess.SubprocessError):
            remaining.append(lease.get("pod_id") or lease.get("name"))
    return remaining


def _watchdog(lease_path):
    """Independent process survives an interrupted HTTP job or dispatcher thread."""
    path = Path(lease_path)
    while path.exists():
        lease = json.loads(path.read_text())
        if lease.get("state") == "deleted":
            return
        if time.time() >= lease["deadline"]:
            try:
                pod_id = lease.get("pod_id") or _find_created_pod(lease["name"])
                deleted = not pod_id or cleanup_owned_pod(pod_id)
            except (CloudError, OSError, subprocess.SubprocessError):
                deleted = False
            if deleted:
                lease.update(state="deleted", deleted_at=time.time(), watchdog_deleted=True)
                _write(path, lease)
                return
            time.sleep(5)
        else:
            time.sleep(min(5, lease["deadline"] - time.time()))


def _start_watchdog(lease_path):
    with Path(lease_path).with_suffix(".watchdog.log").open("a") as log:
        return subprocess.Popen(
            [sys.executable, "-m", "biocad.cloud", "watchdog", str(lease_path)],
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=log,
            start_new_session=True,
        )


def _ssh_base(host, port, key):
    if not host or host.startswith("-") or any(c.isspace() for c in host):
        raise CloudError("Invalid worker address")
    return [
        "ssh",
        "-i",
        str(key),
        "-p",
        str(int(port)),
        "-o",
        "BatchMode=yes",
        "-o",
        "ConnectTimeout=10",
        "-o",
        "StrictHostKeyChecking=accept-new",
        "root@" + host,
    ]


def _connect(pod_id, deadline, cancel_event):
    while time.time() < deadline:
        if cancel_event.is_set():
            raise CloudError("Job cancelled")
        try:
            address = _cli("ssh", "info", pod_id)
            host, port = address.get("ip"), address.get("port")
            if host and port:
                with socket.create_connection((host, int(port)), timeout=3):
                    return host, int(port)
        except (CloudError, OSError, subprocess.SubprocessError):
            pass
        time.sleep(2)
    raise CloudError("GPU startup deadline exceeded")


def _bootstrap(host, port, key, job_id, deadline, cancel_event, log_path):
    remote = "/workspace/cloud_bootstraps/" + uuid.uuid4().hex
    script = "\n".join(
        [
            "set -euo pipefail",
            "mkdir -p " + shlex.quote(remote),
            "cp /workspace/boltzgen_smoke/{bootstrap.sh,run.sh,pdl1.yaml,target.cif} "
            + shlex.quote(remote)
            + "/",
            "cd " + shlex.quote(remote),
            "bash bootstrap.sh warmup",
            # worker.run_job reads this stable pointer. All snapshots restore the same interpreter.
            "cp esm_python.txt /workspace/boltzgen_smoke/esm_python.txt",
        ]
    )
    with Path(log_path).open("w") as log:
        process = subprocess.Popen(
            _ssh_base(host, port, key) + ["bash -s"],
            stdin=subprocess.PIPE,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
        try:
            process.stdin.write(script)
            process.stdin.close()
            while process.poll() is None:
                if cancel_event.is_set() or time.time() >= deadline:
                    process.kill()
                    raise CloudError("Bootstrap cancelled or deadline exceeded")
                time.sleep(1)
            if process.returncode:
                raise CloudError("GPU cache bootstrap failed; inspect warmup log")
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()


def run_managed_job(prepared, output, *, job_id, cancel_event, lease_callback, timeout=7200):
    """Create, warm, execute, download and finally delete one bounded GPU lease.

    Caller must reserve $10 before invocation (maximum 2h at $4/h plus buffer).
    ``lease_callback(pod_id, state='active'|'deleted')`` persists service ownership.
    ``actual_cost_usd`` is a conservative estimate, not finalized provider billing.
    """
    if not 1 <= timeout <= 7200:
        raise ValueError("GPU lease timeout must be between 1 and 7200 seconds")
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    key = Path(os.environ.get("BIOCAD_SSH_KEY", "")).expanduser()
    if not key.is_file():
        raise ValueError("BIOCAD_SSH_KEY must name the registered gateway SSH private key")
    state_root = Path(os.environ.get("BIOCAD_CLOUD_STATE_DIR", str(output.parent / "cloud_leases")))
    lease_path = state_root / (uuid.uuid4().hex + ".json")
    started, deadline = time.time(), time.time() + timeout
    receipt = {
        "status": "failed_or_interrupted",
        "execution_performed": False,
        "job_id": job_id,
        "actual_cost_usd": 0.0,
        "billing_finalized": False,
    }
    pod_id, rate, lease, thread = None, 0.0, None, None
    try:
        if cancel_event.is_set():
            receipt["status"] = "cancelled"
            return receipt
        catalog = _cli("gpu", "list", "--include-unavailable")
        for gpu_id in GPU_PRIORITY:
            match = next((g for g in catalog if g.get("gpuId") == gpu_id), None)
            price = match.get("securePricePerHr") if match else None
            if not isinstance(price, (int, float)) or not 0 < price <= MAX_RATE:
                continue
            if not match.get("available"):
                continue
            gateway = os.environ.get("BIOCAD_GATEWAY_POD_ID", "local")
            if not re.fullmatch(r"[A-Za-z0-9-]{1,40}", gateway):
                raise CloudError("Invalid gateway pod identifier")
            name = "biocad-job-" + gateway + "-" + uuid.uuid4().hex[:12]
            lease = {
                "name": name,
                "state": "provisioning",
                "deadline": deadline,
                "started": started,
                "rate": float(price),
                "job_id": job_id,
            }
            _write(lease_path, lease)
            if not receipt.get("watchdog_armed"):
                _start_watchdog(lease_path)
                receipt["watchdog_armed"] = True
            try:
                pod = _cli(
                    "pod",
                    "create",
                    "--name",
                    name,
                    "--template-id",
                    "runpod-torch-v280",
                    "--gpu-id",
                    gpu_id,
                    "--data-center-ids",
                    DATA_CENTER,
                    "--network-volume-id",
                    VOLUME_ID,
                    "--volume-mount-path",
                    "/workspace",
                    "--min-cuda-version",
                    "13.0",
                    "--container-disk-in-gb",
                    "100",
                    "--ports",
                    "22/tcp",
                    "--ssh",
                    timeout=75,
                )
                pod_id = pod["id"]
                rate = float(pod.get("costPerHr", price))
                break
            except CloudError as exc:
                if exc.pod_id:
                    pod_id, rate = exc.pod_id, float(price)
                    break
                if exc.code != "bad_request":
                    raise
        if not pod_id:
            raise CloudError("No eligible A100/H100 capacity at or below $4/hour")
        lease.update(pod_id=pod_id, state="active", rate=rate)
        _write(lease_path, lease)
        lease_callback(pod_id, state="active")
        receipt.update(pod_id=pod_id, hourly_rate_usd=rate)
        if rate > MAX_RATE:
            raise CloudError("Allocated GPU exceeds maximum hourly rate")
        host, port = _connect(pod_id, min(deadline, time.time() + 300), cancel_event)
        _bootstrap(
            host,
            port,
            key,
            job_id,
            deadline,
            cancel_event,
            output.parent / (output.name + ".warmup.log"),
        )
        result = {}

        def execute():
            try:
                result.update(
                    run_job(
                        prepared,
                        output,
                        host,
                        port,
                        key,
                        timeout=max(1, int(deadline - time.time())),
                    )
                )
            except (RuntimeError, ValueError, KeyError, OSError, subprocess.SubprocessError):
                result.update(status="failed_or_interrupted", error="Worker execution failed")

        receipt["execution_performed"] = True
        thread = threading.Thread(target=execute, daemon=True)
        thread.start()
        while thread.is_alive():
            if cancel_event.is_set() or time.time() >= deadline:
                break
            thread.join(timeout=1)
        receipt.update(result)
        if receipt["status"] not in {"completed", "cancelled"}:
            receipt["worker_status"] = receipt["status"]
            receipt["status"] = "failed_or_interrupted"
        if cancel_event.is_set():
            receipt["status"] = "cancelled"
        elif thread.is_alive():
            receipt.update(status="failed_or_interrupted", error="GPU job deadline exceeded")
    except (CloudError, OSError, subprocess.SubprocessError, ValueError, KeyError):
        receipt.update(
            status="cancelled" if cancel_event.is_set() else "failed_or_interrupted",
            error="Cloud execution failed; inspect gateway and warmup logs",
        )
    finally:
        if not pod_id and lease:
            try:
                pod_id = _find_created_pod(lease["name"])
                if pod_id:
                    rate = lease["rate"]
                    lease_callback(pod_id, state="active")
                    receipt["pod_id"] = pod_id
                else:
                    lease.update(state="deleted", deleted_at=time.time())
                    _write(lease_path, lease)
            except (CloudError, OSError, subprocess.SubprocessError):
                receipt["cleanup_required"] = True
        if pod_id:
            deleted = cleanup_owned_pod(pod_id)
            receipt["pod_shutdown_performed"] = deleted
            receipt["cleanup_required"] = not deleted
            if deleted and thread is not None and thread.is_alive():
                thread.join(timeout=30)
            if thread is not None:
                receipt["download_thread_finished"] = not thread.is_alive()
            receipt["actual_cost_usd"] = (
                math.ceil((time.time() - started) / 3600 * rate * 100) / 100
            )
            receipt["actual_cost_usd"] += 0.05  # conservative local disk/storage allowance
            if deleted:
                if lease:
                    lease.update(state="deleted", deleted_at=time.time())
                    _write(lease_path, lease)
                lease_callback(pod_id, state="deleted")
            else:
                receipt["status"] = "failed_or_interrupted"
        output.mkdir(parents=True, exist_ok=True)
        _write(output / "cloud_status.json", receipt)
    return receipt


if __name__ == "__main__" and len(sys.argv) == 3 and sys.argv[1] == "watchdog":
    _watchdog(sys.argv[2])
