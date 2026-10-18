"""Execute one prepared job on an existing warmed SSH worker; never provision compute."""

import json
import shlex
import subprocess
import uuid
from pathlib import Path

from .boltzgen_adapter import COMMIT
from .design_bundle import export_bundle


def run_job(prepared, output, host, port, key, *, timeout=3600):
    prepared, output, key = (Path(p).resolve() for p in (prepared, output, key))
    plan = json.loads((prepared / "boltzgen_plan.json").read_text())
    if plan["commit"] != COMMIT or plan["protocol"] != "protein-anything":
        raise ValueError("Prepared plan does not match the supported pinned protocol")
    command = plan["command"]
    if command[:2] != ["boltzgen", "run"]:
        raise ValueError("Invalid plan command")
    # Execute a reconstructed bounded command, never arbitrary text from an input plan.
    count, length = int(plan["num_designs"]), int(plan["binder_residues"])
    if not 1 <= count <= 32 or not 25 <= length <= 300:
        raise ValueError("Interactive worker supports 1..32 designs, 25..300 residues")
    budget = plan.get("budget", count)
    if isinstance(budget, bool) or not isinstance(budget, int) or not 1 <= budget <= count:
        raise ValueError("Plan budget must be an integer between 1 and num_designs")
    yaml_name = command[2]
    if Path(yaml_name).name != yaml_name or not (prepared / yaml_name).is_file():
        raise ValueError("Design YAML must be a local prepared basename")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output must be empty")
    output.mkdir(parents=True, exist_ok=True)
    if not host or host.startswith("-") or any(c.isspace() for c in host):
        raise ValueError("Invalid SSH host")
    ssh = [
        "ssh",
        "-i",
        str(key),
        "-p",
        str(int(port)),
        "-o",
        "BatchMode=yes",
        "-o",
        "StrictHostKeyChecking=accept-new",
        "root@" + host,
    ]
    scp = [
        "scp",
        "-i",
        str(key),
        "-P",
        str(int(port)),
        "-o",
        "BatchMode=yes",
        "-o",
        "StrictHostKeyChecking=accept-new",
    ]
    remote = "/workspace/jobs/" + uuid.uuid4().hex
    receipt = {"status": "uploading", "remote_directory": remote, "execution_performed": False}

    def save():
        (output / "worker_status.json").write_text(json.dumps(receipt, indent=2))

    save()
    try:
        subprocess.run(ssh + ["mkdir -p " + shlex.quote(remote)], check=True, timeout=60)
        subprocess.run(
            scp + ["-r", str(prepared) + "/.", "root@" + host + ":" + remote],
            check=True,
            timeout=300,
        )
    except (subprocess.SubprocessError, OSError) as exc:
        receipt.update(status="failed_or_interrupted", error=str(exc), remote_cleanup_required=True)
        save()
        return receipt
    model = [
        "/opt/boltzgen-env/bin/boltzgen",
        "run",
        yaml_name,
        "--output",
        "results",
        "--protocol",
        "protein-anything",
        "--num_designs",
        str(count),
        "--budget",
        str(budget),
        "--diffusion_batch_size",
        "1",
        "--inverse_fold_num_sequences",
        "1",
        "--num_workers",
        "0",
        "--devices",
        "1",
        "--use_kernels",
        "false",
        "--cache",
        "/workspace/cache/boltzgen",
    ]
    # The interpreter pointer is produced by our pinned bootstrap, not the incoming plan.
    remote_script = "\n".join(
        [
            "set -euo pipefail",
            "cd " + shlex.quote(remote),
            "export HF_HOME=/workspace/cache/huggingface UV_CACHE_DIR=/tmp/boltzgen-uv-cache",
            "export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=4",
            "export PYTHONUNBUFFERED=1",
            "export PATH=/opt/boltzgen-env/bin:$PATH",
            "esm_python=$(cat /workspace/boltzgen_smoke/esm_python.txt)",
            shlex.join(model) + ' --esmfold2_python "$esm_python"',
        ]
    )
    receipt.update(status="running", execution_performed=True)
    save()
    error = None
    try:
        with (output / "worker.log").open("w") as log:
            subprocess.run(
                ssh + ["bash -s"],
                input=remote_script,
                text=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=timeout,
            )
    except (subprocess.SubprocessError, OSError) as exc:
        error = str(exc)
    # Always collect partial results. A timeout can leave remote execution alive: report it.
    downloaded = output / "downloaded"
    try:
        subprocess.run(
            scp + ["-r", "root@" + host + ":" + remote, str(downloaded)], check=True, timeout=300
        )
    except (subprocess.SubprocessError, OSError) as exc:
        error = (error + "; " if error else "") + "download failed: " + str(exc)
    if error:
        receipt.update(status="failed_or_interrupted", error=error, remote_cleanup_required=True)
        save()
        return receipt
    try:
        bundle = export_bundle(
            downloaded / "results",
            prepared / "target_state.json",
            output / "bundle",
            mode="live_model_outputs",
        )
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        receipt.update(status="failed_export", error=str(exc), remote_cleanup_required=False)
        save()
        return receipt
    receipt.update(
        status="completed",
        bundle="bundle/bundle.json",
        candidates=len(bundle["components"]),
        pod_shutdown_performed=False,
    )
    save()
    return receipt
