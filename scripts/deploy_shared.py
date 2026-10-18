"""Build a secret-free service bundle; optionally deploy it to an existing SSH pod.

Private config is a JSON mapping of environment names to strings. It is transferred
separately, protected with mode 0600, never printed, and never included in the bundle.
"""

import argparse
import hashlib
import json
import math
import subprocess
import tarfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REMOTE = "/workspace/shared-service"
REMOTE_SECRETS = "/root/.config/biocad-shared"
LINUX_CLI_SHA = "b66d11eca8cacbf18a2de6943c6e18a79ffbdc8a3c30aa13a49dc6cfc15d4239"

BOOTSTRAP = r"""set -euo pipefail
umask 077
cd /workspace/shared-service
mkdir -p runs logs state
mkdir -p /opt/biocad-service/app /opt/biocad-service/bin
mkdir -p /root/.config/biocad-shared
chmod 700 /root/.config/biocad-shared
chmod 600 /root/.config/biocad-shared/environment.json /root/.config/biocad-shared/worker_key
if ! command -v python3 >/dev/null || ! python3 -m venv /tmp/biocad-venv-probe >/dev/null 2>&1; then
  apt-get update -qq
  DEBIAN_FRONTEND=noninteractive apt-get install -y -qq python3 python3-venv openssh-client ca-certificates
fi
python3 - <<'INNER'
import hashlib, pathlib, urllib.request
url = 'https://github.com/runpod/runpodctl/releases/download/v3.0.1/runpodctl-linux-amd64'
p = pathlib.Path('/opt/biocad-service/bin/runpodctl')
if not p.exists():
    urllib.request.urlretrieve(url, p)
assert hashlib.sha256(p.read_bytes()).hexdigest() == 'CLI_SHA', 'CLI checksum mismatch'
p.chmod(0o700)
INNER
rm -rf /tmp/biocad-venv-probe
tar -xzf service.tar.gz -C /opt/biocad-service/app
python3 -m venv /opt/biocad-service/runtime
/opt/biocad-service/runtime/bin/pip install --disable-pip-version-check '/opt/biocad-service/app[design]'
# Stop only our previous wrapper/child on redeploy; retain the lifetime guard.
python3 - <<'INNER'
import os, pathlib, signal, time
for name in ['service.pid','child.pid']:
    p=pathlib.Path(name)
    if not p.exists(): continue
    try:
        pid=int(p.read_text())
        command=pathlib.Path('/proc')/str(pid)/'cmdline'
        text=command.read_bytes() if command.exists() else b''
        if b'start_service.py' in text or b'biocad.service:app' in text:
            os.kill(pid,signal.SIGTERM)
    except (ValueError,ProcessLookupError): pass
time.sleep(12)
for name in ['service.pid','child.pid']:
    p=pathlib.Path(name)
    if not p.exists(): continue
    try:
        pid=int(p.read_text())
        command=pathlib.Path('/proc')/str(pid)/'cmdline'
        text=command.read_bytes() if command.exists() else b''
        if b'start_service.py' in text or b'biocad.service:app' in text:
            os.kill(pid,signal.SIGKILL)
    except (ValueError,ProcessLookupError): pass
INNER
cat > /opt/biocad-service/start_service.py <<'INNER'
import json, os, pathlib, subprocess, time
root=pathlib.Path('/workspace/shared-service')
exec_root=pathlib.Path('/opt/biocad-service')
secret_root=pathlib.Path('/root/.config/biocad-shared')
config=json.loads((secret_root/'environment.json').read_text())
os.environ.update({k:str(v) for k,v in config.items()})
os.environ['BIOCAD_REPO_DIR']=str(exec_root/'app')
os.environ['BIOCAD_RUNS_DIR']=str(root/'runs')
os.environ['BIOCAD_WORKER_KEY']=str(secret_root/'worker_key')
os.environ['BIOCAD_SSH_KEY']=str(secret_root/'worker_key')
os.environ['RUNPODCTL']=str(exec_root/'bin/runpodctl')
os.environ['BIOCAD_CLOUD_STATE_DIR']=str(root/'cloud_leases')
os.environ['PATH']=str(exec_root/'bin')+':'+os.environ['PATH']
# Lifetime starts once, survives process restarts. A stopped/deleted gateway cannot
# execute this guard; the operator must also verify final cloud cleanup externally.
deadline_path=root/'state/gateway_deadline'
if not deadline_path.exists():
    deadline_path.write_text(str(float(config.get('BIOCAD_SERVICE_EXPIRES_AT',time.time()+86400))))
if time.time() >= float(deadline_path.read_text()):
    raise SystemExit('Gateway lifetime expired')
os.chdir(exec_root/'app')
import fcntl, signal
lock=(secret_root/'service.lock').open('w')
try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError: raise SystemExit('Service supervisor already running')
stopping=False
child=None
def stop(signum,frame):
    global stopping
    stopping=True
    if child and child.poll() is None: child.terminate()
signal.signal(signal.SIGTERM,stop)
signal.signal(signal.SIGINT,stop)
while not stopping and time.time()<float(deadline_path.read_text()):
    child=subprocess.Popen([str(exec_root/'runtime/bin/uvicorn'),'biocad.service:app','--host','0.0.0.0','--port','8000','--timeout-graceful-shutdown','10'])
    (root/'child.pid').write_text(str(child.pid))
    child.wait()
    if not stopping: time.sleep(2)
INNER
cat > /opt/biocad-service/lifetime_guard.py <<'INNER'
import json, os, pathlib, subprocess, time
root=pathlib.Path('/workspace/shared-service')
exec_root=pathlib.Path('/opt/biocad-service')
secret_root=pathlib.Path('/root/.config/biocad-shared')
config=json.loads((secret_root/'environment.json').read_text())
os.environ.update({k:str(v) for k,v in config.items()})
import fcntl, signal
lock=(secret_root/'lifetime.lock').open('w')
try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError: raise SystemExit('Lifetime guard already running')
while not (root/'state/gateway_deadline').exists(): time.sleep(1)
deadline=float((root/'state/gateway_deadline').read_text())
while time.time()<deadline: time.sleep(min(30,deadline-time.time()))
os.environ['RUNPODCTL']=str(exec_root/'bin/runpodctl')
os.environ['BIOCAD_CLOUD_STATE_DIR']=str(root/'cloud_leases')
from biocad.cloud import reconcile_owned_leases
# Stop admission/execution before reconciling, including provisioning leases whose
# create request was interrupted before its pod ID returned.
for name in ['service.pid','child.pid']:
    try:
        pid=int((root/name).read_text())
        command=pathlib.Path('/proc')/str(pid)/'cmdline'
        text=command.read_bytes() if command.exists() else b''
        if b'start_service.py' in text or b'biocad.service:app' in text: os.kill(pid,signal.SIGTERM)
    except (FileNotFoundError,ProcessLookupError,ValueError): pass
time.sleep(2)
while reconcile_owned_leases(root/'cloud_leases'):
    time.sleep(10)
pod=config.get('BIOCAD_GATEWAY_POD_ID')
if pod:
    # Run state remains on the mounted network volume after gateway deletion.
    with (root/'logs/lifetime_stop.json').open('w') as out:
        subprocess.run([str(exec_root/'bin/runpodctl'),'pod','delete',pod],stdout=out,stderr=subprocess.STDOUT,check=False)
INNER
nohup /opt/biocad-service/runtime/bin/python /opt/biocad-service/start_service.py > logs/service.log 2>&1 < /dev/null &
echo $! > service.pid
if [ -f lifetime.pid ] && kill -0 "$(cat lifetime.pid)" 2>/dev/null && grep -zq lifetime_guard.py "/proc/$(cat lifetime.pid)/cmdline"; then
  : # Preserve the original independent lifetime guard and its PID.
else
  nohup /opt/biocad-service/runtime/bin/python /opt/biocad-service/lifetime_guard.py > logs/lifetime.log 2>&1 < /dev/null &
  echo $! > lifetime.pid
fi
""".replace("CLI_SHA", LINUX_CLI_SHA)


def build(output):
    """Allowlist source/artifact paths rather than trying to blacklist secret names."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    roots = [
        ROOT / "src",
        ROOT / "pyproject.toml",
        ROOT / "README.md",
        ROOT / "examples/5JDS.cif",
        ROOT / "examples/boltzgen_results",
        ROOT / "examples/boltzgen_smoke",
        ROOT / "scripts/runpod.py",
    ]
    with tarfile.open(output, "w:gz") as archive:
        for root in roots:
            if not root.exists():
                continue
            paths = root.rglob("*") if root.is_dir() else [root]
            for path in paths:
                if not path.is_file() or path.is_symlink() or "__pycache__" in path.parts:
                    continue
                if path.suffix in {".pyc", ".zip", ".tar", ".gz"}:
                    continue
                archive.add(path, arcname=path.relative_to(ROOT), recursive=False)
    return {"path": str(output), "sha256": hashlib.sha256(output.read_bytes()).hexdigest()}


def deploy(bundle, host, port, key, private_config, worker_key, binary_path=None):
    if not host or host.startswith("-") or any(c.isspace() for c in host):
        raise ValueError("Invalid SSH host")
    config_path = Path(private_config)
    config = json.loads(config_path.read_text())
    if not isinstance(config, dict) or any(not isinstance(k, str) for k in config):
        raise ValueError("Private config must be an environment mapping")
    if str(config.get("BIOCAD_MANAGED_WORKERS")) == "1":
        required = [
            "BIOCAD_API_TOKEN",
            "RUNPOD_API_KEY",
            "BIOCAD_GATEWAY_POD_ID",
            "BIOCAD_SERVICE_EXPIRES_AT",
        ]
        if not all(config.get(name) for name in required):
            raise ValueError("Managed deployment requires token, API key, gateway ID and expiry")
        expires = float(config["BIOCAD_SERVICE_EXPIRES_AT"])
        if not math.isfinite(expires) or not time.time() < expires <= time.time() + 86400:
            raise ValueError("Managed gateway expiry must be within the next 24 hours")
    if not 1 <= port <= 65535:
        raise ValueError("Invalid SSH port")
    ssh = [
        "ssh",
        "-i",
        str(key),
        "-p",
        str(port),
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
        str(port),
        "-o",
        "BatchMode=yes",
        "-o",
        "StrictHostKeyChecking=accept-new",
    ]
    subprocess.run(
        ssh
        + [
            "umask 077; mkdir -p " + REMOTE + " " + REMOTE_SECRETS + "; chmod 700 " + REMOTE_SECRETS
        ],
        check=True,
    )
    for source, destination in [
        (bundle, REMOTE + "/service.tar.gz"),
        (private_config, REMOTE_SECRETS + "/environment.json"),
        (worker_key, REMOTE_SECRETS + "/worker_key"),
    ]:
        subprocess.run(scp + [str(source), "root@" + host + ":" + destination], check=True)
    if binary_path:
        if hashlib.sha256(Path(binary_path).read_bytes()).hexdigest() != LINUX_CLI_SHA:
            raise ValueError("Local Linux CLI checksum mismatch")
        subprocess.run(ssh + ["mkdir -p /opt/biocad-service/bin"], check=True)
        subprocess.run(
            scp + [str(binary_path), "root@" + host + ":/opt/biocad-service/bin/runpodctl"],
            check=True,
        )
    subprocess.run(ssh + ["bash -s"], input=BOOTSTRAP, text=True, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="outputs/shared_deploy/service.tar.gz")
    parser.add_argument("--host")
    parser.add_argument("--port", type=int, default=22)
    parser.add_argument("--ssh-key")
    parser.add_argument("--private-config", "--env-file", dest="private_config")
    parser.add_argument("--binary-path")
    parser.add_argument("--worker-key")
    args = parser.parse_args()
    receipt = build(args.output)
    if args.host:
        if not all([args.ssh_key, args.private_config, args.worker_key]):
            parser.error("Deploy requires --ssh-key, --private-config and --worker-key")
        deploy(
            args.output,
            args.host,
            args.port,
            args.ssh_key,
            args.private_config,
            args.worker_key,
            args.binary_path,
        )
        receipt["deployed_to_existing_host"] = args.host
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
