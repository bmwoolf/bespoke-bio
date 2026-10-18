"""Invoke the locally pinned CLI with credentials loaded without echoing them."""

import os
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
env = os.environ.copy()
secret = root / ".env.runpod.local"
if not env.get("RUNPOD_API_KEY") and secret.exists():
    line = secret.read_text().strip()
    name, sep, value = line.partition("=")
    if name != "RUNPOD_API_KEY" or not sep or not value or "\n" in value:
        raise SystemExit("Invalid local credential file")
    env["RUNPOD_API_KEY"] = value
if not sys.argv[1:]:
    raise SystemExit("Usage: python scripts/runpod.py user | gpu list | pod ...")
if sys.argv[1] == "config":
    raise SystemExit("Use scripts/configure_runpod.py to keep secrets out of command arguments.")
result = subprocess.run(
    [str(root / "vendor/runpod/runpodctl"), *sys.argv[1:]], env=env, cwd=root, check=False
)
raise SystemExit(result.returncode)
