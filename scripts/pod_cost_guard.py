"""Stop only a task-named pod at a deadline if its termination receipt is absent."""

import argparse
import subprocess
import time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("pod_id")
p.add_argument("output")
p.add_argument("--seconds", type=int, default=10800)
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
output = Path(a.output).resolve()
if not output.is_relative_to(root / "outputs"):
    raise SystemExit("Guard receipts must be under workspace outputs")
deadline = time.monotonic() + a.seconds
while time.monotonic() < deadline:
    if (output / "terminated.json").exists():
        raise SystemExit(0)
    time.sleep(min(30, max(0, deadline - time.monotonic())))
if not (output / "terminated.json").exists():
    with (output / "cost_guard_stop.json").open("w") as f:
        subprocess.run(
            ["python3", str(root / "scripts/runpod.py"), "pod", "stop", a.pod_id],
            cwd=root,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=False,
        )
