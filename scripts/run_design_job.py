"""Run one prepared job on an existing, warmed SSH worker; never provisions compute.

The operator owns the pod lifetime/cost guard. Failed jobs still download their logs.
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from biocad.worker import run_job


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared")
    p.add_argument("--output", required=True)
    p.add_argument("--host", required=True)
    p.add_argument("--port", type=int, required=True)
    p.add_argument("--key", default=str(ROOT / "vendor/runpod/ssh/smoke"))
    p.add_argument("--timeout", type=int, default=3600)
    a = p.parse_args()
    result = run_job(a.prepared, a.output, a.host, a.port, a.key, timeout=a.timeout)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
