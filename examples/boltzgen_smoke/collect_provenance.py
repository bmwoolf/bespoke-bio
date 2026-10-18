"""Record worker identity and model-cache references without copying credentials."""

import importlib.metadata
import json
import subprocess
from pathlib import Path

import torch

root = Path(__file__).resolve().parent
files = []
seen = set()
for folder in [Path("/workspace/cache/boltzgen"), Path("/workspace/cache/huggingface")]:
    for path in folder.rglob("*"):
        if path.is_file() and path.suffix in {".ckpt", ".pt", ".safetensors"}:
            actual = path.resolve()
            if actual in seen:
                continue
            seen.add(actual)
            files.append(
                {
                    "path": str(path),
                    "cache_blob_reference": actual.name,
                    "bytes": actual.stat().st_size,
                }
            )
report = {
    "boltzgen_commit": "a07fea1153111ef453aaea7a15697ee05de9291f",
    "packages": {
        k: importlib.metadata.version(k) for k in ["boltzgen", "torch", "numpy", "biotite", "gemmi"]
    },
    "cuda": torch.version.cuda,
    "gpu": torch.cuda.get_device_name(0),
    "hardware": subprocess.check_output(
        ["nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"],
        text=True,
    ).strip(),
    "model_cache_references": files,
    "note": "Cache blob references are provenance, not separately rehashed weight verification.",
}
(root / "worker_provenance.json").write_text(json.dumps(report, indent=2) + "\n")
