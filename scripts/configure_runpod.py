"""Run this yourself in a terminal; the API key is hidden and never printed."""

import getpass
import os
from pathlib import Path

root = Path(__file__).resolve().parents[1]
path = root / ".env.runpod.local"
if path.exists():
    raise SystemExit(
        "Local Runpod credentials already exist. Remove that file yourself to replace them."
    )
key = getpass.getpass("Runpod API key (hidden): ").strip()
if len(key) < 20 or any(c.isspace() for c in key):
    raise SystemExit("Expected a nonempty API key without whitespace.")
fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
with os.fdopen(fd, "w") as f:
    f.write("RUNPOD_API_KEY=" + key + "\n")
print("Saved the API key locally with owner-only permissions. The file is ignored by Git.")
