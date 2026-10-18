"""Build a portable source/demo handoff from an explicit public-file allowlist."""

import hashlib
import json
import zipfile
from pathlib import Path

root = Path(__file__).resolve().parents[1]
files = [
    root / name
    for name in (
        "README.md",
        "pyproject.toml",
        "uv.lock",
        ".gitignore",
        "biocad_molecular_design_implementation_prompt.txt",
        "example_context.txt",
    )
]
for directory in ("src", "tests", "docs", "scripts", "examples"):
    files.extend(
        p
        for p in (root / directory).rglob("*")
        if p.is_file() and "__pycache__" not in p.parts and p.suffix not in {".pyc", ".pyo"}
    )
manifest = {}
with zipfile.ZipFile(root / "biocad_handoff.zip", "w", zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(set(files)):
        if path.is_symlink():
            raise ValueError(f"Unexpected symlink: {path.relative_to(root)}")
        name = str(path.relative_to(root))
        content = path.read_bytes()
        manifest[name] = hashlib.sha256(content).hexdigest()
        archive.writestr(name, content)
    archive.writestr("handoff_manifest.json", json.dumps(manifest, indent=2))
print(
    json.dumps(
        {
            "archive": str(root / "biocad_handoff.zip"),
            "files": len(manifest),
            "bytes": (root / "biocad_handoff.zip").stat().st_size,
        },
        indent=2,
    )
)
