"""Strict JSON, content identities and durable publication on a local POSIX filesystem."""

import fcntl
import hashlib
import json
import math
import os
from contextlib import contextmanager
from pathlib import Path

from .contracts import Artifact, Intake


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _constant(value):
    raise ValueError(f"non-finite JSON number: {value}")


def _float(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite JSON number: {value}")
    return result


def read_json(path: Path):
    """Read JSON without accepting duplicate keys or JavaScript numeric constants."""
    return decode_json(path.read_text())


def decode_json(text: str):
    """Decode one JSON value with the same strict rules used for manifests."""
    return json.loads(text, object_pairs_hook=_pairs, parse_constant=_constant, parse_float=_float)


def canonical(value) -> bytes:
    """Encode deterministic UTF-8 JSON; reject numbers outside the JSON domain."""
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")


def digest(value) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def file_hash(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def artifact(path: Path, role: str, format: str, **kwargs) -> Artifact:
    return Artifact(
        role=role,
        location=str(path.resolve()),
        bytes=path.stat().st_size,
        sha256=file_hash(path),
        format=format,
        **kwargs,
    )


def verify(item: Artifact) -> Path:
    """Verify a staged local artifact at registration, transfer, restore or publication."""
    path = Path(item.location)
    if not path.is_file() or path.stat().st_size != item.bytes or file_hash(path) != item.sha256:
        raise ValueError(f"artifact integrity failure: {item.role}")
    return path


def source_hash(root: Path) -> str:
    """Hash delivered source bytes, including uncommitted changes, in sorted relative paths."""
    files = sorted(p for p in root.rglob("*.py") if "__pycache__" not in p.parts)
    if not files:
        raise ValueError("source tree is empty")
    return digest({str(p.relative_to(root)): file_hash(p) for p in files})


def run_identity(
    intake: Intake,
    *,
    reference: str,
    config: str,
    tools: str,
    adapter: str,
    step2_source: str,
    contract: str,
) -> str:
    """Locations do not change biological identity; exact content and code always do."""
    payload = intake.model_dump(mode="json")
    for item in payload["reads"]:
        item["artifact"].pop("location")
    payload["reads"].sort(key=lambda x: (x["library"], x["lane"], x["mate"]))
    return digest(
        {
            "intake": payload,
            "reference": reference,
            "config": config,
            "tools": tools,
            "adapter": adapter,
            "step2_source": step2_source,
            "contract": contract,
        }
    )


def sync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_new(path: Path, data: bytes) -> None:
    """Create and flush an immutable file; fail instead of replacing existing evidence."""
    with path.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    sync_directory(path.parent)


def publish_ready(attempt: Path, destination: Path, value) -> None:
    """Link a flushed manifest atomically without replacement, on the same filesystem."""
    staging = attempt / "manifest.json"
    write_new(staging, canonical(value) + b"\n")
    os.link(staging, destination)  # EXDEV and EEXIST deliberately fail closed.
    sync_directory(destination.parent)


@contextmanager
def run_lock(path: Path):
    """Hold one nonblocking OS lock for the entire coordinator lifetime."""
    with path.open("a+") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def write_jsonl(path: Path, records) -> None:
    """Write records in bounded memory and flush before allowing publication."""
    with path.open("xb") as stream:
        for record in records:
            stream.write(canonical(record) + b"\n")
        stream.flush()
        os.fsync(stream.fileno())
    sync_directory(path.parent)
