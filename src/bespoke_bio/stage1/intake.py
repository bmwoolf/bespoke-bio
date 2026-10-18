"""Acquire pinned files and validate complete paired FASTQs with the native fq parser."""

import hashlib
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlsplit

from .contracts import Artifact, Intake, ReadInput
from .io import canonical, read_json, sync_directory, verify, write_new


def hashes(path: Path, algorithm: str) -> tuple[str, str]:
    """Compute the publisher checksum and SHA-256 together in one bounded read pass."""
    publisher, sha = hashlib.new(algorithm), hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(16 * 1024 * 1024):
            publisher.update(block)
            sha.update(block)
    return publisher.hexdigest(), sha.hexdigest()


def fetch_file(source: dict, directory: Path) -> Artifact:
    """Download with three transient-error attempts; a checksum mismatch never retries."""
    filename = source["filename"]
    if Path(filename).name != filename or urlsplit(source["url"]).scheme != "https":
        raise ValueError("source requires a plain filename and HTTPS URL")
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / filename
    registration = directory / f"{filename}.registered.json"
    if registration.exists():
        record = Artifact.model_validate(read_json(registration))
        verify(record)
        return record
    if destination.exists():
        raise FileExistsError(f"unregistered destination exists: {filename}")
    partial = directory / (filename + ".partial")
    started = time.monotonic()
    command = [
        "curl",
        "--fail",
        "--location",
        "--proto",
        "=https",
        "--tlsv1.2",
        "--retry",
        "2",
        "--retry-delay",
        "5",
        "--connect-timeout",
        "30",
        "--speed-time",
        "180",
        "--speed-limit",
        "1024",
        "--continue-at",
        "-",
        "--output",
        str(partial),
        source["url"],
    ]
    with (directory / (filename + ".transfer.log")).open("ab") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)
    if source.get("bytes") is not None and partial.stat().st_size != source["bytes"]:
        raise ValueError(f"source size mismatch: {filename}")
    publisher, sha = hashes(partial, source["checksum_algorithm"])
    if publisher != source["checksum"]:
        raise ValueError(f"source checksum mismatch: {filename}")
    with partial.open("rb") as stream:
        os.fsync(stream.fileno())
    os.link(partial, destination)
    destination.chmod(0o444)
    partial.unlink()
    sync_directory(directory)
    item = Artifact(
        role=source["role"],
        location=str(destination.resolve()),
        bytes=destination.stat().st_size,
        sha256=sha,
        format=filename.split(".")[-1],
    )
    write_new(registration, canonical(item.model_dump(mode="json")) + b"\n")
    write_new(
        directory / (filename + ".transfer.json"),
        canonical(
            {
                "source": source["url"],
                "seconds": time.monotonic() - started,
                "bytes": item.bytes,
                "sha256": sha,
                "command": command,
            }
        )
        + b"\n",
    )
    return item


def fetch_locked(lock: dict, root: Path, workers: int = 4) -> list[Artifact]:
    """Stage sources once; truth files occupy a separate directory from processing inputs."""
    if not 1 <= workers <= 4:
        raise ValueError("download concurrency must be between one and four")

    def fetch(source):
        role = source["role"]
        folder = (
            "benchmarks"
            if role.startswith("truth_")
            else ("inputs" if role.startswith(("dna_", "rna_")) else "references")
        )
        return fetch_file(source, root / folder)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(fetch, lock["files"]))


def pilot_intake(lock: dict, artifacts: list[Artifact], case_id: str = "HG002-pilot") -> Intake:
    """Describe the locked delivery without claiming matched aliquots or original lane IDs."""
    by_role = {a.role: a for a in artifacts}
    data = lock["dataset"]
    reads = []
    for assay, prefix, protocol in [
        ("germline_wgs", "dna", "Illumina-PCR-free-WGS"),
        ("bulk_rna", "rna", "Illumina-polyA-2x150"),
    ]:
        for mate in (1, 2):
            reads.append(
                ReadInput(
                    donor=data["donor"],
                    specimen=data[f"{prefix}_specimen"],
                    library=data[f"{prefix}_library"],
                    assay=assay,
                    lane="delivery-1",
                    mate=mate,
                    protocol=protocol,
                    passage=None,
                    aliquot=None,
                    artifact=by_role[f"{prefix}_r{mate}"],
                )
            )
    return Intake(case_id=case_id, sex=data["sex"], reads=reads)


def validate_pairs(intake: Intake, output: Path, executable: list[str]) -> list[dict]:
    """Run fq lint with its paired-name P001 check enabled; retain the full native report."""
    output.mkdir(parents=True, exist_ok=True)
    results = []
    groups = sorted({(r.library, r.lane) for r in intake.reads})
    for number, (library, lane) in enumerate(groups):
        pair = sorted(
            [r for r in intake.reads if (r.library, r.lane) == (library, lane)],
            key=lambda r: r.mate,
        )
        paths = [str(verify(r.artifact)) for r in pair]
        command = [*executable, "lint", *paths]
        log = output / f"pair-{number}.fq-lint.txt"
        with log.open("xb") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=False)
        entry = {
            "library": library,
            "lane": lane,
            "command": command,
            "exit_code": result.returncode,
            "log": str(log),
            "paired_name_check": "P001",
        }
        results.append(entry)
        if result.returncode:
            raise ValueError(f"FASTQ pairing/format failure: {library}/{lane}; see {log}")
    write_new(output / "pair-checks.json", json.dumps(results, indent=2).encode() + b"\n")
    return results
