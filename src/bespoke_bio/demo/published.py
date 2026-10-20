"""Explicit, checksummed snapshots of completed jobs for network-free public serving."""

import argparse
import json
import re
import shutil
import tempfile
from hashlib import sha256
from pathlib import Path
from urllib.parse import unquote, urlsplit
from urllib.request import Request, build_opener

from .client import NoRedirect

JOB_ID = re.compile(r"^[a-f0-9]{32}$")
COMPLETE = {"completed", "cached_demo_completed"}


def artifact_key(url, job_id):
    if not isinstance(url, str):
        raise TypeError("Invalid published artifact URL")
    path = unquote(url)
    prefix = f"/api/jobs/{job_id}/"
    if (
        not path.startswith(prefix)
        or (path != prefix + "download" and not path.startswith(prefix + "artifacts/"))
        or any(part in {"", ".", ".."} for part in path[1:].split("/"))
        or any(char in path for char in "\\?#\r\n")
    ):
        raise ValueError("Artifact does not belong to this published job")
    return path


def result_urls(result):
    urls = [result["download_url"]]
    if result.get("bundle_url"):
        urls.append(result["bundle_url"])
    for candidate in result["candidates"]:
        urls.extend(candidate[field] for field in ("binder_structure_url", "complex_structure_url"))
        if candidate.get("target_structure_url"):
            urls.append(candidate["target_structure_url"])
        for comparison in candidate.get("comparison_scores", []):
            urls.extend(
                comparison[field]
                for field in (
                    "target_structure_url",
                    "binder_structure_url",
                    "complex_structure_url",
                )
            )
        for prediction in candidate.get("predictions", []):
            urls.extend(
                prediction[field]
                for field in (
                    "target_structure_url",
                    "binder_structure_url",
                    "complex_structure_url",
                    "evidence_url",
                )
            )
    return urls


def export_snapshot(destination, get_json, get_bytes, job_ids):
    """Publish only explicitly selected completed jobs; never copy the private runs directory."""
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".published-", dir=destination.parent))
    try:
        catalog = get_json("/api/catalog")
        manifest = {
            "version": 1,
            "catalog_sha256": catalog["manifest_sha256"],
            "jobs": {},
            "files": {},
        }
        (temporary / "files").mkdir()
        for job_id in dict.fromkeys(job_ids):
            if not JOB_ID.fullmatch(job_id):
                raise ValueError("Invalid job ID")
            job = get_json(f"/api/jobs/{job_id}")
            if job["job_id"] != job_id or job["status"] not in COMPLETE:
                raise ValueError("Only completed jobs can be published")
            result = get_json(f"/api/jobs/{job_id}/results")
            if result["target_id"] != job["target_id"]:
                raise ValueError("Result target differs from its job")
            # Public jobs carry biological/result identity, never operator errors or worker progress.
            job = {
                key: job.get(key)
                for key in (
                    "job_id",
                    "target_id",
                    "status",
                    "mode",
                    "execution_performed",
                    "dispatch_attempted",
                )
            }
            job.update(progress=None, error=None)
            manifest["jobs"][job_id] = {"job": job, "results": result}
            for url in result_urls(result):
                key = artifact_key(url, job_id)
                if key in manifest["files"]:
                    continue
                data = get_bytes(url)
                digest = sha256(data).hexdigest()
                (temporary / "files" / digest).write_bytes(data)
                manifest["files"][key] = {"sha256": digest, "size": len(data)}
        (temporary / "manifest.json").write_text(json.dumps(manifest, allow_nan=False) + "\n")
        PublishedSnapshot(temporary, catalog["manifest_sha256"])
        temporary.rename(destination)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


class PublishedSnapshot:
    def __init__(self, root, catalog_hash):
        self.root = Path(root).resolve()
        manifest = json.loads((self.root / "manifest.json").read_text())
        if manifest["version"] != 1 or manifest["catalog_sha256"] != catalog_hash:
            raise ValueError("Published snapshot does not match this catalog")
        self.jobs, self.files = manifest["jobs"], manifest["files"]
        for job_id, entry in self.jobs.items():
            if not JOB_ID.fullmatch(job_id) or entry["job"]["job_id"] != job_id:
                raise ValueError("Invalid published job identity")
            if entry["job"]["status"] not in COMPLETE:
                raise ValueError("Published job is incomplete")
            if entry["results"]["target_id"] != entry["job"]["target_id"]:
                raise ValueError("Published result belongs to another target")
            design = entry.get("design")
            if design is not None:
                if not isinstance(design.get("name"), str) or not design["name"].strip():
                    raise ValueError("Published design needs its source target name")
                if design.get("target_id") != entry["job"]["target_id"]:
                    raise ValueError("Published design metadata belongs to another target")
                sequence = design.get("amino_acid_sequence")
                if not isinstance(sequence, str) or not sequence:
                    raise ValueError("Published design needs its complete input sequence")
                if sha256(sequence.encode()).hexdigest() != design.get("sequence_sha256"):
                    raise ValueError("Published design sequence checksum mismatch")
            for url in result_urls(entry["results"]):
                if artifact_key(url, job_id) not in self.files:
                    raise ValueError("Published artifact is missing")
        for url, entry in self.files.items():
            job_id = url.split("/")[3]
            if job_id not in self.jobs or artifact_key(url, job_id) != url:
                raise ValueError("Unlisted published job")
            if not re.fullmatch(r"[a-f0-9]{64}", entry["sha256"]):
                raise ValueError("Invalid artifact checksum")
            path = self.root / "files" / entry["sha256"]
            if path.is_symlink() or path.stat().st_size != entry["size"]:
                raise ValueError("Published artifact size or location changed")
            if sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
                raise ValueError("Published artifact checksum mismatch")

    def file(self, url):
        entry = self.files[url]
        return self.root / "files" / entry["sha256"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway", default="http://127.0.0.1:8787")
    parser.add_argument("--job", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    origin = urlsplit(args.gateway)
    if origin.scheme != "http" or origin.hostname not in {"127.0.0.1", "localhost", "::1"}:
        parser.error("Export from the local operator gateway through loopback or an SSH tunnel")
    if origin.path not in {"", "/"} or origin.query or origin.fragment or origin.username:
        parser.error("Use a plain loopback origin")
    opener = build_opener(NoRedirect())

    def fetch(path, limit):
        if not path.startswith("/api/") or path.startswith("//"):
            raise ValueError("Invalid gateway path")
        with opener.open(Request(args.gateway.rstrip("/") + path), timeout=60) as response:
            data = response.read(limit + 1)
        if len(data) > limit:
            raise ValueError("Published response exceeds the size limit")
        return data

    export_snapshot(
        args.output,
        lambda path: json.loads(fetch(path, 8 * 1024 * 1024)),
        lambda path: fetch(path, 128 * 1024 * 1024),
        args.job,
    )
    print(f"Published {len(set(args.job))} completed jobs to {args.output}")


if __name__ == "__main__":
    main()
