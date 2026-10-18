"""Small HTTP client for the upstream-package to binder-design handoff."""

import json
import time
from hashlib import sha256
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from bespoke_bio.api import Protein, TranslateResponse


class DesignAPIError(RuntimeError):
    pass


def bundle_payload(bundle, *, binder_length=75, num_designs=2, chain_selections=None):
    """Embed only declared, contained structure files from an upstream bundle."""
    root = Path(bundle).resolve()
    provenance = json.loads((root / "provenance.json").read_text())
    proteins = [
        Protein.model_validate_json(line)
        for line in (root / "targets.jsonl").read_text().splitlines()
        if line.strip()
    ]
    if not proteins or len({p.target_id for p in proteins}) != len(proteins):
        raise ValueError("Upstream package needs nonempty, unique target IDs")
    artifacts = {}
    for protein in proteins:
        if (
            protein.sequence_hash
            != "sha256:" + sha256(protein.amino_acid_sequence.encode()).hexdigest()
        ):
            raise ValueError(f"Sequence hash mismatch: {protein.target_id}")
        structure = protein.structure
        if structure and structure.status == "complete" and structure.file:
            path = (root / structure.file).resolve()
            if Path(structure.file).is_absolute() or not path.is_relative_to(root):
                raise ValueError("Structure paths must stay relative and inside the package")
            artifacts[structure.file] = path.read_text()
    response = TranslateResponse(
        case_id=provenance["case_id"],
        reference=provenance["reference"],
        targets=proteins,
        target_count=len(proteins),
        unique_sequence_count=len({p.sequence_hash for p in proteins}),
    )
    return {
        **response.model_dump(mode="json"),
        "artifact_texts": artifacts,
        "provenance": provenance,
        "chain_selections": chain_selections or {},
        "binder_length": binder_length,
        "num_designs": num_designs,
    }


class DesignClient:
    def __init__(self, base_url="http://127.0.0.1:8000", *, timeout=120, opener=urlopen):
        parsed = urlsplit(base_url)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ValueError("Use an http(s) API base URL")
        self.base_url = base_url.rstrip("/")
        self.timeout, self.opener = timeout, opener

    def _request(self, path, payload=None, *, raw=False):
        # Returned artifact URLs must remain on the explicitly selected API origin.
        if not path.startswith("/") or path.startswith("//") or urlsplit(path).netloc:
            raise ValueError("Expected a relative API path")
        request = Request(
            self.base_url + path,
            data=None if payload is None else json.dumps(payload).encode(),
            headers={"Content-Type": "application/json"},
        )
        try:
            with self.opener(request, timeout=self.timeout) as response:
                data = response.read()
        except HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise DesignAPIError(f"API returned {exc.code}: {detail}") from exc
        return data if raw else json.loads(data)

    def submit_bundle(self, bundle, **options):
        return self._request("/v1/designs/upstream", bundle_payload(bundle, **options))

    def status(self, run_id):
        return self._request(self._run_path(run_id))

    def start(self, run_id):
        return self._request(self._run_path(run_id) + "/run", {})

    def results(self, run_id):
        return self._request(self._run_path(run_id) + "/results")

    def download(self, run_id, destination):
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(self._request(self._run_path(run_id) + "/download", raw=True))
        return destination

    def wait(self, run_id, *, timeout=3600, poll_seconds=5):
        if timeout < 0 or poll_seconds <= 0:
            raise ValueError("Timeout must be nonnegative and polling interval positive")
        deadline = time.monotonic() + timeout
        while True:
            result = self.status(run_id)
            if result["status"] not in {"queued", "running"}:
                return result
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Timed out waiting for {run_id}; remote job was not canceled")
            time.sleep(min(poll_seconds, max(0, deadline - time.monotonic())))

    @staticmethod
    def _run_path(run_id):
        if len(run_id) != 32 or any(c not in "0123456789abcdef" for c in run_id):
            raise ValueError("Invalid run ID")
        return "/v1/designs/" + run_id
