"""Loopback-only demo gateway; serve through an SSH tunnel, with one process."""

import json
import os
import re
import threading
from hashlib import sha256
from pathlib import Path
from urllib.parse import quote, urlsplit
from uuid import uuid4

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field
from starlette.middleware.trustedhost import TrustedHostMiddleware

from .catalog import Catalog
from .client import DesignService, ServiceError

STATIC = Path(__file__).parent / "static"
JOB_ID = re.compile(r"^[a-f0-9]{32}$")


class Selection(BaseModel):
    model_config = ConfigDict(extra="forbid")
    target_id: str = Field(min_length=1, max_length=200)
    binder_length: int = Field(default=75, ge=25, le=300)
    num_designs: int = Field(default=2, ge=1, le=32)


def sync_directory(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_json(path, value):
    temporary = path.with_suffix(".partial")
    with temporary.open("w") as stream:
        stream.write(json.dumps(value, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    sync_directory(path.parent)


def create_app(index=None, runs=None, service=None, *, allow_run=False):
    app = FastAPI(title="Bespoke Bio protein demo", version="1.0.0")
    app.add_middleware(
        TrustedHostMiddleware, allowed_hosts=["localhost", "127.0.0.1", "[::1]", "testserver"]
    )
    catalog = Catalog(Path(index)) if index else None
    run_root = Path(runs).resolve() if runs else None
    lock = threading.Lock()

    @app.middleware("http")
    async def same_origin(request: Request, call_next):
        if request.method not in {"GET", "HEAD", "OPTIONS"}:
            origin = request.headers.get("origin")
            if request.headers.get("x-bespoke-request") != "1" or (
                origin and urlsplit(origin).netloc != request.headers.get("host")
            ):
                return JSONResponse({"detail": "Use the local demo client"}, status_code=403)
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Cache-Control"] = "no-store"
        return response

    @app.exception_handler(ServiceError)
    async def upstream_error(request, exc):
        return JSONResponse({"detail": exc.message}, status_code=exc.status)

    def require_catalog():
        if not catalog or not catalog.path.is_file():
            raise HTTPException(503, "Build and configure the protein catalog first")
        return catalog

    def require_service():
        if service is None or run_root is None:
            raise HTTPException(503, "Configure the design API and a private runs directory")
        run_root.mkdir(parents=True, exist_ok=True, mode=0o700)
        return service

    def target(target_id):
        try:
            return require_catalog().get(target_id)
        except KeyError as exc:
            raise HTTPException(404, "Unknown target") from exc

    def receipt(job_id):
        if not JOB_ID.fullmatch(job_id) or run_root is None:
            raise HTTPException(404, "Unknown local job")
        path = run_root / (job_id + ".json")
        if not path.is_file():
            raise HTTPException(404, "Unknown local job")
        record = json.loads(path.read_text())
        if record["status"] == "submitting" and not record.get("remote_id"):
            record["status"] = "submission_unknown"
        return path, record

    def remote_path(record):
        remote_id = record.get("remote_id", "")
        if not JOB_ID.fullmatch(remote_id):
            raise HTTPException(409, "Submission outcome is unknown; do not automatically resubmit")
        return "/v1/designs/" + remote_id

    def envelope(record):
        return {
            k: record.get(k)
            for k in (
                "job_id",
                "target_id",
                "status",
                "mode",
                "execution_performed",
                "progress",
                "candidates",
                "native_filter_passed",
                "error",
                "dispatch_attempted",
            )
        }

    def update(record, remote):
        if remote.get("run_id") != record.get("remote_id"):
            raise ServiceError(502, "Response belongs to a different design job")
        for key in (
            "status",
            "mode",
            "execution_performed",
            "progress",
            "candidates",
            "native_filter_passed",
        ):
            if key in remote:
                record[key] = remote[key]
        return record

    def submit(payload, *, cached=False):
        upstream = require_service()
        fingerprint = sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        if cached:
            fingerprint = "cached-demo"
        marker = run_root / (fingerprint + ".request")
        with lock:
            if marker.exists():
                _, saved = receipt(marker.read_text().strip())
                return envelope(saved)
            job_id = uuid4().hex
            record = {
                "job_id": job_id,
                "target_id": None if cached else payload["target_id"],
                "status": "submitting",
                "mode": "cached_demo" if cached else "prepare",
                "execution_performed": False,
                "submitted_record": payload,
                "manifest_sha256": None
                if cached
                else require_catalog().metadata()["manifest_sha256"],
            }
            path = run_root / (job_id + ".json")
            write_json(path, record)
            # Persist before POST; a timeout never silently causes a duplicate submission.
            with marker.open("x") as stream:
                stream.write(job_id + "\n")
                stream.flush()
                os.fsync(stream.fileno())
            sync_directory(run_root)
            try:
                response = upstream.request(
                    "/v1/designs/demo" if cached else "/v1/designs", payload
                )
                if not JOB_ID.fullmatch(response.get("run_id", "")):
                    raise ServiceError(502, "Design service returned an invalid run ID")
                record["remote_id"] = response["run_id"]
                update(record, response)
            except ServiceError as exc:
                definite = exc.status in {401, 403, 404, 422}
                record.update(
                    status="rejected" if definite else "submission_unknown", error=exc.message
                )
                write_json(path, record)
                if definite:
                    marker.unlink()
                raise
            write_json(path, record)
            return envelope(record)

    @app.get("/api/catalog")
    def catalog_info():
        return {
            **require_catalog().metadata(),
            "service_configured": service is not None,
            "gpu_execution_enabled": allow_run,
        }

    @app.get("/api/targets")
    def search(
        q: str = Query("", max_length=200),
        changed: bool = False,
        limit: int = Query(25, ge=1, le=100),
        offset: int = Query(0, ge=0),
    ):
        return require_catalog().search(q.strip(), changed, limit, offset)

    @app.get("/api/targets/{target_id}")
    def target_record(target_id: str):
        return target(target_id)

    @app.post("/api/prepare", status_code=201)
    def prepare(selection: Selection):
        payload = {
            **target(selection.target_id),
            "binder_length": selection.binder_length,
            "num_designs": selection.num_designs,
        }
        return submit(payload)

    @app.post("/api/demo", status_code=201)
    def cached_demo():
        return submit({}, cached=True)

    @app.get("/api/jobs")
    def jobs():
        if run_root is None or not run_root.exists():
            return []
        return [envelope(json.loads(p.read_text())) for p in sorted(run_root.glob("*.json"))]

    @app.get("/api/jobs/{job_id}")
    def status(job_id: str):
        with lock:
            path, record = receipt(job_id)
            if record.get("remote_id"):
                update(record, require_service().request(remote_path(record)))
                write_json(path, record)
            return envelope(record)

    @app.post("/api/jobs/{job_id}/run")
    def run(job_id: str):
        if not allow_run:
            raise HTTPException(403, "GPU execution is disabled for integration testing")
        with lock:
            path, record = receipt(job_id)
            if record["status"] != "prepared_pending_gpu_execution" or record.get(
                "dispatch_attempted"
            ):
                raise HTTPException(409, "Only a prepared job can be started once")
            record["status"] = "dispatch_requested"
            record["dispatch_attempted"] = True
            write_json(path, record)
            update(record, require_service().request(remote_path(record) + "/run", {}))
            write_json(path, record)
            return envelope(record)

    @app.post("/api/jobs/{job_id}/cancel")
    def cancel(job_id: str):
        with lock:
            path, record = receipt(job_id)
            update(record, require_service().request(remote_path(record) + "/cancel", {}))
            write_json(path, record)
            return envelope(record)

    @app.get("/api/jobs/{job_id}/results")
    def results(job_id: str):
        _, record = receipt(job_id)
        prefix = remote_path(record)
        result = require_service().request(prefix + "/results")
        if result.get("run_id") != record["remote_id"]:
            raise ServiceError(502, "Result belongs to a different design job")

        def local_url(value):
            if not isinstance(value, str) or not value.startswith(prefix + "/artifacts/"):
                raise ServiceError(502, "Unexpected structure artifact location")
            suffix = value[len(prefix + "/artifacts/") :]
            return f"/api/jobs/{job_id}/artifacts/" + suffix

        if not isinstance(result.get("candidates"), list):
            raise ServiceError(502, "Design service returned invalid candidates")
        for candidate in result["candidates"]:
            if not isinstance(candidate, dict):
                raise ServiceError(502, "Design service returned an invalid candidate")
            for field in ("binder_structure_url", "complex_structure_url"):
                candidate[field] = local_url(candidate.get(field))
        if result.get("bundle_url"):
            result["bundle_url"] = local_url(result["bundle_url"])
        result["download_url"] = f"/api/jobs/{job_id}/download"
        result["target_id"] = record["target_id"]
        return result

    @app.get("/api/jobs/{job_id}/artifacts/{artifact_path:path}")
    def artifact(job_id: str, artifact_path: str):
        if (
            any(part in {"", ".", ".."} for part in artifact_path.split("/"))
            or "\\" in artifact_path
        ):
            raise HTTPException(400, "Invalid artifact path")
        _, record = receipt(job_id)
        data = require_service().request(
            remote_path(record) + "/artifacts/" + quote(artifact_path, safe="/"),
            binary=True,
        )
        return Response(data, media_type="text/plain")

    @app.get("/api/jobs/{job_id}/download")
    def download(job_id: str):
        _, record = receipt(job_id)
        data = require_service().request(remote_path(record) + "/download", binary=True)
        return Response(
            data,
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="{job_id}.zip"'},
        )

    app.mount("/static", StaticFiles(directory=STATIC), name="static")

    @app.get("/", include_in_schema=False)
    def index_page():
        return FileResponse(STATIC / "index.html")

    return app


def from_environment():
    base = os.environ.get("BIOCAD_BASE_URL")
    token_file = os.environ.get("BIOCAD_API_TOKEN_FILE")
    token = (
        Path(token_file).read_text().strip() if token_file else os.environ.get("BIOCAD_API_TOKEN")
    )
    return create_app(
        os.environ.get("BESPOKE_DEMO_INDEX"),
        os.environ.get("BESPOKE_DEMO_RUNS"),
        DesignService(base, token) if base else None,
        allow_run=os.environ.get("BESPOKE_DEMO_ALLOW_RUN") == "1",
    )


app = from_environment()
