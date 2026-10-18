"""Local teammate API: prepare new designs and replay actual cached demo artifacts.

Run: uvicorn biocad.service:app --host 127.0.0.1 --port 8000
Compute must already be configured; only the explicit /run endpoint starts inference.
"""

import io
import json
import os
import shutil
import threading
import uuid
import zipfile
from pathlib import Path
from typing import Literal

from fastapi import BackgroundTasks, FastAPI, HTTPException
from fastapi.responses import FileResponse, StreamingResponse
from pydantic import BaseModel, ConfigDict, Field, model_validator

from bespoke_bio.api import TranslateResponse
from bespoke_bio.api import app as upstream_app

from .boltzgen_adapter import prepare_target
from .design_bundle import export_bundle
from .pipeline import save
from .upstream import contained_path, import_targets

REPO = Path(__file__).resolve().parents[2]


class DesignInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    sequence: str | None = Field(default=None, max_length=10000)
    structure_text: str | None = Field(default=None, max_length=20_000_000)
    format: Literal["pdb", "cif", "mmcif"] = "cif"
    chain: str | None = None
    mode: Literal["prepare"] = "prepare"
    binder_length: int = Field(default=75, ge=25, le=300)
    num_designs: int = Field(default=2, ge=1, le=10000)

    @model_validator(mode="after")
    def one_input(self):
        if (self.sequence is None) == (self.structure_text is None):
            raise ValueError("Provide exactly one of sequence or structure_text")
        return self


class UpstreamInput(TranslateResponse):
    artifact_texts: dict[str, str] = Field(default_factory=dict)
    provenance: dict = Field(default_factory=dict)
    chain_selections: dict[str, list[str]] = Field(default_factory=dict)
    binder_length: int = Field(default=75, ge=25, le=300)
    num_designs: int = Field(default=2, ge=1, le=10000)


class RunInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    target_index: int | None = Field(default=None, ge=0)


def create_app(runs_dir=None, *, runner=None):
    app = FastAPI(
        title="Bespoke Bio design pipeline",
        version="0.1.0",
        description="Prepare BoltzGen jobs from sequence, structure, or upstream targets. Demo replays actual cached outputs. GPU execution requires an explicit /run request and a configured existing worker.",
    )
    app.include_router(upstream_app.router)
    claim_lock = threading.Lock()
    worker_lock = threading.Lock()

    def root():
        path = Path(
            runs_dir or os.environ.get("BIOCAD_RUNS_DIR", REPO / "outputs/service")
        ).resolve()
        path.mkdir(parents=True, exist_ok=True)
        return path

    def allocate():
        run_id = uuid.uuid4().hex
        folder = root() / run_id
        folder.mkdir()
        return run_id, folder

    def locate(run_id):
        if len(run_id) != 32 or any(c not in "0123456789abcdef" for c in run_id):
            raise HTTPException(404, "Unknown run")
        folder = root() / run_id
        if not folder.is_dir():
            raise HTTPException(404, "Unknown run")
        return folder

    def finish(run_id, folder, payload):
        response = {
            "run_id": run_id,
            "cloud_launched": False,
            **payload,
            "status_url": f"/v1/designs/{run_id}",
            "artifact_base_url": f"/v1/designs/{run_id}/artifacts/",
        }
        temporary = folder / ".status.tmp.json"
        save(temporary, response)
        temporary.replace(folder / "status.json")
        return response

    @app.post("/v1/designs", status_code=201)
    def prepare_design(request: DesignInput):
        run_id, folder = allocate()
        source = folder / (
            "input.fasta" if request.sequence is not None else "input." + request.format
        )
        source.write_text(
            ">target\n" + request.sequence.strip() + "\n"
            if request.sequence is not None
            else request.structure_text
        )
        try:
            plan = prepare_target(
                source,
                folder / "prepared",
                chain=request.chain,
                binder_length=request.binder_length,
                num_designs=request.num_designs,
            )
        except (ValueError, RuntimeError, OSError) as exc:
            finish(run_id, folder, {"status": "failed_preparation", "error": str(exc)})
            raise HTTPException(422, str(exc)) from exc
        return finish(
            run_id,
            folder,
            {
                "status": "prepared_pending_gpu_execution",
                "mode": "prepare",
                "plan": plan,
                "prepared_directory": "prepared",
                "plan_path": "prepared/boltzgen_plan.json",
                "target_state_path": "prepared/target_state.json",
                "design_yaml": "prepared/" + plan["command"][2],
                "structural_validation": "pending_model_execution"
                if request.sequence is not None
                else "input_coordinates_prepared",
                "execution_performed": False,
            },
        )

    @app.post("/v1/designs/demo", status_code=201)
    def demo_design():
        run_id, folder = allocate()
        plan = prepare_target(
            REPO / "examples/5JDS.cif",
            folder / "prepared",
            chain="A",
            name="pdl1",
            evidence="experimental",
        )
        bundle = export_bundle(
            REPO / "examples/boltzgen_results",
            folder / "prepared/target_state.json",
            folder / "bundle",
        )
        return finish(
            run_id,
            folder,
            {
                "status": "cached_demo_completed",
                "mode": "actual_cached_replay",
                "execution_performed": False,
                "target": "PD-L1 observed chain A from 5JDS",
                "plan": plan,
                "prepared_directory": "prepared",
                "plan_path": "prepared/boltzgen_plan.json",
                "target_state_path": "prepared/target_state.json",
                "design_yaml": "prepared/" + plan["command"][2],
                "bundle": "bundle/bundle.json",
                "candidates": len(bundle["components"]),
                "native_filter_passed": sum(
                    e["status"] == "pass"
                    for e in bundle["evaluations"]
                    if e["id"].endswith("_filters")
                ),
                "binding_evidence": "computational only; no measured binding or affinity",
            },
        )

    @app.post("/v1/designs/upstream", status_code=201)
    def upstream_design(request: UpstreamInput):
        run_id, folder = allocate()
        package = folder / "upstream_package"
        package.mkdir()
        save(package / "provenance.json", {"source": "teammate HTTP intake", **request.provenance})
        (package / "targets.jsonl").write_text(
            "".join(p.model_dump_json() + "\n" for p in request.targets)
        )
        try:
            for relative, text in request.artifact_texts.items():
                destination = contained_path(package, relative)
                if destination == package or relative in {"provenance.json", "targets.jsonl"}:
                    raise ValueError("Reserved artifact path")
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(text)
            imported = import_targets(package, folder / "imported", request.chain_selections)
            proteins = {p.target_id: p for p in request.targets}
            states = {s.id: s for s in imported.target_states}
            records = []
            for index, record in enumerate(imported.records):
                result = dict(record)
                if record["status"] == "failed":
                    records.append(result)
                    continue
                prepared = folder / "prepared" / str(index)
                if record["status"] == "imported":
                    state = states[record["target_state_id"]]
                    source = folder / "imported" / state.prepared
                    chain = record["author_chain"]
                else:
                    source = folder / f"unresolved_{index}.fasta"
                    source.write_text(
                        ">target\n"
                        + proteins[record["upstream_target_id"]].amino_acid_sequence
                        + "\n"
                    )
                    chain = None
                try:
                    plan = prepare_target(
                        source,
                        prepared,
                        chain=chain,
                        binder_length=request.binder_length,
                        num_designs=request.num_designs,
                    )
                    prepared_state_path = prepared / "target_state.json"
                    prepared_state = json.loads(prepared_state_path.read_text())
                    protein = proteins[record["upstream_target_id"]]
                    prepared_state["context"]["upstream_record"] = protein.model_dump(mode="json")
                    prepared_state["context"]["upstream_provenance"] = imported.provenance
                    if record["status"] == "imported":
                        prepared_state["context"]["upstream_imported_state"] = state.model_dump()
                        prepared_state["context"]["upstream_target_state_id"] = state.id
                        prepared_state["measurements"] = state.measurements
                        prepared_state["variants"] = state.variants
                    else:
                        prepared_state["variants"] = protein.source_variant_ids
                        prepared_state["measurements"] = [
                            {
                                "kind": "measured_RNA",
                                "units": "TPM",
                                "value": protein.expression.tpm,
                                "raw_count": protein.expression.raw_count,
                                "transcript_id": protein.transcript_id,
                                "protein_abundance": None,
                            }
                        ]
                    save(prepared_state_path, prepared_state)
                    result.update(
                        status="prepared_pending_gpu_execution",
                        plan_path=str((prepared / "boltzgen_plan.json").relative_to(folder)),
                        prepared_directory=str(prepared.relative_to(folder)),
                        target_state_path=str((prepared / "target_state.json").relative_to(folder)),
                        design_yaml=str((prepared / plan["command"][2]).relative_to(folder)),
                        target_input_mode=plan["target_input_mode"],
                    )
                except (ValueError, RuntimeError, OSError) as exc:
                    result.update(status="failed_preparation", reason=str(exc))
                if result["status"] == "prepared_pending_gpu_execution":
                    child_id, child_folder = allocate()
                    shutil.copytree(prepared, child_folder / "prepared")
                    finish(
                        child_id,
                        child_folder,
                        {
                            "status": "prepared_pending_gpu_execution",
                            "mode": "prepare",
                            "execution_performed": False,
                            "parent_run_id": run_id,
                            "upstream_target_id": record["upstream_target_id"],
                            "plan": plan,
                            "prepared_directory": "prepared",
                            "plan_path": "prepared/boltzgen_plan.json",
                            "target_state_path": "prepared/target_state.json",
                            "design_yaml": "prepared/" + plan["command"][2],
                            "structural_validation": "pending_model_execution"
                            if plan["target_input_mode"] == "sequence_only"
                            else "input_coordinates_prepared",
                        },
                    )
                    result.update(
                        child_run_id=child_id,
                        child_status_url=f"/v1/designs/{child_id}",
                        child_run_url=f"/v1/designs/{child_id}/run",
                        child_download_url=f"/v1/designs/{child_id}/download",
                    )
                records.append(result)
        except (ValueError, RuntimeError, OSError) as exc:
            finish(run_id, folder, {"status": "failed_intake", "error": str(exc)})
            raise HTTPException(422, str(exc)) from exc
        ready = sum(r["status"] == "prepared_pending_gpu_execution" for r in records)
        return finish(
            run_id,
            folder,
            {
                "status": "prepared_pending_gpu_execution" if ready else "failed_preparation",
                "mode": "upstream_prepare",
                "execution_performed": False,
                "prepared_targets": ready,
                "records": records,
                "imported_targets": "imported/targets.json",
            },
        )

    def execute_live(run_id, folder, original, prepared, configuration):
        with worker_lock:
            active = {
                **original,
                "status": "running",
                "mode": "live_model_outputs",
                "execution_performed": True,
                "execution_attempted": True,
            }
            finish(run_id, folder, active)
            try:
                selected_runner = runner
                if selected_runner is None:
                    from .worker import run_job

                    selected_runner = run_job
                receipt = selected_runner(prepared, folder / "live", **configuration, timeout=3600)
                completed = receipt.get("status") == "completed"
                payload = {
                    **active,
                    "status": "completed" if completed else "failed",
                    "worker_receipt": receipt,
                    "cloud_launched": False,
                }
                if completed:
                    bundle_path = contained_path(
                        folder / "live", receipt.get("bundle", "bundle/bundle.json")
                    )
                    if not bundle_path.is_file():
                        raise FileNotFoundError(
                            "Worker reported completion without a downloadable bundle"
                        )
                    payload["bundle"] = str(bundle_path.relative_to(folder))
                finish(run_id, folder, payload)
            except Exception as exc:  # noqa: BLE001 - background failures must persist a terminal job state
                finish(
                    run_id,
                    folder,
                    {
                        **active,
                        "status": "failed",
                        "error": str(exc),
                        "remote_cleanup_required": True,
                    },
                )

    @app.post("/v1/designs/{run_id}/run", status_code=202)
    def run_design(run_id: str, background_tasks: BackgroundTasks, request: RunInput | None = None):
        """Explicitly execute against an operator-configured existing worker; never provision a pod."""
        request = request or RunInput()
        folder = locate(run_id)
        host = os.environ.get("BIOCAD_WORKER_HOST")
        key = os.environ.get("BIOCAD_WORKER_KEY")
        port = os.environ.get("BIOCAD_WORKER_PORT")
        if not all((host, key, port)):
            raise HTTPException(
                503,
                "Live execution disabled: configure BIOCAD_WORKER_HOST, BIOCAD_WORKER_PORT and BIOCAD_WORKER_KEY",
            )
        try:
            port_number = int(port)
            if not 1 <= port_number <= 65535 or not Path(key).is_file():
                raise ValueError("Invalid port or missing SSH key file")
        except ValueError as exc:
            raise HTTPException(503, "Live worker port/key configuration is invalid") from exc
        with claim_lock:
            status_path = folder / "status.json"
            if not status_path.is_file():
                raise HTTPException(409, "Run has no prepared job")
            current = json.loads(status_path.read_text())
            if current.get("status") != "prepared_pending_gpu_execution":
                raise HTTPException(409, "Job is already submitted or is not ready for execution")
            if current.get("mode") == "upstream_prepare":
                raise HTTPException(
                    422,
                    {
                        "message": "Execute each prepared upstream child run independently using child_run_url",
                        "child_runs": [
                            r.get("child_run_url")
                            for r in current.get("records", [])
                            if r.get("child_run_url")
                        ],
                    },
                )
            if request.target_index is not None:
                raise HTTPException(
                    422, "Use the selected upstream record's child_run_url instead of target_index"
                )
            selected = current
            prepared = contained_path(folder, selected["prepared_directory"])
            plan = json.loads((prepared / "boltzgen_plan.json").read_text())
            if not 1 <= plan.get("num_designs", 0) <= 32:
                raise HTTPException(422, "Live API jobs are limited to 1..32 designs")
            queued = finish(
                run_id,
                folder,
                {
                    **current,
                    "status": "queued",
                    "mode": "live_model_outputs",
                    "execution_performed": False,
                    "selected_target_index": request.target_index,
                    "worker_scope": "existing operator-configured worker; one job at a time per API process",
                },
            )
            background_tasks.add_task(
                execute_live,
                run_id,
                folder,
                queued,
                prepared,
                {"host": host, "port": port_number, "key": key},
            )
        return queued

    @app.get("/v1/designs/{run_id}")
    def status(run_id: str):
        folder = locate(run_id)
        path = folder / "status.json"
        if not path.is_file():
            raise HTTPException(404, "Run status unavailable")
        payload = json.loads(path.read_text())
        if payload.get("mode") == "upstream_prepare":
            child_statuses = []
            for record in payload.get("records", []):
                child_id = record.get("child_run_id")
                if child_id:
                    child = json.loads((locate(child_id) / "status.json").read_text())
                    record["status"] = child["status"]
                    record["execution_performed"] = child.get("execution_performed", False)
                    child_statuses.append(child["status"])
            if child_statuses and all(s == "completed" for s in child_statuses):
                payload["status"] = "completed"
            elif any(s == "running" for s in child_statuses):
                payload["status"] = "running"
            elif any(s == "queued" for s in child_statuses):
                payload["status"] = "queued"
            elif child_statuses and all(s in {"completed", "failed"} for s in child_statuses):
                payload["status"] = "completed_with_failures"
            payload["execution_performed"] = any(
                r.get("execution_performed", False) for r in payload.get("records", [])
            )
        payload["artifacts"] = [
            str(p.relative_to(folder)) for p in sorted(folder.rglob("*")) if p.is_file()
        ]
        return payload

    @app.get("/v1/designs/{run_id}/download")
    def download(run_id: str):
        """Download only this run's contained regular artifacts as one ZIP."""
        folder = locate(run_id)
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as output:
            for path in sorted(folder.rglob("*")):
                if path.is_symlink() or not path.is_file():
                    continue
                resolved = path.resolve()
                if not resolved.is_relative_to(folder.resolve()):
                    continue
                output.write(resolved, arcname=str(path.relative_to(folder)))
        archive.seek(0)
        return StreamingResponse(
            archive,
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="{run_id}.zip"'},
        )

    @app.get("/v1/designs/{run_id}/artifacts/{artifact_path:path}")
    def artifact(run_id: str, artifact_path: str):
        folder = locate(run_id)
        try:
            path = contained_path(folder, artifact_path)
        except ValueError as exc:
            raise HTTPException(404, "Artifact unavailable") from exc
        if not path.is_file():
            raise HTTPException(404, "Artifact unavailable")
        return FileResponse(path, filename=path.name)

    return app


app = create_app()
