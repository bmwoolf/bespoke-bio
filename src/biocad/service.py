"""Local teammate API: prepare new designs and replay actual cached demo artifacts.

Run: uvicorn biocad.service:app --host 127.0.0.1 --port 8000
Only the explicit /run endpoint starts inference, using a configured static or managed worker.
"""

import asyncio
import io
import json
import os
import re
import secrets
import shutil
import threading
import time
import uuid
import zipfile
from contextlib import asynccontextmanager
from hashlib import sha256
from pathlib import Path
from typing import Annotated, Literal
from urllib.parse import quote

from fastapi import BackgroundTasks, Depends, FastAPI, HTTPException, Request, Security
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel, ConfigDict, Field, model_validator

from bespoke_bio.api import TranslateResponse
from bespoke_bio.api import app as upstream_app

from .boltzgen_adapter import prepare_target
from .design_bundle import DesignBundle, export_bundle
from .jobs import BudgetExceeded, JobConflict, JobStore, ManagedQueue
from .pipeline import save
from .upstream import contained_path, import_targets

REPO = Path(os.environ.get("BIOCAD_REPO_DIR", Path(__file__).resolve().parents[2])).resolve()


class DesignInput(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
        json_schema_extra={
            "examples": [
                {
                    "sequence": "AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNA",
                    "binder_length": 75,
                    "num_designs": 2,
                    "mode": "prepare",
                }
            ]
        },
    )
    sequence: str | None = Field(
        default=None,
        max_length=10000,
        description="Target protein amino-acid sequence, using canonical one-letter codes. Supply sequence OR structure_text. Sequence input has no prior folded coordinates; its structural prediction is pending model execution.",
    )
    amino_acid_sequence: str | None = Field(
        default=None,
        max_length=10000,
        description="Upstream alias for sequence. Submitted residues are preserved; no automatic signal-peptide trimming.",
    )
    target_id: str | None = None
    donor: str | None = None
    gene_id: str | None = None
    gene_symbol: str | None = None
    transcript_id: str | None = None
    phase_scenario: str | None = None
    sequence_sha256: str | None = Field(
        default=None,
        description="Optional SHA-256 of the exact submitted sequence, bare hex or sha256: prefixed. Mismatches are rejected.",
    )
    reference_sequence_sha256: str | None = None
    expression: dict | float | None = None
    structure_text: str | None = Field(
        default=None,
        max_length=20_000_000,
        description="Complete PDB or mmCIF file contents as a JSON string, not a filename or URL. Supply this instead of sequence. The original file and author-to-normalized residue mapping are retained.",
    )
    format: Literal["pdb", "cif", "mmcif"] = Field(
        default="cif", description="Format of structure_text; ignored for sequence input."
    )
    chain: str | None = Field(
        default=None,
        description="Author chain ID to select from structure_text, e.g. A. Required for structures containing multiple protein chains. Omit for a single protein chain or sequence input.",
    )
    mode: Literal["prepare"] = Field(
        default="prepare",
        description="Preparation creates downloadable BoltzGen inputs. It does not start inference; explicitly POST the returned status_url + /run to use a configured GPU worker.",
    )
    binder_length: int = Field(
        default=75,
        ge=25,
        le=300,
        description="Number of amino-acid residues to generate in each binder. Start with 75. No epitope or hotspot is specified.",
    )
    num_designs: int = Field(
        default=2,
        ge=1,
        le=10000,
        description="Number of candidate designs to request. Start with 2 for a smoke test. Live HTTP execution is limited to 32 designs per job; generation does not guarantee candidates pass filters or bind experimentally.",
    )

    @model_validator(mode="after")
    def one_input(self):
        if self.amino_acid_sequence is not None:
            if self.sequence is not None and self.sequence != self.amino_acid_sequence:
                raise ValueError("sequence and amino_acid_sequence disagree")
            self.sequence = self.amino_acid_sequence
        if self.sequence_sha256 is not None and (
            self.sequence is None
            or sha256(self.sequence.encode()).hexdigest()
            != self.sequence_sha256.removeprefix("sha256:").lower()
        ):
            raise ValueError("sequence_sha256 does not match the submitted sequence")
        if (self.sequence is None) == (self.structure_text is None):
            raise ValueError("Provide exactly one of sequence or structure_text")
        return self


class UpstreamInput(TranslateResponse):
    artifact_texts: dict[str, str] = Field(default_factory=dict)
    provenance: dict = Field(default_factory=dict)
    chain_selections: dict[str, list[str]] = Field(default_factory=dict)
    binder_length: int = Field(default=75, ge=25, le=300)
    num_designs: int = Field(default=2, ge=1, le=10000)


class JobResponse(BaseModel):
    """Common job envelope; extra stage-specific fields are preserved."""

    model_config = ConfigDict(
        extra="allow",
        json_schema_extra={
            "examples": [
                {
                    "run_id": "0123456789abcdef0123456789abcdef",
                    "status": "prepared_pending_gpu_execution",
                    "mode": "prepare",
                    "execution_performed": False,
                    "cloud_launched": False,
                    "status_url": "/v1/designs/0123456789abcdef0123456789abcdef",
                    "artifact_base_url": "/v1/designs/0123456789abcdef0123456789abcdef/artifacts/",
                    "prepared_directory": "prepared",
                    "design_yaml": "prepared/design.yaml",
                }
            ]
        },
    )
    run_id: str = Field(
        description="Unique run identifier. Upstream batches also return independent child runs for each prepared target."
    )
    status: str = Field(
        description="Job progress: prepared_pending_gpu_execution, queued, running, completed, failed, or cached_demo_completed. Parent batches aggregate child progress."
    )
    mode: str | None = Field(
        default=None,
        description="prepare, upstream_prepare, actual_cached_replay, or live_model_outputs. Cached replay is never a new GPU experiment.",
    )
    execution_performed: bool = Field(
        default=False,
        description="Whether the worker reports actual remote model invocation. A queued job, cached demo, or failed connection is not model execution.",
    )
    execution_attempted: bool = Field(
        default=False,
        description="Whether live dispatch was attempted. This does not guarantee successful connection, inference, or binding.",
    )
    cloud_launched: bool = Field(
        default=False,
        description="False for preparation, cached replay and static workers. Managed execution reports whether its receipt contains a newly provisioned pod ID.",
    )
    status_url: str = Field(
        description="GET this relative URL to poll progress and list downloadable artifacts."
    )
    artifact_base_url: str = Field(
        description="Append a returned artifact path to download one file. GET status_url + /download downloads this run as a ZIP."
    )
    bundle: str | None = Field(
        default=None,
        description="Relative artifact path to the typed binder bundle when available. Contains TargetState, Component and Evaluation records, sequence/structure links and provenance.",
    )
    candidates: int | None = Field(
        default=None,
        description="Number of designed candidates exported. Does not imply experimental binding.",
    )
    native_filter_passed: int | None = Field(
        default=None,
        description="Candidates passing native computational filters; separate from interface confidence or measured affinity.",
    )
    plan: dict | None = Field(
        default=None,
        description="BoltzGen inputs, command, sequence hashes, unconstrained site policy and execution stages.",
    )
    records: list[dict] | None = Field(
        default=None,
        description="Upstream target results, including failures and independent child_run_url/child_status_url/child_download_url for prepared targets.",
    )
    artifacts: list[str] | None = Field(
        default=None,
        description="Relative file paths available via artifact_base_url. Returned by the status endpoint.",
    )
    progress: dict | None = Field(
        default=None,
        description="Latest native pipeline stage observed in the last 1 MiB of the worker log: stage_index, total_stages, stage. Stage entry does not imply stage or pipeline completion.",
    )
    prepared_directory: str | None = Field(
        default=None,
        description="Run-relative directory containing the complete prepared model inputs.",
    )
    plan_path: str | None = None
    target_state_path: str | None = None
    design_yaml: str | None = Field(
        default=None,
        description="Run-relative BoltzGen YAML input. Target is fixed-sequence or fixed-structure; binder length is generated, with no specified epitope.",
    )


class CandidateResult(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    sequence: str
    target_state_id: str | None
    native_filter_passed: bool | None = Field(
        description="Native computational filter result, not a binding claim."
    )
    failed_filters: list[str]
    upstream_filter_flags: dict[str, str]
    native_filter_metrics: dict[str, float | int | None]
    native_filter_metric_units: dict[str, str]
    confidence_metrics: dict[str, float | int | None]
    confidence_units: dict[str, str]
    confidence_evidence: str = Field(
        description="Origin and limits of the predicted interface confidence."
    )
    binder_structure_url: str
    complex_structure_url: str
    evidence: dict


class ResultsResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")
    run_id: str
    status: str
    mode: str | None
    target_state_ids: list[str]
    candidate_count: int
    candidates: list[CandidateResult]
    bundle_url: str | None
    download_url: str
    records: list[dict] | None = Field(
        default=None,
        description="Upstream parent results are an index of independent child jobs, not flattened binder results. Follow child_results_url after each child completes.",
    )
    evidence: str


class RunInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    target_index: int | None = Field(
        default=None,
        ge=0,
        json_schema_extra={"deprecated": True},
        description="Use each upstream record child_run_url instead. Direct jobs accept an empty object or no body.",
    )


def create_app(runs_dir=None, *, runner=None, managed_runner=None, managed_cleanup=None):
    managed = os.environ.get("BIOCAD_MANAGED_WORKERS") == "1"
    bearer = HTTPBearer(auto_error=False, scheme_name="DesignServiceToken")

    def authenticate(
        request: Request,
        credentials: Annotated[HTTPAuthorizationCredentials | None, Security(bearer)] = None,
    ):
        if request.url.path == "/healthz":
            return
        token = os.environ.get("BIOCAD_API_TOKEN")
        if managed and not token:
            raise HTTPException(503, "Managed service requires BIOCAD_API_TOKEN")
        if token and (
            credentials is None or not secrets.compare_digest(credentials.credentials, token)
        ):
            raise HTTPException(
                401,
                "Provide the service bearer token using Authorize",
                headers={"WWW-Authenticate": "Bearer"},
            )

    app = FastAPI(
        title="Bespoke Bio design pipeline",
        dependencies=[Depends(authenticate)],
        version="0.1.0",
        description="Prepare BoltzGen jobs from sequence, structure, or upstream targets. Demo replays actual cached outputs. GPU execution requires an explicit /run request. Managed mode reserves the phase budget and provisions a disposable worker; static mode uses a configured existing worker.",
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

    @app.post(
        "/v1/designs",
        status_code=201,
        response_model=JobResponse,
        response_model_exclude_unset=True,
    )
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
        upstream_fields = {
            "amino_acid_sequence",
            "target_id",
            "donor",
            "gene_id",
            "gene_symbol",
            "transcript_id",
            "phase_scenario",
            "sequence_sha256",
            "reference_sequence_sha256",
            "expression",
        }
        if request.model_fields_set & upstream_fields:
            record = request.model_dump(include=upstream_fields, exclude_unset=True)
            save(folder / "upstream_record.json", record)
            state_path = folder / "prepared/target_state.json"
            state = json.loads(state_path.read_text())
            state.setdefault("context", {})["upstream_record"] = record
            save(state_path, state)
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

    @app.post(
        "/v1/designs/demo",
        status_code=201,
        response_model=JobResponse,
        response_model_exclude_unset=True,
    )
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

    @app.post(
        "/v1/designs/upstream",
        status_code=201,
        response_model=JobResponse,
        response_model_exclude_unset=True,
    )
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
                "execution_performed": False,
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
                    "execution_performed": receipt.get("execution_performed", completed),
                    "execution_attempted": True,
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
                    exported = json.loads(bundle_path.read_text())
                    payload["candidates"] = len(exported.get("components", []))
                    payload["native_filter_passed"] = sum(
                        e["status"] == "pass"
                        for e in exported.get("evaluations", [])
                        if e["id"].endswith("_filters")
                    )
                    payload["structural_validation"] = "model_predictions_exported"
                    payload["binding_evidence"] = (
                        "computational only; no measured binding or affinity"
                    )
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

    @app.post(
        "/v1/designs/{run_id}/run",
        status_code=202,
        response_model=JobResponse,
        response_model_exclude_unset=True,
    )
    def run_design(run_id: str, background_tasks: BackgroundTasks, request: RunInput | None = None):
        """Explicitly submit to the configured static worker or durable managed GPU queue."""
        request = request or RunInput()
        folder = locate(run_id)
        host = os.environ.get("BIOCAD_WORKER_HOST")
        key = os.environ.get("BIOCAD_WORKER_KEY")
        port = os.environ.get("BIOCAD_WORKER_PORT")
        if not managed and not all((host, key, port)):
            raise HTTPException(
                503,
                "Live execution disabled: configure BIOCAD_WORKER_HOST, BIOCAD_WORKER_PORT and BIOCAD_WORKER_KEY",
            )
        port_number = None
        if managed:
            queue = getattr(app.state, "managed_queue", None)
            if queue is None or not queue.ready:
                raise HTTPException(
                    503, "Managed queue is not ready; recovery or cleanup is required"
                )
            if queue.expires_at and time.time() + 7200 > queue.expires_at:
                raise HTTPException(
                    503, "Shared service expires in less than two hours; new GPU jobs are disabled"
                )
        else:
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
            if managed:
                try:
                    queue.store.enqueue(run_id, prepared, folder / "live")
                except JobConflict as exc:
                    raise HTTPException(409, str(exc)) from exc
                except BudgetExceeded as exc:
                    raise HTTPException(402, str(exc)) from exc
                return finish(
                    run_id,
                    folder,
                    {
                        **current,
                        "status": "queued",
                        "mode": "live_model_outputs",
                        "execution_performed": False,
                        "execution_attempted": False,
                        "managed_worker": True,
                        "budget": queue.store.budget(),
                        "worker_scope": "durable managed queue; 10 USD reserved, maximum 2 hours at 4 USD/hour; pod deleted after result download",
                    },
                )
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

    @app.post(
        "/v1/designs/{run_id}/cancel", response_model=JobResponse, response_model_exclude_unset=True
    )
    def cancel_design(run_id: str):
        folder = locate(run_id)
        queue = getattr(app.state, "managed_queue", None)
        if not managed or queue is None:
            raise HTTPException(409, "Durable cancellation is available for managed jobs only")
        try:
            queue.cancel(run_id)
        except KeyError as exc:
            raise HTTPException(409, "Job has not been queued") from exc
        except JobConflict as exc:
            raise HTTPException(409, str(exc)) from exc
        return json.loads((folder / "status.json").read_text())

    @app.get("/v1/service", tags=["Service"])
    def service_status():
        queue = getattr(app.state, "managed_queue", None)
        return {
            "managed_workers": managed,
            "ready": bool(queue and queue.ready) if managed else True,
            "active_job_id": queue.active_id if queue else None,
            "budget": queue.store.budget() if queue else None,
            "expires_at": queue.expires_at if queue else None,
        }

    @app.get("/v1/designs/{run_id}", response_model=JobResponse, response_model_exclude_unset=True)
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
        log_path = folder / "live/worker.log"
        if log_path.is_file() and log_path.resolve().is_relative_to(folder.resolve()):
            with log_path.open("rb") as handle:
                handle.seek(0, 2)
                handle.seek(max(0, handle.tell() - 1024 * 1024))
                tail = handle.read(1024 * 1024).decode("utf-8", errors="replace")
            matches = re.findall(r"Pipeline step (\d+) of (\d+): ([A-Za-z0-9_]+)", tail)
            valid = [
                (int(index), int(total), stage)
                for index, total, stage in matches
                if 1 <= int(index) <= int(total)
            ]
            if valid:
                index, total, stage = valid[-1]
                payload["progress"] = {"stage_index": index, "total_stages": total, "stage": stage}
        payload["artifacts"] = [
            str(p.relative_to(folder)) for p in sorted(folder.rglob("*")) if p.is_file()
        ]
        return payload

    @app.get(
        "/v1/designs/{run_id}/results",
        response_model=ResultsResponse,
        response_model_exclude_none=True,
    )
    def results(run_id: str):
        """Normalized candidate handoff from the actual completed bundle, or an upstream child index."""
        folder = locate(run_id)
        current = status(run_id)
        base = {
            "run_id": run_id,
            "status": current["status"],
            "mode": current.get("mode"),
            "download_url": f"/v1/designs/{run_id}/download",
            "evidence": "Computational designs and predictions only; no measured binding or affinity.",
        }
        if current.get("mode") == "upstream_prepare":
            records = []
            targets = []
            for record in current.get("records", []):
                entry = dict(record)
                child_id = entry.get("child_run_id")
                if child_id:
                    entry["child_results_url"] = f"/v1/designs/{child_id}/results"
                    child = json.loads((locate(child_id) / "status.json").read_text())
                    target_id = child.get("plan", {}).get("target_state_id")
                    if target_id:
                        targets.append(target_id)
                records.append(entry)
            return {
                **base,
                "target_state_ids": sorted(set(targets)),
                "candidate_count": 0,
                "candidates": [],
                "bundle_url": None,
                "records": records,
            }
        if current["status"] not in {"completed", "cached_demo_completed"} or not current.get(
            "bundle"
        ):
            raise HTTPException(
                409,
                {
                    "message": "Candidate results are not available until the run completes. Poll status_url; no cached candidates are substituted.",
                    "status": current["status"],
                    "status_url": current["status_url"],
                },
            )
        bundle_path = contained_path(folder, current["bundle"])
        if not bundle_path.is_file():
            raise HTTPException(409, "Completed job has no downloadable bundle")
        bundle = DesignBundle.model_validate_json(bundle_path.read_text())

        def artifact_url(relative):
            artifact_path = contained_path(bundle_path.parent, relative)
            if not artifact_path.is_file():
                raise HTTPException(409, "A candidate structure artifact is unavailable")
            return current["artifact_base_url"] + quote(
                str(artifact_path.relative_to(folder)), safe="/"
            )

        candidates = []
        for component in bundle.components:
            evaluations = [e for e in bundle.evaluations if e.subject_id == component.id]
            native = next((e for e in evaluations if e.id.endswith("_filters")), None)
            confidence = next((e for e in evaluations if e.id.endswith("_binding")), None)
            metrics = confidence.metrics if confidence else {}
            candidates.append(
                {
                    "id": component.id,
                    "sequence": component.sequence,
                    "target_state_id": component.evidence.get("target_state_id"),
                    "native_filter_passed": native.status == "pass" if native else None,
                    "failed_filters": native.details.get("failed_filters", []) if native else [],
                    "upstream_filter_flags": native.details.get("flags", {}) if native else {},
                    "native_filter_metrics": native.metrics if native else {},
                    "native_filter_metric_units": {
                        metric: native.units.get(
                            metric,
                            "angstrom"
                            if "rmsd" in metric
                            else "dimensionless"
                            if metric
                            in {
                                "design_ptm",
                                "design_to_target_iptm",
                                "esmfold2_ipsae_min",
                                "quality_score",
                            }
                            else "unspecified",
                        )
                        for metric in (native.metrics if native else {})
                    },
                    "confidence_metrics": metrics,
                    "confidence_units": {
                        metric: confidence.units.get(metric, "dimensionless") for metric in metrics
                    },
                    "confidence_evidence": confidence.evidence
                    if confidence
                    else "Interface confidence missing; not assessed.",
                    "binder_structure_url": artifact_url(component.structure),
                    "complex_structure_url": artifact_url(component.evidence["complex_structure"]),
                    "evidence": component.evidence,
                }
            )
        return {
            **base,
            "target_state_ids": [target.id for target in bundle.target_states],
            "candidate_count": len(candidates),
            "candidates": candidates,
            "bundle_url": current["artifact_base_url"] + quote(current["bundle"], safe="/"),
        }

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

    def persist_managed(job, state, receipt):
        folder = Path(job["output"]).parent
        with claim_lock:
            latest = app.state.managed_queue.store.get(job["id"])
            state = latest["status"]
            receipt = latest.get("receipt") or receipt
            current = json.loads((folder / "status.json").read_text())
            payload = {
                **current,
                "status": state,
                "mode": "live_model_outputs",
                "managed_worker": True,
                "execution_attempted": state not in {"queued", "cancelled"}
                or current.get("execution_attempted", False),
                "execution_performed": receipt.get(
                    "execution_performed", current.get("execution_performed", False)
                ),
                "cloud_launched": bool(receipt.get("pod_id")),
                "budget": app.state.managed_queue.store.budget(),
            }
            if receipt:
                payload["worker_receipt"] = receipt
            if state == "completed":
                bundle_path = contained_path(
                    Path(job["output"]), receipt.get("bundle", "bundle/bundle.json")
                )
                if bundle_path.is_file():
                    payload["bundle"] = str(bundle_path.relative_to(folder))
                    exported = json.loads(bundle_path.read_text())
                    payload["candidates"] = len(exported.get("components", []))
                    payload["native_filter_passed"] = sum(
                        e["status"] == "pass"
                        for e in exported.get("evaluations", [])
                        if e["id"].endswith("_filters")
                    )
                    payload["structural_validation"] = "model_predictions_exported"
                else:
                    payload.update(
                        status="failed",
                        error="Managed worker completed without a downloaded bundle",
                    )
                    app.state.managed_queue.store.finish(job["id"], "failed", payload)
            finish(job["id"], folder, payload)

    @asynccontextmanager
    async def lifespan(_app):
        if managed:
            if not os.environ.get("BIOCAD_API_TOKEN"):
                raise RuntimeError("BIOCAD_API_TOKEN is required for managed shared service")
            selected_runner, selected_cleanup = managed_runner, managed_cleanup
            startup_reconcile = None
            if selected_runner is None:
                from .cloud import reconcile_owned_leases

                startup_reconcile = reconcile_owned_leases
            if selected_runner is None or selected_cleanup is None:
                from .cloud import cleanup_owned_pod, run_managed_job

                selected_runner = selected_runner or run_managed_job
                selected_cleanup = selected_cleanup or cleanup_owned_pod
            queue = ManagedQueue(
                JobStore(Path(os.environ.get("BIOCAD_JOB_DB", root() / "jobs.sqlite3"))),
                selected_runner,
                selected_cleanup,
                persist_managed,
                startup_reconcile=startup_reconcile,
                expires_at=float(os.environ["BIOCAD_SERVICE_EXPIRES_AT"])
                if os.environ.get("BIOCAD_SERVICE_EXPIRES_AT")
                else None,
            )
            app.state.managed_queue = queue
            await asyncio.to_thread(queue.start)
        try:
            yield
        finally:
            if managed:
                await asyncio.to_thread(app.state.managed_queue.stop)

    app.router.lifespan_context = lifespan
    return app


app = create_app()
