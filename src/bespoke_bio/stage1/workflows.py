"""Prepare reproducible native workflows without embedding machine or cloud choices."""

import csv
import io
import os
import re
import subprocess
from pathlib import Path

from pydantic import Field

from .contracts import Intake, Record
from .io import canonical, read_json, write_new
from .runtime import Runtime


class ReferenceAssets(Record):
    """Exact, prebuilt assets shared by both workflows; no implicit reference downloads."""

    fasta: Path
    dictionary: Path
    gtf: Path
    transcript_fasta: Path
    bwa_index: Path
    star_index: Path
    vep_cache: Path
    bundle_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


def sample_sheet(intake: Intake, assay: str) -> bytes:
    """One row per delivery pair. Distinct DNA specimens and RNA libraries remain distinct."""
    output = io.StringIO(newline="")
    writer = csv.writer(output)
    dna = assay == "germline_wgs"
    writer.writerow(
        ["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"]
        if dna
        else ["sample", "fastq_1", "fastq_2", "strandedness"]
    )
    reads = [r for r in intake.reads if r.assay == assay]
    for library, lane in sorted({(r.library, r.lane) for r in reads}):
        pair = sorted(
            [r for r in reads if (r.library, r.lane) == (library, lane)], key=lambda r: r.mate
        )
        if len(pair) != 2 or [r.mate for r in pair] != [1, 2]:
            raise ValueError("invalid workflow mate pair")
        first = pair[0]
        paths = [r.artifact.location for r in pair]
        writer.writerow(
            [first.donor, intake.sex, 0, first.specimen, lane, *paths]
            if dna
            else [library, *paths, "auto"]
        )
    return output.getvalue().encode()


def parameters(branch: str, sheet: Path, output: Path, assets: ReferenceAssets) -> dict:
    """The qualified biological baseline. Changing settings creates a different run."""
    common = {
        "input": str(sheet),
        "outdir": str(output),
        "fasta": str(assets.fasta),
        "igenomes_ignore": True,
        "custom_config_base": "",
        "save_reference": True,
    }
    if branch == "dna":
        return common | {
            "tools": "deepvariant,vep",
            "aligner": "bwa-mem2",
            "wes": False,
            "trim_fastq": False,
            "split_fastq": 0,
            "skip_tools": "baserecalibrator",
            "normalize_vcfs": True,
            "filter_vcfs": False,
            "dict": str(assets.dictionary),
            "fasta_fai": str(assets.fasta) + ".fai",
            "bwamem2": str(assets.bwa_index),
            "save_mapped": True,
            "save_output_as_bam": True,
            "vep_cache": str(assets.vep_cache),
            "vep_species": "homo_sapiens",
            "vep_genome": "GRCh38",
            "vep_cache_version": "116",
            "vep_include_fasta": True,
            "vep_custom_args": (
                "--everything --total_length --offline --format vcf "
                "--transcript_version --allele_number --check_ref"
            ),
        }
    if branch == "rna":
        return common | {
            "aligner": "star_salmon",
            "gencode": True,
            "gtf": str(assets.gtf),
            "transcript_fasta": str(assets.transcript_fasta),
            "star_index": str(assets.star_index),
            "trimmer": "trimgalore",
            "skip_trimming": False,
            "save_align_intermeds": True,
            "extra_salmon_quant_args": "",
            "stranded_threshold": 0.8,
            "unstranded_threshold": 0.1,
        }
    raise ValueError(f"unknown branch: {branch}")


def groovy_string(value: str) -> str:
    """Quote a literal without allowing interpolation in generated Nextflow configuration."""
    if any(c in value for c in "\n\r\0"):
        raise ValueError("configuration literal contains a control character")
    return "'" + value.replace("\\", "\\\\").replace("'", "\\'") + "'"


def prepare(
    intake: Intake, assets: ReferenceAssets, runtime: Runtime, output: Path, run_id: str
) -> None:
    """Write samplesheets and explicit branch budgets; do not launch incomplete assets."""
    runtime.validate_budgets()
    if not re.fullmatch(r"[a-zA-Z0-9_-]{1,80}", run_id):
        raise ValueError("run label is not safe")
    output.mkdir(parents=True, exist_ok=False)
    for branch, assay in [("dna", "germline_wgs"), ("rna", "bulk_rna")]:
        allocation = getattr(runtime, branch)
        directory = output / branch
        directory.mkdir()
        sheet = directory / "samples.csv"
        write_new(sheet, sample_sheet(intake, assay))
        write_new(
            directory / "params.json",
            canonical(parameters(branch, sheet, directory / "results", assets)) + b"\n",
        )
        settings = f"""executor.$local.cpus = {allocation.cpus}
executor.$local.memory = '{allocation.memory_gib} GB'
process.resourceLimits = [cpus: {allocation.cpus}, memory: '{allocation.memory_gib} GB', time: '{allocation.max_hours}h']
docker.runOptions = '--label bespoke.run={run_id} --network none'
"""
        if runtime.task_scratch is not None:
            settings += f"process.scratch = {groovy_string(str(runtime.task_scratch))}\n"
        write_new(directory / "runtime.config", settings.encode())
    write_new(output / "reference-assets.json", canonical(assets.model_dump(mode="json")) + b"\n")
    write_new(output / "runtime.json", canonical(runtime.model_dump(mode="json")) + b"\n")
    write_new(output / "run.json", canonical({"run_id": run_id}) + b"\n")


def environment(runtime: Runtime, *, offline: bool) -> dict[str, str]:
    return os.environ | {
        "NXF_HOME": str(runtime.nxf_home),
        "NXF_OFFLINE": str(offline).lower(),
        "NXF_DISABLE_CHECK_LATEST": "true",
        "NXF_ANSI_LOG": "false",
    }


def configuration_arguments(repository: Path, directory: Path, branch: str) -> list[str]:
    return [
        "-profile",
        "docker",
        "-c",
        str(repository / "workflows/stage1/base.config"),
        "-c",
        str(repository / f"workflows/stage1/{branch}.config"),
        "-c",
        str(directory / "runtime.config"),
    ]


def inspect(runtime: Runtime, repository: Path, root: Path, branch: str) -> dict:
    """Compile the selected workflow and record its reachable container declarations."""
    directory = root / branch
    pipeline = runtime.sarek if branch == "dna" else runtime.rnaseq
    command = [
        str(runtime.nextflow),
        "-log",
        str(directory / "inspect.nextflow.log"),
        "inspect",
        str(pipeline),
        "-format",
        "json",
        "-params-file",
        str(directory / "params.json"),
        *configuration_arguments(repository, directory, branch),
    ]
    with (
        (directory / "inspect.json").open("xb") as out,
        (directory / "inspect.stderr").open("xb") as err,
    ):
        subprocess.run(
            command,
            cwd=directory,
            env=environment(runtime, offline=False),
            stdout=out,
            stderr=err,
            check=True,
        )
    return read_json(directory / "inspect.json")


def launch_command(
    runtime: Runtime, repository: Path, root: Path, branch: str, attempt: Path, resume: bool
) -> list[str]:
    """Run local pinned sources, with explicit reports, work cache and container digests."""
    directory = root / branch
    pipeline = runtime.sarek if branch == "dna" else runtime.rnaseq
    command = [
        str(runtime.nextflow),
        "-log",
        str(attempt / "nextflow.log"),
        "run",
        str(pipeline),
        "-params-file",
        str(directory / "params.json"),
        *configuration_arguments(repository, directory, branch),
        "-c",
        str(directory / "containers.config"),
        "-work-dir",
        str(runtime.scratch / branch),
        "-with-trace",
        str(attempt / "trace.tsv"),
        "-with-report",
        str(attempt / "report.html"),
        "-with-timeline",
        str(attempt / "timeline.html"),
        "-ansi-log",
        "false",
    ]
    if resume:
        command.append("-resume")
    return command


def validate_prepared(root: Path, runtime: Runtime, run_id: str) -> None:
    """Reject unfinished assets and orphaned jobs before starting another coordinator.

    This checks availability, not biological qualification. Content identities and native
    validation evidence are checked separately before a handoff can become ready.
    """
    prepared = Runtime.model_validate(read_json(root / "runtime.json"))
    if read_json(root / "run.json") != {"run_id": run_id}:
        raise ValueError("run identity differs from the prepared container labels")
    if prepared != runtime:
        raise ValueError("runtime differs from the prepared workflow budgets")
    assets = ReferenceAssets.model_validate(read_json(root / "reference-assets.json"))
    required = [
        assets.fasta,
        Path(str(assets.fasta) + ".fai"),
        assets.dictionary,
        assets.gtf,
        assets.transcript_fasta,
        assets.vep_cache / "homo_sapiens/116_GRCh38/info.txt",
    ]
    required.extend(assets.star_index / name for name in ("Genome", "SA", "SAindex"))
    for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
        matches = list(assets.bwa_index.glob("*" + suffix))
        if len(matches) != 1:
            raise ValueError(f"expected one BWA-MEM2 index component: {suffix}")
        required.extend(matches)
    for branch in ("dna", "rna"):
        directory = root / branch
        required.extend(
            directory / name
            for name in ("params.json", "samples.csv", "runtime.config", "containers.config")
        )
        with (directory / "samples.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                required.extend(Path(row[key]) for key in ("fastq_1", "fastq_2"))
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"required workflow input is missing or empty: {path}")
    active = subprocess.run(
        ["docker", "ps", "-q", "--filter", f"label=bespoke.run={run_id}"],
        capture_output=True,
        text=True,
        check=True,
    )
    if active.stdout.strip():
        raise RuntimeError("previous run containers are still active; reconcile them before resume")


def launch(
    runtime: Runtime,
    repository: Path,
    root: Path,
    run_id: str,
    attempt_id: str,
    resume: bool = False,
) -> None:
    """Coordinate both branches. Failure or interruption cancels work and keeps its logs."""
    import signal
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from .io import run_lock
    from .runtime import Events, execute

    runtime.validate_budgets()
    if any(not re.fullmatch(r"[a-zA-Z0-9_-]{1,80}", value) for value in (run_id, attempt_id)):
        raise ValueError("unsafe run or attempt identifier")
    cancelled = threading.Event()
    root.mkdir(parents=True, exist_ok=True)
    with run_lock(root / "coordinator.lock"):
        validate_prepared(root, runtime, run_id)
        attempt = root / "attempts" / attempt_id
        attempt.mkdir(parents=True, exist_ok=False)
        previous = {}
        for signum in (signal.SIGINT, signal.SIGTERM):
            previous[signum] = signal.signal(signum, lambda *_: cancelled.set())

        def run(branch):
            directory = attempt / branch
            directory.mkdir()
            execute(
                launch_command(runtime, repository, root, branch, directory, resume),
                cwd=root / branch,
                log=directory / "launcher.log",
                events=Events(directory / "events.jsonl", run_id, attempt_id),
                stage=branch,
                disks=runtime.disks,
                timeout_seconds=getattr(runtime, branch).max_hours * 3600,
                poll_seconds=runtime.poll_seconds,
                environment=environment(runtime, offline=True),
                docker_cleanup=True,
                cancelled=cancelled,
            )

        errors = []
        try:
            with ThreadPoolExecutor(max_workers=2) as pool:
                futures = [pool.submit(run, branch) for branch in ("dna", "rna")]
                for future in as_completed(futures):
                    error = future.exception()
                    if error is not None:
                        cancelled.set()
                        errors.append(error)
            if errors:
                raise ExceptionGroup(
                    "native workflow execution failed; inspect attempt logs", errors
                )
            if cancelled.is_set():
                raise InterruptedError("workflow execution interrupted")
            write_new(
                attempt / "execution-complete.json",
                canonical(
                    {
                        "run_id": run_id,
                        "attempt_id": attempt_id,
                        "branches": ["dna", "rna"],
                        "handoff_ready": False,
                    }
                )
                + b"\n",
            )
        finally:
            for signum, handler in previous.items():
                signal.signal(signum, handler)
