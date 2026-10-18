"""Bounded subprocess execution with durable logs, cancellation and disk headroom checks."""

import os
import signal
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path
from threading import Event

from pydantic import Field

from .contracts import Record
from .io import canonical


class DiskLimit(Record):
    """Stop before a mounted filesystem crosses its free-byte or free-inode floor."""

    path: Path
    min_free_bytes: int = Field(ge=0)
    min_free_inodes: int = Field(ge=0)


class Allocation(Record):
    """One workflow's enforced container resources and independent executor budget."""

    cpus: int = Field(ge=1)
    memory_gib: int = Field(ge=1)
    max_hours: float = Field(gt=0)


class Runtime(Record):
    """Private runtime configuration; paths, machines and spending never enter biology."""

    root: Path
    # Durable Nextflow work outputs and cache dependencies, retained across retries.
    scratch: Path
    # Optional node-local task storage. Completed outputs return to the durable work area.
    task_scratch: Path | None = None
    nextflow: Path
    nxf_home: Path
    sarek: Path
    rnaseq: Path
    dna: Allocation
    rna: Allocation
    host_cpus: int = Field(ge=1)
    host_memory_gib: int = Field(ge=1)
    disks: list[DiskLimit] = Field(min_length=1)
    poll_seconds: float = Field(default=10, ge=0.1, le=60)

    def validate_budgets(self) -> None:
        if self.dna.cpus + self.rna.cpus > self.host_cpus:
            raise ValueError("combined DNA/RNA CPU budgets exceed the host allocation")
        if self.dna.memory_gib + self.rna.memory_gib > self.host_memory_gib:
            raise ValueError("combined DNA/RNA memory budgets exceed the host allocation")
        if self.task_scratch is not None and not any(
            self.task_scratch.is_relative_to(limit.path) for limit in self.disks
        ):
            raise ValueError("task scratch needs an independent disk headroom monitor")


class Events:
    """One writer per task log. Records contain identities and operations, never sequences."""

    def __init__(self, path: Path, run_id: str, attempt_id: str):
        self.path, self.run_id, self.attempt_id = path, run_id, attempt_id

    def write(self, stage: str, status: str, **fields) -> None:
        record = {
            "utc": datetime.now(UTC).isoformat(),
            "run_id": self.run_id,
            "attempt_id": self.attempt_id,
            "stage": stage,
            "status": status,
            **fields,
        }
        with self.path.open("ab") as stream:
            stream.write(canonical(record) + b"\n")
            stream.flush()
            os.fsync(stream.fileno())


def disk_status(limits: list[DiskLimit]) -> list[dict]:
    result = []
    for limit in limits:
        value = os.statvfs(limit.path)
        free = value.f_bavail * value.f_frsize
        if free < limit.min_free_bytes or value.f_favail < limit.min_free_inodes:
            raise RuntimeError(f"disk headroom floor crossed: {limit.path}")
        result.append({"path": str(limit.path), "free_bytes": free, "free_inodes": value.f_favail})
    return result


def stop_process(process: subprocess.Popen, grace_seconds: float = 20) -> None:
    """Terminate the process group, then reap it. Escalate only after the grace period."""
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        process.wait()
        return
    try:
        process.wait(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def stop_containers(run_id: str) -> None:
    """Docker containers outlive shell groups; stop only containers carrying this run's label."""
    result = subprocess.run(
        ["docker", "ps", "-q", "--filter", f"label=bespoke.run={run_id}"],
        capture_output=True,
        text=True,
        check=True,
    )
    containers = result.stdout.split()
    if containers:
        subprocess.run(["docker", "stop", "--time", "20", *containers], check=True)


def execute(
    command: list[str],
    *,
    cwd: Path,
    log: Path,
    events: Events,
    stage: str,
    disks: list[DiskLimit],
    timeout_seconds: float,
    poll_seconds: float = 10,
    environment: dict[str, str] | None = None,
    docker_cleanup: bool = False,
    cancelled: Event | None = None,
) -> None:
    """Run once. Any command, deadline or storage failure suppresses downstream publication."""
    disk_status(disks)
    started = time.monotonic()
    events.write(stage, "started", command=command, log=str(log))
    with log.open("xb") as stream:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            env=environment,
            stdout=stream,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        events.write(stage, "running", pid=process.pid)
        try:
            while process.poll() is None:
                if cancelled is not None and cancelled.is_set():
                    raise InterruptedError(f"{stage} cancelled")
                elapsed = time.monotonic() - started
                if elapsed >= timeout_seconds:
                    raise TimeoutError(f"{stage} exceeded its allocation deadline")
                status = disk_status(disks)
                if int(elapsed) % 60 < poll_seconds:
                    events.write(stage, "resources", seconds=elapsed, disks=status)
                time.sleep(poll_seconds)
            if process.returncode:
                raise subprocess.CalledProcessError(process.returncode, command)
            disk_status(disks)
        except BaseException as error:
            stop_process(process)
            if docker_cleanup:
                try:
                    stop_containers(events.run_id)
                except (OSError, subprocess.CalledProcessError) as cleanup_error:
                    error.add_note(f"container cleanup also failed: {cleanup_error}")
                    events.write(stage, "cleanup_failed", error=str(cleanup_error))
            events.write(
                stage,
                "failed",
                seconds=time.monotonic() - started,
                exit_code=process.returncode,
                error_category=type(error).__name__,
                error=str(error),
                log=str(log),
            )
            raise
    events.write(stage, "completed", seconds=time.monotonic() - started, exit_code=0, log=str(log))
