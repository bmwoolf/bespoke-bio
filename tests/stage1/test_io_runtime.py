import json
import os
import subprocess
import sys

import pytest

from bespoke_bio.stage1.contracts import Artifact
from bespoke_bio.stage1.io import artifact, publish_ready, read_json, run_lock, verify
from bespoke_bio.stage1.runtime import DiskLimit, Events, disk_status, execute


@pytest.mark.parametrize("text", ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e999}'])
def test_strict_json(text, tmp_path):
    path = tmp_path / "input.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_json(path)


def test_corruption_is_not_a_cache_hit(tmp_path):
    path = tmp_path / "data.txt"
    path.write_text("original")
    item = artifact(path, "fixture", "text")
    assert verify(item) == path
    path.write_text("changed!")
    with pytest.raises(ValueError, match="integrity"):
        verify(item)


def test_published_manifest_is_never_replaced(tmp_path):
    first, second = tmp_path / "a", tmp_path / "b"
    first.mkdir()
    second.mkdir()
    ready = tmp_path / "ready.json"
    publish_ready(first, ready, {"attempt": "first"})
    with pytest.raises(FileExistsError):
        publish_ready(second, ready, {"attempt": "second"})
    assert read_json(ready) == {"attempt": "first"}


def test_duplicate_coordinator_cannot_acquire_lock(tmp_path):
    path = tmp_path / "run.lock"
    with run_lock(path), pytest.raises(BlockingIOError), run_lock(path):
        pytest.fail("duplicate coordinator entered")


@pytest.mark.parametrize("url", ["https://host/file?sig=secret", "https://user:pass@host/file"])
def test_signed_locations_rejected(url):
    with pytest.raises(ValueError):
        Artifact(role="x", location=url, bytes=1, sha256="0" * 64, format="bin")


def test_disk_floor_prevents_start(tmp_path):
    limit = DiskLimit(path=tmp_path, min_free_bytes=2**63, min_free_inodes=0)
    with pytest.raises(RuntimeError, match="headroom"):
        disk_status([limit])


def test_native_failure_retains_error_and_never_completes(tmp_path):
    events = Events(tmp_path / "events.jsonl", "run", "attempt")
    with pytest.raises(subprocess.CalledProcessError):
        execute(
            [sys.executable, "-c", "import sys; print('native diagnostic'); sys.exit(7)"],
            cwd=tmp_path,
            log=tmp_path / "native.log",
            events=events,
            stage="fixture",
            disks=[],
            timeout_seconds=10,
            poll_seconds=0.01,
        )
    rows = [json.loads(s) for s in events.path.read_text().splitlines()]
    assert rows[-1]["status"] == "failed"
    assert rows[-1]["exit_code"] == 7
    assert "native diagnostic" in (tmp_path / "native.log").read_text()
    assert not any(row["status"] == "completed" for row in rows)


def test_deadline_kills_and_reaps_process(tmp_path):
    events = Events(tmp_path / "events.jsonl", "run", "attempt")
    with pytest.raises(TimeoutError):
        execute(
            [sys.executable, "-c", "import time; time.sleep(60)"],
            cwd=tmp_path,
            log=tmp_path / "native.log",
            events=events,
            stage="fixture",
            disks=[],
            timeout_seconds=0.05,
            poll_seconds=0.01,
        )
    rows = [json.loads(s) for s in events.path.read_text().splitlines()]
    pid = next(row["pid"] for row in rows if row["status"] == "running")
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


def test_cancelled_branch_is_reaped_and_cannot_publish_success(tmp_path):
    from threading import Event, Timer

    cancelled = Event()
    events = Events(tmp_path / "events.jsonl", "run", "attempt")
    timer = Timer(0.03, cancelled.set)
    timer.start()
    try:
        with pytest.raises(InterruptedError):
            execute(
                [sys.executable, "-c", "import time; time.sleep(60)"],
                cwd=tmp_path,
                log=tmp_path / "native.log",
                events=events,
                stage="fixture",
                disks=[],
                timeout_seconds=10,
                poll_seconds=0.01,
                cancelled=cancelled,
            )
    finally:
        timer.join()
    rows = [json.loads(line) for line in events.path.read_text().splitlines()]
    assert rows[-1]["status"] == "failed"
    assert rows[-1]["error_category"] == "InterruptedError"
    assert not any(row["status"] == "completed" for row in rows)


def test_ephemeral_scratch_requires_its_own_headroom_monitor(tmp_path):
    from bespoke_bio.stage1.runtime import Runtime

    settings = {
        "root": tmp_path / "durable",
        "scratch": tmp_path / "durable/work",
        "task_scratch": tmp_path / "ephemeral/tasks",
        "nextflow": tmp_path / "nextflow",
        "nxf_home": tmp_path / "nxf",
        "sarek": tmp_path / "sarek",
        "rnaseq": tmp_path / "rnaseq",
        "dna": {"cpus": 2, "memory_gib": 4, "max_hours": 1},
        "rna": {"cpus": 2, "memory_gib": 4, "max_hours": 1},
        "host_cpus": 4,
        "host_memory_gib": 8,
        "disks": [{"path": tmp_path / "durable", "min_free_bytes": 1, "min_free_inodes": 1}],
    }
    with pytest.raises(ValueError, match="scratch"):
        Runtime(**settings).validate_budgets()
    settings["disks"].append(
        {"path": tmp_path / "ephemeral", "min_free_bytes": 1, "min_free_inodes": 1}
    )
    Runtime(**settings).validate_budgets()
