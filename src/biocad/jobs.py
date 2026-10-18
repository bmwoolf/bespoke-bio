"""Durable, single-worker design queue with conservative persistent spend reservations."""

import fcntl
import json
import sqlite3
import threading
import time
from contextlib import contextmanager
from pathlib import Path


class JobConflict(ValueError):
    pass


class BudgetExceeded(ValueError):
    pass


class JobStore:
    def __init__(self, path, *, cap_usd=100, prior_spend_usd=3, service_reserve_usd=10):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.connect() as db:
            db.executescript("""
                PRAGMA journal_mode=DELETE;
                CREATE TABLE IF NOT EXISTS budget (
                    id INTEGER PRIMARY KEY CHECK(id=1), cap REAL NOT NULL,
                    prior REAL NOT NULL, service REAL NOT NULL);
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY, prepared TEXT NOT NULL, output TEXT NOT NULL,
                    status TEXT NOT NULL, reserved REAL NOT NULL, created REAL NOT NULL,
                    updated REAL NOT NULL, receipt TEXT);
                CREATE TABLE IF NOT EXISTS leases (
                    pod_id TEXT PRIMARY KEY, job_id TEXT NOT NULL,
                    state TEXT NOT NULL, updated REAL NOT NULL);
            """)
            db.execute(
                "INSERT OR IGNORE INTO budget VALUES(1,?,?,?)",
                (cap_usd, prior_spend_usd, service_reserve_usd),
            )

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.path, timeout=30)
        db.row_factory = sqlite3.Row
        try:
            with db:
                yield db
        finally:
            db.close()

    def budget(self):
        with self.connect() as db:
            row = db.execute("SELECT * FROM budget WHERE id=1").fetchone()
            reserved = db.execute("SELECT COALESCE(SUM(reserved),0) FROM jobs").fetchone()[0]
        return {
            "cap_usd": row["cap"],
            "prior_spend_usd": row["prior"],
            "service_storage_reserve_usd": row["service"],
            "job_reserved_usd": reserved,
            "remaining_unreserved_usd": row["cap"] - row["prior"] - row["service"] - reserved,
            "reservation_policy": "10 USD per job retained after execution; queued cancellations release their reservation",
        }

    def enqueue(self, job_id, prepared, output, *, reservation_usd=10):
        if reservation_usd < 10:
            raise ValueError("Managed jobs require at least a 10 USD reservation")
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            if db.execute("SELECT id FROM jobs WHERE id=?", (job_id,)).fetchone():
                raise JobConflict("Job already submitted")
            budget = db.execute("SELECT * FROM budget WHERE id=1").fetchone()
            reserved = db.execute("SELECT COALESCE(SUM(reserved),0) FROM jobs").fetchone()[0]
            if reserved + reservation_usd + budget["prior"] + budget["service"] > budget["cap"]:
                raise BudgetExceeded(
                    "Persistent 100 USD phase cap has insufficient unreserved budget"
                )
            now = time.time()
            db.execute(
                "INSERT INTO jobs VALUES(?,?,?,?,?,?,?,NULL)",
                (job_id, str(prepared), str(output), "queued", reservation_usd, now, now),
            )
        return self.get(job_id)

    def get(self, job_id):
        with self.connect() as db:
            row = db.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
        if row is None:
            raise KeyError(job_id)
        result = dict(row)
        result["receipt"] = json.loads(result["receipt"]) if result["receipt"] else None
        return result

    def claim(self):
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            row = db.execute(
                "SELECT * FROM jobs WHERE status='queued' ORDER BY created,id LIMIT 1"
            ).fetchone()
            if row is None:
                return None
            db.execute(
                "UPDATE jobs SET status='running',updated=? WHERE id=?", (time.time(), row["id"])
            )
        return self.get(row["id"])

    def finish(self, job_id, status, receipt):
        if status not in {"completed", "failed", "cancelled", "interrupted"}:
            raise ValueError("Invalid terminal status")
        with self.connect() as db:
            db.execute(
                "UPDATE jobs SET status=?,receipt=?,updated=? WHERE id=?",
                (status, json.dumps(receipt), time.time(), job_id),
            )
        return self.get(job_id)

    def cancel(self, job_id):
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            row = db.execute("SELECT status FROM jobs WHERE id=?", (job_id,)).fetchone()
            if row is None:
                raise KeyError(job_id)
            if row["status"] == "queued":
                db.execute(
                    "UPDATE jobs SET status='cancelled',reserved=0,updated=? WHERE id=?",
                    (time.time(), job_id),
                )
            elif row["status"] == "running":
                db.execute(
                    "UPDATE jobs SET status='cancel_requested',updated=? WHERE id=?",
                    (time.time(), job_id),
                )
            elif row["status"] != "cancel_requested":
                raise JobConflict("Job is already terminal")
        return self.get(job_id)

    def lease(self, job_id, pod_id, *, state="active"):
        if state not in {"active", "deleted"}:
            raise ValueError("Invalid pod lease state")
        with self.connect() as db:
            existing = db.execute("SELECT job_id FROM leases WHERE pod_id=?", (pod_id,)).fetchone()
            if existing and existing["job_id"] != job_id:
                raise JobConflict("Pod lease belongs to another job")
            db.execute(
                "INSERT INTO leases VALUES(?,?,?,?) ON CONFLICT(pod_id) DO UPDATE SET state=excluded.state,updated=excluded.updated",
                (pod_id, job_id, state, time.time()),
            )

    def active_leases(self):
        with self.connect() as db:
            return [dict(r) for r in db.execute("SELECT * FROM leases WHERE state='active'")]

    def interrupt_running(self):
        with self.connect() as db:
            rows = list(
                db.execute("SELECT id FROM jobs WHERE status IN ('running','cancel_requested')")
            )
        return [
            self.finish(
                r["id"],
                "interrupted",
                {
                    "status": "interrupted",
                    "error": "Service restarted; job is not automatically resumed",
                    "execution_performed": False,
                    "execution_attempted": True,
                },
            )
            for r in rows
        ]


class ManagedQueue:
    """One serial dispatcher; SQLite owns claims, reservations, cancellation and leases."""

    def __init__(
        self,
        store,
        runner,
        cleanup,
        on_status,
        *,
        poll_seconds=0.5,
        expires_at=None,
        startup_reconcile=None,
    ):
        self.store, self.runner, self.cleanup, self.on_status = store, runner, cleanup, on_status
        self.poll_seconds = poll_seconds
        self.stop_event = threading.Event()
        self.active_cancel = None
        self.active_id = None
        self.thread = None
        self.ready = False
        self.expires_at = expires_at
        self.dispatcher_lock = None
        self.startup_reconcile = startup_reconcile

    def start(self):
        self.dispatcher_lock = self.store.path.with_suffix(".dispatcher.lock").open("a+")
        try:
            try:
                fcntl.flock(self.dispatcher_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError("Only one managed API worker may dispatch this queue") from exc
            # Own the controller lock before inspecting or deleting ANY cloud lease.
            # A second API process must never reconcile the first process's active pod.
            if self.startup_reconcile:
                unresolved = self.startup_reconcile()
                if unresolved:
                    raise RuntimeError(
                        "Owned cloud leases require cleanup before accepting jobs: "
                        + ",".join(unresolved)
                    )
            for job in self.store.interrupt_running():
                self.on_status(job, "interrupted", job["receipt"])
            for lease in self.store.active_leases():
                if not self.cleanup(lease["pod_id"]):
                    raise RuntimeError(
                        "Unconfirmed owned-pod cleanup blocks managed execution: " + lease["pod_id"]
                    )
                self.store.lease(lease["job_id"], lease["pod_id"], state="deleted")
            self.ready = True
            self.thread = threading.Thread(
                target=self._loop, name="biocad-design-queue", daemon=True
            )
            self.thread.start()
        except BaseException:
            self.ready = False
            self.dispatcher_lock.close()
            raise

    def stop(self):
        self.ready = False
        self.stop_event.set()
        if self.active_cancel:
            self.active_cancel.set()
        if self.thread:
            self.thread.join(timeout=15)
        if self.dispatcher_lock and (self.thread is None or not self.thread.is_alive()):
            self.dispatcher_lock.close()

    def cancel(self, job_id):
        job = self.store.cancel(job_id)
        if self.active_id == job_id and self.active_cancel:
            self.active_cancel.set()
        self.on_status(job, job["status"], job.get("receipt") or {})
        return job

    def _loop(self):
        while not self.stop_event.is_set():
            job = self.store.claim()
            if job is None:
                self.stop_event.wait(self.poll_seconds)
                continue
            if self.expires_at and time.time() + 7200 > self.expires_at:
                receipt = {
                    "status": "failed",
                    "error": "Service lifetime cannot cover a complete two-hour job",
                    "execution_performed": False,
                    "execution_attempted": False,
                }
                job = self.store.finish(job["id"], "failed", receipt)
                self.on_status(job, "failed", receipt)
                continue
            self.active_id = job["id"]
            self.active_cancel = threading.Event()
            if self.store.get(job["id"])["status"] == "cancel_requested":
                self.active_cancel.set()
            try:
                self.on_status(job, "running", {})
                receipt = self.runner(
                    Path(job["prepared"]),
                    Path(job["output"]),
                    job_id=job["id"],
                    cancel_event=self.active_cancel,
                    lease_callback=lambda pod_id, state="active", job_id=job["id"]: (
                        self.store.lease(job_id, pod_id, state=state)
                    ),
                    timeout=7200,
                )
                state = (
                    "completed"
                    if receipt.get("status") == "completed"
                    else "cancelled"
                    if receipt.get("status") in {"cancelled", "canceled"}
                    else "failed"
                )
            except Exception as exc:  # noqa: BLE001 - persist every background failure
                receipt = {
                    "status": "failed",
                    "error": str(exc),
                    "execution_performed": False,
                    "execution_attempted": True,
                    "remote_cleanup_required": True,
                }
                state = "failed"
            job = self.store.finish(job["id"], state, receipt)
            try:
                self.on_status(job, state, receipt)
            except Exception:  # noqa: BLE001 - status persistence failure must block further paid work
                self.ready = False
                break
            self.active_id, self.active_cancel = None, None
            if (
                self.store.active_leases()
                or receipt.get("cleanup_required")
                or receipt.get("remote_cleanup_required")
            ):
                # Never launch another paid job while deletion of an owned pod is unconfirmed.
                self.ready = False
                break
