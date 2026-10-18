import threading
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from biocad.jobs import BudgetExceeded, JobConflict, JobStore, ManagedQueue


def test_atomic_budget_reservations_persist_and_queued_cancellation_releases(tmp_path):
    store = JobStore(tmp_path / "jobs.sqlite3")

    def enqueue(index):
        try:
            store.enqueue(str(index), tmp_path / "prepared", tmp_path / str(index))
            return True
        except BudgetExceeded:
            return False

    with ThreadPoolExecutor(max_workers=12) as pool:
        accepted = list(pool.map(enqueue, range(12)))
    assert sum(accepted) == 8
    assert store.budget()["remaining_unreserved_usd"] == 7
    reopened = JobStore(store.path)
    assert reopened.budget()["job_reserved_usd"] == 80
    queued = next(str(i) for i, accepted in enumerate(accepted) if accepted)
    reopened.cancel(queued)
    assert reopened.budget()["job_reserved_usd"] == 70
    reopened.enqueue("replacement", tmp_path / "prepared", tmp_path / "output")
    with pytest.raises(JobConflict):
        reopened.enqueue("replacement", tmp_path, tmp_path)
    with store.connect() as db:
        assert db.execute("PRAGMA journal_mode").fetchone()[0] == "delete"


def test_restart_interrupts_inflight_and_cleans_owned_lease_without_resume(tmp_path):
    store = JobStore(tmp_path / "jobs.sqlite3")
    store.enqueue("old", tmp_path, tmp_path / "old")
    store.claim()
    store.lease("old", "owned-pod")
    cleaned, states, called = [], [], []
    queue = ManagedQueue(
        store,
        lambda *a, **k: called.append(True),
        lambda pod: cleaned.append(pod) or True,
        lambda job, state, receipt: states.append(state),
        poll_seconds=0.01,
    )
    queue.start()
    queue.stop()
    assert cleaned == ["owned-pod"] and not called
    assert store.get("old")["status"] == "interrupted"
    assert states == ["interrupted"] and store.active_leases() == []
    assert store.budget()["job_reserved_usd"] == 10


def test_unconfirmed_cleanup_blocks_queue_and_second_dispatcher(tmp_path):
    store = JobStore(tmp_path / "jobs.sqlite3")
    store.enqueue("old", tmp_path, tmp_path / "out")
    store.lease("old", "owned-pod")
    blocked = ManagedQueue(store, lambda *a, **k: None, lambda pod: False, lambda *a: None)
    with pytest.raises(RuntimeError, match="cleanup blocks"):
        blocked.start()
    assert not blocked.ready
    blocked.stop()
    store.lease("old", "owned-pod", state="deleted")
    store.cancel("old")
    first = ManagedQueue(store, lambda *a, **k: None, lambda pod: True, lambda *a: None)
    second = ManagedQueue(store, lambda *a, **k: None, lambda pod: True, lambda *a: None)
    first.start()
    try:
        with pytest.raises(RuntimeError, match="one managed API worker"):
            second.start()
    finally:
        first.stop()


def test_active_cancellation_keeps_reservation_and_dispatch_deadline_rejects(tmp_path):
    store = JobStore(tmp_path / "jobs.sqlite3")
    entered, done = threading.Event(), threading.Event()

    def runner(prepared, output, *, job_id, cancel_event, lease_callback, timeout):
        lease_callback("owned-active")
        entered.set()
        assert cancel_event.wait(2)
        lease_callback("owned-active", state="deleted")
        return {"status": "cancelled", "execution_performed": True}

    def status(job, state, receipt):
        if state == "cancelled":
            done.set()

    store.enqueue("active", tmp_path, tmp_path / "output")
    queue = ManagedQueue(store, runner, lambda pod: True, status, poll_seconds=0.01)
    queue.start()
    try:
        assert entered.wait(2)
        queue.cancel("active")
        assert done.wait(2)
    finally:
        queue.stop()
    assert store.get("active")["status"] == "cancelled"
    assert store.budget()["job_reserved_usd"] == 10
    store.enqueue("too-late", tmp_path, tmp_path / "out2")
    rejected = threading.Event()
    late = ManagedQueue(
        store,
        lambda *a, **k: pytest.fail("must not run"),
        lambda pod: True,
        lambda job, state, receipt: rejected.set(),
        expires_at=time.time() + 100,
    )
    late.start()
    try:
        assert rejected.wait(2)
    finally:
        late.stop()
    assert store.get("too-late")["status"] == "failed"
    assert store.get("too-late")["receipt"]["execution_performed"] is False


def test_controller_lock_precedes_cloud_reconciliation(tmp_path):
    store = JobStore(tmp_path / "jobs.sqlite3")
    first_reconciled, second_reconciled = [], []
    first = ManagedQueue(
        store,
        lambda *a, **k: None,
        lambda pod: True,
        lambda *a: None,
        startup_reconcile=lambda: first_reconciled.append(True) or [],
    )
    second = ManagedQueue(
        store,
        lambda *a, **k: None,
        lambda pod: pytest.fail("must not clean an active controller pod"),
        lambda *a: None,
        startup_reconcile=lambda: second_reconciled.append(True) or [],
    )
    first.start()
    try:
        with pytest.raises(RuntimeError, match="one managed API worker"):
            second.start()
        assert first_reconciled == [True]
        assert second_reconciled == []
        assert second.dispatcher_lock.closed
    finally:
        first.stop()


def test_unresolved_cloud_preflight_blocks_startup_and_releases_controller_lock(tmp_path):
    store = JobStore(tmp_path / "jobs.sqlite3")
    blocked = ManagedQueue(
        store,
        lambda *a, **k: pytest.fail("must not dispatch"),
        lambda pod: True,
        lambda *a: None,
        startup_reconcile=lambda: ["owned-pod"],
    )
    with pytest.raises(RuntimeError, match="Owned cloud leases require cleanup"):
        blocked.start()
    assert not blocked.ready and blocked.dispatcher_lock.closed
    recovered = ManagedQueue(
        store, lambda *a, **k: None, lambda pod: True, lambda *a: None, startup_reconcile=list
    )
    recovered.start()
    recovered.stop()
