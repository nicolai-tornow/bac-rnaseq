import json
import os
import socket
import time
import pytest
from engine.python import runlock as L


def _write_lock(d, **kw):
    info = {"host": socket.gethostname(), "pid": os.getpid(), "started": "x",
            "heartbeat": time.time(), **kw}
    (d / L.LOCK_NAME).write_text(json.dumps(info))
    return info


def test_acquire_writes_lock_and_release_removes_it(tmp_path):
    lk = L.RunLock(tmp_path)
    assert lk.acquire() is None
    info = json.loads((tmp_path / L.LOCK_NAME).read_text())
    assert info["host"] == socket.gethostname() and info["pid"] == os.getpid()
    assert {"started", "heartbeat"} <= set(info)
    assert L.lock_status(tmp_path)[0] == "live"
    lk.release()
    assert not (tmp_path / L.LOCK_NAME).exists()
    assert L.lock_status(tmp_path) == ("none", None)


def test_second_acquire_is_refused_while_live(tmp_path):
    lk = L.RunLock(tmp_path)
    lk.acquire()
    try:
        with pytest.raises(L.RunLocked, match="in use"):
            L.RunLock(tmp_path).acquire()
    finally:
        lk.release()


def test_fresh_lock_from_another_host_is_live(tmp_path):
    _write_lock(tmp_path, host="otherhost", pid=1)
    assert L.lock_status(tmp_path)[0] == "live"
    with pytest.raises(L.RunLocked, match="otherhost"):
        L.RunLock(tmp_path).acquire()


def test_old_heartbeat_from_another_host_is_stale_and_cleared(tmp_path):
    old = _write_lock(tmp_path, host="otherhost", pid=1, heartbeat=time.time() - 601)
    assert L.lock_status(tmp_path)[0] == "stale"
    lk = L.RunLock(tmp_path)
    assert lk.acquire() == old
    lk.release()


def test_dead_pid_on_this_host_is_stale(tmp_path):
    _write_lock(tmp_path, pid=99999999)
    assert L.lock_status(tmp_path)[0] == "stale"
    lk = L.RunLock(tmp_path)
    assert lk.acquire()["pid"] == 99999999
    lk.release()


def test_heartbeat_advances(tmp_path):
    lk = L.RunLock(tmp_path, heartbeat_s=0.05)
    lk.acquire()
    try:
        hb1 = json.loads((tmp_path / L.LOCK_NAME).read_text())["heartbeat"]
        time.sleep(0.3)
        hb2 = json.loads((tmp_path / L.LOCK_NAME).read_text())["heartbeat"]
        assert hb2 > hb1
    finally:
        lk.release()


def test_release_leaves_someone_elses_lock(tmp_path):
    lk = L.RunLock(tmp_path)
    lk.acquire()
    _write_lock(tmp_path, host="otherhost", pid=1)      # taken over meanwhile
    lk.release()
    assert json.loads((tmp_path / L.LOCK_NAME).read_text())["host"] == "otherhost"


def test_unreadable_lock_judged_by_age(tmp_path):
    p = tmp_path / L.LOCK_NAME
    p.write_text("")
    assert L.lock_status(tmp_path)[0] == "live"
    os.utime(p, (time.time() - 601, time.time() - 601))
    assert L.lock_status(tmp_path)[0] == "stale"


def test_heartbeat_survives_a_write_error(tmp_path, monkeypatch):
    lk = L.RunLock(tmp_path, heartbeat_s=0.05)
    lk.acquire()
    real, fails = os.replace, []

    def flaky(a, b):
        if not fails:
            fails.append(1)
            raise OSError(5, "Input/output error")
        return real(a, b)

    monkeypatch.setattr(L.os, "replace", flaky)
    try:
        hb1 = json.loads((tmp_path / L.LOCK_NAME).read_text())["heartbeat"]
        time.sleep(0.4)
        hb2 = json.loads((tmp_path / L.LOCK_NAME).read_text())["heartbeat"]
        assert fails and hb2 > hb1
    finally:
        lk.release()


def test_refusal_does_not_invite_deleting_the_lock(tmp_path):
    _write_lock(tmp_path, host="otherhost", pid=1)
    with pytest.raises(L.RunLocked) as e:
        L.RunLock(tmp_path).acquire()
    assert "delete" not in str(e.value).lower() and "clears itself" in str(e.value)
