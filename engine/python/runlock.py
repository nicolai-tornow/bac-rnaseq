"""One process at a time per run folder: `out/<run>/.lock` holds host, PID and a
heartbeat refreshed every minute. On a shared NFS folder a PID cannot be checked
from another host, so a lock is live while its heartbeat is fresh and, on this
host, while its PID exists."""
from __future__ import annotations
import json
import os
import socket
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

LOCK_NAME = ".lock"
HEARTBEAT_S = 60
STALE_AFTER_S = 600


class RunLocked(RuntimeError):
    pass


def _host() -> str:
    return socket.gethostname()


def _pid_alive(pid) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except (OverflowError, ValueError, TypeError):
        return False
    return True


def read_lock(run_dir) -> dict | None:
    p = Path(run_dir) / LOCK_NAME
    if not p.exists():
        return None
    try:
        d = json.loads(p.read_text())
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}                       # being written, or damaged


def lock_status(run_dir, now: float | None = None) -> tuple[str, dict | None]:
    """("none" | "live" | "stale", lock content)."""
    p = Path(run_dir) / LOCK_NAME
    try:
        mtime = p.stat().st_mtime
    except FileNotFoundError:
        return "none", None
    info = read_lock(run_dir) or {}
    now = time.time() if now is None else now
    hb = info.get("heartbeat")
    if not isinstance(hb, (int, float)):
        hb = mtime
    if now - hb >= STALE_AFTER_S:
        return "stale", info
    if info.get("host") == _host() and not _pid_alive(info.get("pid")):
        return "stale", info
    return "live", info


def describe(info: dict | None, now: float | None = None) -> str:
    info = info or {}
    hb = info.get("heartbeat")
    age = f"{(now or time.time()) - hb:.0f} s ago" if isinstance(hb, (int, float)) else "unknown"
    return (f"host {info.get('host', '?')}, PID {info.get('pid', '?')}, "
            f"started {info.get('started', '?')}, last heartbeat {age}")


class RunLock:
    def __init__(self, run_dir, heartbeat_s: float = HEARTBEAT_S):
        self.dir = Path(run_dir)
        self.path = self.dir / LOCK_NAME
        self.heartbeat_s = heartbeat_s
        self.info: dict | None = None
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def acquire(self) -> dict | None:
        """Take the lock. Returns the content of a stale lock it cleared, else None."""
        cleared = None
        for _ in range(3):
            try:
                fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            except FileExistsError:
                status, info = lock_status(self.dir)
                if status == "live":
                    raise RunLocked(f"{self.dir} is in use by another run ({describe(info)}). "
                                    f"If that run is certainly gone, delete {self.path}.")
                if status == "stale":
                    cleared = info
                    self.path.unlink(missing_ok=True)
                continue
            self.info = {"host": _host(), "pid": os.getpid(),
                         "started": datetime.now(timezone.utc).isoformat(),
                         "heartbeat": time.time()}
            with os.fdopen(fd, "w") as fh:
                json.dump(self.info, fh)
            self._stop.clear()
            self._thread = threading.Thread(target=self._beat, daemon=True)
            self._thread.start()
            return cleared
        raise RunLocked(f"could not create {self.path}")

    def _ours(self) -> bool:
        cur = read_lock(self.dir) or {}
        return self.info is not None and all(cur.get(k) == self.info[k]
                                             for k in ("host", "pid", "started"))

    def _beat(self):
        while not self._stop.wait(self.heartbeat_s):
            if not self._ours():
                return
            self.info["heartbeat"] = time.time()
            tmp = self.dir / f"{LOCK_NAME}.{os.getpid()}.tmp"
            tmp.write_text(json.dumps(self.info))
            os.replace(tmp, self.path)

    def release(self):
        if self.info is None:
            return
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        if self._ours():
            self.path.unlink(missing_ok=True)
        self.info = None

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, *exc):
        self.release()
