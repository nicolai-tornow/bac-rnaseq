"""Tool processes the pipeline can stop: a failed sample or a SIGTERM kills every
child still running, so no tool keeps writing into files that are being removed."""
from __future__ import annotations
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path


class Cancelled(Exception):
    """The run is being cancelled; no new tool process is started."""


class Terminated(BaseException):
    """SIGTERM arrived (a BaseException, like KeyboardInterrupt)."""


@dataclass
class Result:
    returncode: int
    stdout: str = ""
    stderr: str = ""


@dataclass
class PipeResult:
    rc1: int
    rc2: int
    stdout: str = ""    # of the second command
    stderr: str = ""    # of the second command
    stderr1: str = ""   # of the first command, unless it went to a file


def _text(b) -> str:
    return b.decode(errors="replace") if isinstance(b, bytes) else (b or "")


class ProcRunner:
    """subprocess.Popen with a registry of live children."""

    def __init__(self, grace: float = 10.0):
        self.grace = grace
        self.cancelled = threading.Event()
        self._live: set = set()
        self._lock = threading.Lock()

    def _start(self, cmd, **kw):
        with self._lock:
            if self.cancelled.is_set():
                raise Cancelled(cmd[0])
            p = subprocess.Popen(cmd, **kw)
            self._live.add(p)
        return p

    def _done(self, *ps):
        with self._lock:
            for p in ps:
                self._live.discard(p)

    def _stop(self, *ps):
        """SIGTERM, then SIGKILL after `grace`, and reap: used when the waiting thread
        is interrupted (Ctrl-C, SIGTERM), so no child outlives its call."""
        for p in ps:
            if p.poll() is None:
                try:
                    p.terminate()
                except OSError:
                    pass
        deadline = time.monotonic() + self.grace
        for p in ps:
            try:
                p.wait(max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait()

    def run(self, cmd) -> Result:
        p = self._start(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            out, err = p.communicate()
        except BaseException:
            self._stop(p)
            raise
        finally:
            self._done(p)
        return Result(p.returncode, _text(out), _text(err))

    def pipe(self, cmd1, cmd2, stderr1=None) -> PipeResult:
        """cmd1 | cmd2, both exit codes kept. stderr1: file for cmd1's stderr."""
        err1 = open(stderr1, "wb") if stderr1 else tempfile.TemporaryFile()
        p1 = p2 = None
        try:
            p1 = self._start(cmd1, stdout=subprocess.PIPE, stderr=err1)
            try:
                p2 = self._start(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            except BaseException:
                p1.kill()
                p1.wait()
                raise
            p1.stdout.close()          # cmd1 gets SIGPIPE if cmd2 exits early
            try:
                out, err = p2.communicate()
                rc1 = p1.wait()
            except BaseException:
                self._stop(p1, p2)
                raise
            s1 = ""
            if not stderr1:
                err1.seek(0)
                s1 = _text(err1.read())
            return PipeResult(rc1, p2.returncode, _text(out), _text(err), s1)
        finally:
            err1.close()
            self._done(*(p for p in (p1, p2) if p is not None))

    def cancel(self):
        """Stop starting processes; SIGTERM every live child, SIGKILL after `grace`."""
        self.cancelled.set()
        with self._lock:
            live = list(self._live)
        for p in live:
            if p.poll() is None:
                try:
                    p.terminate()
                except OSError:
                    pass
        deadline = time.monotonic() + self.grace
        while any(p.poll() is None for p in live) and time.monotonic() < deadline:
            time.sleep(0.05)
        for p in live:
            if p.poll() is None:
                p.kill()


def _as_result(r) -> Result:
    return Result(getattr(r, "returncode", 0), getattr(r, "stdout", "") or "",
                  getattr(r, "stderr", "") or "")


class CallableRunner:
    """Adapter for a subprocess.run-like callable (the test fakes)."""

    def __init__(self, fn):
        self.fn = fn
        self.cancelled = threading.Event()

    def _guard(self, cmd):
        if self.cancelled.is_set():
            raise Cancelled(cmd[0])

    def run(self, cmd) -> Result:
        self._guard(cmd)
        return _as_result(self.fn(cmd, capture_output=True, text=True))

    def pipe(self, cmd1, cmd2, stderr1=None) -> PipeResult:
        r1 = self.run(cmd1)
        self._guard(cmd2)
        r2 = _as_result(self.fn(cmd2, capture_output=True, text=True, input=r1.stdout))
        if stderr1:
            Path(stderr1).write_text(r1.stderr)
        return PipeResult(r1.returncode, r2.returncode, r2.stdout, r2.stderr,
                          "" if stderr1 else r1.stderr)

    def cancel(self):
        self.cancelled.set()


def as_runner(r):
    return r if hasattr(r, "pipe") else CallableRunner(r)
