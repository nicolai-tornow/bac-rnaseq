import threading
import time
import pytest
from engine.python import procs


def test_run_captures_output_and_exit_code():
    r = procs.ProcRunner().run(["sh", "-c", "echo hi; echo err >&2; exit 3"])
    assert (r.returncode, r.stdout, r.stderr) == (3, "hi\n", "err\n")


def test_pipe_reports_both_exit_codes():
    p = procs.ProcRunner().pipe(["sh", "-c", "printf 'a\\nb\\n'; exit 4"], ["wc", "-l"])
    assert (p.rc1, p.rc2, p.stdout.strip()) == (4, 0, "2")


def test_pipe_writes_first_stderr_to_file(tmp_path):
    log = tmp_path / "bt2.log"
    p = procs.ProcRunner().pipe(["sh", "-c", "echo 96.5% overall >&2; echo x"], ["cat"],
                                stderr1=str(log))
    assert p.rc1 == 0 and p.stdout == "x\n"
    assert log.read_text() == "96.5% overall\n"


def test_cancel_kills_running_children():
    runner = procs.ProcRunner(grace=2)
    res = {}
    t = threading.Thread(target=lambda: res.update(r=runner.run(["sleep", "30"])))
    t.start()
    time.sleep(0.3)
    start = time.monotonic()
    runner.cancel()
    t.join(5)
    assert not t.is_alive() and time.monotonic() - start < 5
    assert res["r"].returncode != 0


def test_no_new_process_after_cancel():
    runner = procs.ProcRunner()
    runner.cancel()
    with pytest.raises(procs.Cancelled):
        runner.run(["true"])
    with pytest.raises(procs.Cancelled):
        runner.pipe(["true"], ["true"])


def test_callable_runner_pipes_output_into_second_command(tmp_path):
    seen = []

    def fake(cmd, input=None, **kw):
        seen.append((cmd[0], input))
        return procs.Result(5 if cmd[0] == "a" else 0, "out-" + cmd[0], "err-" + cmd[0])

    log = tmp_path / "log"
    p = procs.CallableRunner(fake).pipe(["a"], ["b"], stderr1=str(log))
    assert (p.rc1, p.rc2, p.stdout) == (5, 0, "out-b")
    assert seen == [("a", None), ("b", "out-a")]
    assert log.read_text() == "err-a"


def test_as_runner_wraps_callables_only():
    real = procs.ProcRunner()
    assert procs.as_runner(real) is real
    assert isinstance(procs.as_runner(lambda cmd, **kw: None), procs.CallableRunner)
