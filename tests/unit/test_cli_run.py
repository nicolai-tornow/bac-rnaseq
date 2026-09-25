import yaml
from engine.python import cli
from engine.python import run as R


def test_cli_run_invokes_pipeline(tmp_path, monkeypatch):
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "t", "reference": {"species": "mabs"},
                                   "contrasts": {"explicit": [{"name": "a", "numerator": "SCFM2", "denominator": "7H9"}]}}))
    (tmp_path / "s.tsv").write_text("sample_id\tfastq_r1\tcondition\ns1\tA.fq.gz\t7H9\n")
    monkeypatch.setattr(R, "run_pipeline", lambda *a, **k: {"status": "ok", "run_name": "t"})
    rc = cli.main(["run", str(cfg), "--work-dir", str(tmp_path),
                   "--samplesheet", str(tmp_path / "s.tsv")])
    assert rc == 0


def _cli_env(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "t", "reference": {"species": "mabs"},
                                   "contrasts": {"explicit": [{"name": "a", "numerator": "SCFM2", "denominator": "7H9"}]}}))
    (tmp_path / "s.tsv").write_text("sample_id\tfastq_r1\tcondition\ns1\tA.fq.gz\t7H9\n")
    return ["run", str(cfg), "--work-dir", str(tmp_path), "--samplesheet", str(tmp_path / "s.tsv")]


import pytest  # noqa: E402
from engine.python.procs import Terminated  # noqa: E402
from engine.python.runlock import RunLocked  # noqa: E402


@pytest.mark.parametrize("exc,rc,text", [
    (RunLocked("in use by another run"), 2, "in use by another run"),
    (R.SampleError("s1", "bowtie2", "exited 1"), 1, "s1: bowtie2 failed"),
    (Terminated("signal 15"), 143, "stopped"),
    (KeyboardInterrupt(), 130, "stopped")])
def test_cli_run_failures_exit_cleanly(tmp_path, monkeypatch, capsys, exc, rc, text):
    def boom(*a, **k):
        raise exc
    monkeypatch.setattr(R, "run_pipeline", boom)
    assert cli.main(_cli_env(tmp_path)) == rc
    assert text in capsys.readouterr().err
