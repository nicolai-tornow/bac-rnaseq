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
