import subprocess
import sys
import yaml
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def test_validate_ok(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "r", "reference": {"species": "mabs"}}))
    r = subprocess.run([sys.executable, "-m", "engine.python.cli", "validate", str(cfg)],
                       cwd=REPO, capture_output=True)
    assert r.returncode == 0


def test_validate_bad(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "r", "reference": {"species": "custom"}}))
    r = subprocess.run([sys.executable, "-m", "engine.python.cli", "validate", str(cfg)],
                       cwd=REPO, capture_output=True)
    assert r.returncode == 2


def _bin(args, tmp_path, **env):
    import os
    e = {**os.environ, "XDG_CONFIG_HOME": str(tmp_path / "xdg"), **env}
    # Run from an unrelated directory: the launcher must not depend on the cwd.
    return subprocess.run([sys.executable, str(REPO / "bin" / "bac-rnaseq"), *args],
                          cwd=tmp_path, capture_output=True, text=True, env=e)


def test_launcher_validates_samplesheet_from_any_dir(tmp_path):
    r1 = tmp_path / "a.fq.gz"
    r1.write_text("")
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "r", "reference": {"species": "mabs"},
                                   "contrasts": {"explicit": [{"name": "B_vs_A", "numerator": "B",
                                                               "denominator": "A"}]}}))
    ss = tmp_path / "s.tsv"
    ss.write_text(f"sample_id\tfastq_r1\tcondition\na1\t{r1}\tA\nb1\t{r1}\tB\n")
    ok = _bin(["validate", str(cfg), "--samplesheet", str(ss)], tmp_path)
    assert ok.returncode == 0, ok.stderr
    assert "B_vs_A (B vs A)" in ok.stdout
    ss.write_text(f"sample_id\tfastq_r1\tcondition\na1\t{r1}\tA\nb1\t{r1}\tC\n")
    bad = _bin(["validate", str(cfg), "--samplesheet", str(ss)], tmp_path)
    assert bad.returncode == 2 and "condition 'B'" in bad.stderr


def test_doctor_does_not_overwrite_saved_threads(tmp_path):
    assert _bin(["doctor", "--save-threads", "7"], tmp_path).returncode == 0
    _bin(["doctor"], tmp_path)                      # report only
    site = (tmp_path / "xdg" / "bac-rnaseq" / "site.yaml").read_text()
    assert "threads: 7" in site


def test_cleanup_dry_run_refusal_and_yes(tmp_path):
    from tests.unit.test_cleanup import make_run, _files
    out, _ = make_run(tmp_path)
    before = _files(tmp_path)
    dry = _bin(["cleanup", str(out)], tmp_path)
    assert dry.returncode == 0, dry.stderr
    assert "trimmed" in dry.stdout and "nothing was deleted" in dry.stdout
    assert _files(tmp_path) == before
    yes = _bin(["cleanup", str(out), "--yes"], tmp_path)
    assert yes.returncode == 0, yes.stderr + yes.stdout
    assert "deleted 6 files" in yes.stdout
    assert not (out / "02_trimmed/s1_R1.fq.gz").exists() and (out / "04_align/s1.bam").exists()
    again = _bin(["cleanup", str(out)], tmp_path)                 # manifest just written
    assert again.returncode == 2 and "REFUSED" in again.stdout


def test_cleanup_refused_for_unfinished_run(tmp_path):
    from tests.unit.test_cleanup import make_run
    out, _ = make_run(tmp_path, status="qc_fail")
    r = _bin(["cleanup", str(out), "--yes"], tmp_path)
    assert r.returncode == 2 and "qc_fail" in r.stdout
    assert (out / "02_trimmed/s1_R1.fq.gz").exists()
