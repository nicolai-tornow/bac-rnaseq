"""Unaligned reads are captured in the main alignment and kept only for weak samples."""
import pytest
from engine.python import commands as C
from engine.python import run as R
from engine.python import unaligned as U
from engine.python.procs import CallableRunner
from engine.python.samplesheet import Sample
from tests.unit.fake_tools import FakeTools, write_fastq_gz
from tests.unit.test_run_failures import env, _run  # noqa: F401

B = {"index_prefix": "i", "stamp": {}}


def _sample(tmp_path, paired=True):
    raw = tmp_path / "raw"
    raw.mkdir(exist_ok=True)
    s = Sample("s1", write_fastq_gz(raw / "A_R1.fq.gz", 10), "c",
               write_fastq_gz(raw / "A_R2.fq.gz", 10) if paired else None)
    out = tmp_path / "out" / "t"
    for d in ("02_trimmed", "04_align"):
        (out / d).mkdir(parents=True, exist_ok=True)
    return s, out


@pytest.mark.parametrize("paired,flag", [(True, "--un-conc-gz"), (False, "--un-gz")])
def test_bowtie2_captures_unaligned_reads(tmp_path, paired, flag):
    s, out = _sample(tmp_path, paired)
    tools = FakeTools()
    R._process_sample(s, paired, out, B, 2, CallableRunner(tools))
    bt = tools.called("bowtie2")[0]
    assert ".s1.partial" in bt[bt.index(flag) + 1]
    marker = (out / "04_align/s1.done.json").read_text()
    assert "--un" not in marker                      # output-only flag: not in the fingerprint
    assert C.bowtie2_cmd("i", "r1", 1) == C.bowtie2_cmd("i", "r1", 1, un=None)


def test_unaligned_reads_kept_only_below_95_percent(tmp_path):
    s, out = _sample(tmp_path)
    R._process_sample(s, True, out, B, 2, CallableRunner(FakeTools(align_pct=80.0)))
    kept = sorted(p.name for p in (out / "qc/unaligned/s1").iterdir())
    assert kept == ["s1.unaligned_R1.fq.gz", "s1.unaligned_R2.fq.gz"]
    (out / "04_align/s1.done.json").unlink()          # re-process; now it aligns well
    R._process_sample(s, True, out, B, 2, CallableRunner(FakeTools(align_pct=97.0)))
    assert not (out / "qc/unaligned/s1").exists()


def test_gc_and_count(tmp_path):
    a = write_fastq_gz(tmp_path / "a.fq.gz", 2, seq="GGCCAATT")
    b = write_fastq_gz(tmp_path / "b.fq.gz", 1, seq="GGNN")
    assert U.gc_and_count([a]) == (2, 50.0)
    assert U.gc_and_count([a, b]) == (3, pytest.approx(100 * 10 / 18))


def test_weak_sample_gets_diagnostics_and_a_note_when_resumed(env, monkeypatch):
    monkeypatch.setattr(R, "_collect_qc", lambda *a, **k: {
        "s0": {"verdict": "WARN", "reasons": [], "alignment_pct": 80.0},
        "s1": {"verdict": "PASS", "reasons": [], "alignment_pct": 97.0}})
    rep = _run(env, FakeTools(align_pct=80.0))
    u = rep["samples_qc"]["s0"]["unaligned"]
    assert u["reads"] == 6 and u["gc_pct"] == 50.0
    assert u["rrna_like_frac"] is None and "rRNA" in u["note"]   # test bundle has no rRNA
    assert "unaligned" not in rep["samples_qc"]["s1"]
    import shutil
    shutil.rmtree(env[0] / "out/t/qc/unaligned/s0")
    rep = _run(env, FakeTools(align_pct=80.0))                    # s0 resumes: nothing captured
    assert rep["samples_qc"]["s0"]["unaligned"] is None
    assert "done.json" in rep["samples_qc"]["s0"]["unaligned_note"]
