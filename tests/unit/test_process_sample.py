"""One sample: fastp -> bowtie2 | samtools sort, staged, promoted, marked complete."""
import json
from pathlib import Path
import pytest
from engine.python import run as R
from engine.python.procs import CallableRunner, Cancelled
from engine.python.samplesheet import Sample
from tests.unit.fake_tools import FakeTools, write_fastq_gz

BUNDLE = {"index_prefix": "/idx/ref", "stamp": {"fasta_md5": "abc"}}


def _setup(tmp_path, n=10):
    raw = tmp_path / "raw"
    raw.mkdir()
    r1 = write_fastq_gz(raw / "A_R1.fq.gz", n)
    r2 = write_fastq_gz(raw / "A_R2.fq.gz", n)
    out = tmp_path / "out" / "t"
    for d in ("02_trimmed", "04_align"):
        (out / d).mkdir(parents=True)
    return Sample("s1", r1, "c", r2), out


def _run(s, out, tools, paired=True, bundle=BUNDLE, cleaned=None, runner=None):
    return R._process_sample(s, paired, out, bundle, 4, runner or CallableRunner(tools),
                             cleaned=cleaned)


def _leftovers(root):
    return [p for p in Path(root).rglob("*") if p.name.endswith(".partial") or p.suffix == ".sam"]


def test_fresh_sample_writes_outputs_and_marker(tmp_path):
    s, out = _setup(tmp_path)
    tools = FakeTools()
    rep = _run(s, out, tools)
    assert rep["resumed"] is False and rep["reads_passed"] == 20
    for f in ("02_trimmed/s1_R1.fq.gz", "02_trimmed/s1_R2.fq.gz", "02_trimmed/s1.json",
              "02_trimmed/s1.html", "04_align/s1.bam", "04_align/s1.bam.bai",
              "04_align/s1.bowtie2.log", "04_align/s1.done.json"):
        assert (out / f).exists(), f
    m = json.loads((out / "04_align/s1.done.json").read_text())
    assert (m["bam_records"], m["reads_passed"], m["inputs"]["paired"]) == (100, 20, True)
    assert m["bam_md5"] and m["inputs"]["reference_fasta_md5"] == "abc"
    assert "overall alignment rate" in (out / "04_align/s1.bowtie2.log").read_text()
    sort = tools.called("samtools sort")[0]
    assert ".s1.partial" in sort[sort.index("-T") + 1] and sort[-1] == "-"
    assert "-S" not in tools.called("bowtie2")[0]           # no SAM file, piped
    assert _leftovers(tmp_path) == []


def test_rerun_resumes_without_trimming(tmp_path):
    s, out = _setup(tmp_path)
    _run(s, out, FakeTools())
    tools = FakeTools()
    assert _run(s, out, tools)["resumed"] is True
    assert tools.called("fastp") == [] and tools.called("bowtie2") == []


def test_resumes_when_trimmed_reads_were_cleaned_up(tmp_path):
    s, out = _setup(tmp_path)
    _run(s, out, FakeTools())
    for f in (out / "02_trimmed").glob("*.fq.gz"):
        f.unlink()
    tools = FakeTools()
    assert _run(s, out, tools)["resumed"] is True
    assert tools.called("fastp") == []


def test_missing_bam_is_retrimmed_and_realigned(tmp_path):
    s, out = _setup(tmp_path)
    _run(s, out, FakeTools())
    (out / "04_align/s1.bam").unlink()
    (out / "04_align/s1.bam.bai").unlink()
    tools = FakeTools()
    rep = _run(s, out, tools)
    assert rep["resumed"] is False and "BAM" in rep["resume_skipped"]
    assert tools.called("fastp") and tools.called("bowtie2")
    assert (out / "04_align/s1.bam").exists()


@pytest.mark.parametrize("change", ["bam", "raw_fastq", "layout", "reference", "no_marker"])
def test_changes_invalidate_the_marker(tmp_path, change):
    s, out = _setup(tmp_path)
    _run(s, out, FakeTools())
    paired, bundle = True, BUNDLE
    if change == "bam":
        (out / "04_align/s1.bam").write_text("FAKEBAM 99\n")
    elif change == "raw_fastq":
        write_fastq_gz(s.fastq_r1, 12)
        write_fastq_gz(s.fastq_r2, 12)
    elif change == "layout":
        s, paired = Sample("s1", s.fastq_r1, "c", None), False
    elif change == "reference":
        bundle = {**BUNDLE, "stamp": {"fasta_md5": "other"}}
    else:
        (out / "04_align/s1.done.json").unlink()        # made by the 0.2.0 engine
    tools = FakeTools()
    assert _run(s, out, tools, paired=paired, bundle=bundle)["resumed"] is False
    assert tools.called("fastp")


def test_missing_bai_is_rebuilt_on_resume(tmp_path):
    s, out = _setup(tmp_path)
    _run(s, out, FakeTools())
    (out / "04_align/s1.bam.bai").unlink()
    tools = FakeTools()
    assert _run(s, out, tools)["resumed"] is True
    assert tools.called("samtools index") and (out / "04_align/s1.bam.bai").exists()


@pytest.mark.parametrize("tool,step", [("bowtie2", "bowtie2"), ("samtools sort", "sort")])
def test_alignment_failure_leaves_nothing_behind(tmp_path, tool, step):
    s, out = _setup(tmp_path)
    cleaned = []
    with pytest.raises(R.SampleError) as e:
        _run(s, out, FakeTools(fail={tool: 1}), cleaned=cleaned)
    assert (e.value.sample_id, e.value.step) == ("s1", step)
    assert _leftovers(tmp_path) == []
    assert not list((out / "04_align").iterdir()) and not list((out / "02_trimmed").iterdir())
    assert sorted(cleaned) == ["02_trimmed/.s1.partial", "04_align/.s1.partial"]


def test_stale_sam_of_the_sample_is_removed_when_reprocessing(tmp_path):
    s, out = _setup(tmp_path)
    (out / "04_align/s1.sam").write_text("old engine leftover")
    _run(s, out, FakeTools())
    assert _leftovers(tmp_path) == []


def test_cancel_mid_sample_removes_staging(tmp_path):
    s, out = _setup(tmp_path)
    tools = FakeTools()
    runner = CallableRunner(tools)
    tools.hooks["bowtie2"] = lambda cmd: runner.cancel()
    with pytest.raises(Cancelled):
        _run(s, out, tools, runner=runner)
    assert _leftovers(tmp_path) == [] and not (out / "04_align/s1.done.json").exists()
