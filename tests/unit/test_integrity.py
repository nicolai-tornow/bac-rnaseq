"""After all samples and before counting: trimmed reads and BAMs match the marker."""
import gzip
import pytest
from engine.python import run as R
from engine.python.procs import CallableRunner, Result
from engine.python.samplesheet import Sample
from tests.unit.fake_tools import FakeTools, write_fastq_gz
from tests.unit.test_run_failures import env, _run, _report, _leftovers  # noqa: F401


def _sample(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    s = Sample("s1", write_fastq_gz(raw / "A_R1.fq.gz", 10), "c",
               write_fastq_gz(raw / "A_R2.fq.gz", 10))
    out = tmp_path / "out" / "t"
    for d in ("02_trimmed", "04_align"):
        (out / d).mkdir(parents=True)
    R._process_sample(s, True, out, {"index_prefix": "i", "stamp": {}}, 2,
                      CallableRunner(FakeTools()))
    return out


def _check(out, tools=None):
    return R._sample_integrity("s1", True, out, CallableRunner(tools or FakeTools()), 2)


def test_intact_sample_passes(tmp_path):
    assert _check(_sample(tmp_path)) == []


def test_corrupt_trimmed_fastq(tmp_path):
    out = _sample(tmp_path)
    data = (out / "02_trimmed/s1_R1.fq.gz").read_bytes()
    (out / "02_trimmed/s1_R1.fq.gz").write_bytes(data[: len(data) // 2])
    assert any("s1_R1.fq.gz" in p for p in _check(out))


def test_trimmed_read_count_differs_from_fastp(tmp_path):
    out = _sample(tmp_path)
    write_fastq_gz(out / "02_trimmed/s1_R2.fq.gz", 5)
    probs = _check(out)
    assert len(probs) == 1 and "s1_R2.fq.gz" in probs[0] and "5" in probs[0]


def test_cleaned_trimmed_reads_are_not_required(tmp_path):
    out = _sample(tmp_path)
    for f in (out / "02_trimmed").glob("*.fq.gz"):
        f.unlink()
    tools = FakeTools()
    assert _check(out, tools) == [] and tools.called("gzip") == []
    assert tools.called("samtools view")


def test_bam_record_count_differs_from_marker(tmp_path):
    out = _sample(tmp_path)
    tools = FakeTools(hooks={"samtools view": lambda cmd: Result(0, "7\n")})
    probs = _check(out, tools)
    assert len(probs) == 1 and "s1.bam" in probs[0]


def test_failing_sample_is_reprocessed_once(env, monkeypatch):
    real, seen = R._sample_integrity, []

    def flaky(sid, *a):
        seen.append(sid)
        return ["simulated"] if seen.count("s0") == 1 and sid == "s0" else real(sid, *a)

    monkeypatch.setattr(R, "_sample_integrity", flaky)
    tools = FakeTools()
    rep = _run(env, tools)
    assert rep["status"] == "ok"
    assert rep["provenance"]["reads"]["s0"]["reprocessed_after_integrity"] == ["simulated"]
    assert sum("s0_R1" in c[c.index("--in1") + 1] for c in tools.called("fastp")) == 2


def test_second_integrity_failure_fails_the_run(env, monkeypatch):
    monkeypatch.setattr(R, "_sample_integrity",
                        lambda sid, *a: ["always bad"] if sid == "s1" else [])
    with pytest.raises(R.SampleError, match="always bad"):
        _run(env, FakeTools())
    rep = _report(env)
    assert rep["status"] == "failed"
    assert (rep["failure"]["stage"], rep["failure"]["step"], rep["failure"]["sample"]) == \
        ("integrity", "integrity", "s1")
    assert _leftovers(env[0]) == []
