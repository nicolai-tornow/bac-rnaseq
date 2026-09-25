"""A failed, interrupted or refused run leaves no partial files and says where it stopped."""
import json
import os
import signal
import time
from pathlib import Path
import pytest
from engine.python.config import load_config
from engine.python import run as R
from engine.python.procs import CallableRunner, Result, Terminated
from engine.python.runlock import LOCK_NAME, RunLock, RunLocked
from tests.unit.fake_tools import FakeTools, write_fastq_gz


def _cfg(**extra):
    base = {"run_name": "t", "reference": {"species": "mabs"},
            "resources": {"threads": 4, "parallel_samples": 1},
            "contrasts": {"explicit": [{"name": "B_vs_A", "numerator": "B", "denominator": "A"}]}}
    return load_config({**base, **extra})


@pytest.fixture
def env(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    raw = tmp_path / "raw"
    raw.mkdir()
    rows = ["sample_id\tfastq_r1\tfastq_r2\tcondition"]
    for i in range(3):
        r1 = write_fastq_gz(raw / f"s{i}_R1.fq.gz", 10)
        r2 = write_fastq_gz(raw / f"s{i}_R2.fq.gz", 10)
        rows.append(f"s{i}\t{r1}\t{r2}\t{'A' if i % 2 == 0 else 'B'}")
    ss = tmp_path / "ss.tsv"
    ss.write_text("\n".join(rows) + "\n")
    monkeypatch.setattr(R, "build_bundle", lambda *a, **k: {
        "saf": str(tmp_path / "labels.saf"), "index_prefix": str(tmp_path / "ref"),
        "n_features": 1, "seqids": ["c"], "fasta": "x", "gene_ids": [], "extra_ids": [],
        "structural": {}, "stamp": {"fasta_md5": "abc"}})
    monkeypatch.setattr(R, "_ncrna_fracs", lambda *a, **k: {})
    monkeypatch.setattr(R, "_collect_qc", lambda *a, **k: {"s0": {"verdict": "PASS", "reasons": []}})
    return tmp_path / "work", ss


def _run(env, tools, cfg=None, runner=None):
    work, ss = env
    return R.run_pipeline(cfg or _cfg(), work, None, ss, runner=runner or tools)


def _out(env):
    return env[0] / "out" / "t"


def _report(env):
    return json.loads((_out(env) / "00_run_report.json").read_text())


def _leftovers(root):
    return [str(p) for p in Path(root).rglob("*")
            if p.name.endswith(".partial") or p.suffix == ".sam"
            or p.name.startswith(("temp-core-", "temp-sort-"))]


def _assert_clean_failure(env, stage, step=None, sample=None):
    rep = _report(env)
    assert rep["status"] == "failed", rep
    f = rep["failure"]
    assert f["stage"] == stage and f["error"]
    if step:
        assert f["step"] == step
    if sample:
        assert f["sample"] == sample
    assert _leftovers(env[0]) == []
    assert not (_out(env) / LOCK_NAME).exists()
    return f


def test_ok_run_promotes_everything_and_releases_the_lock(env):
    rep = _run(env, FakeTools())
    out = _out(env)
    assert rep["status"] == "ok" and rep["schema_version"] == "1.2" and rep["warnings"] == []
    assert (out / "06_deseq/results/B_vs_A.tsv").exists() and (out / "06_deseq/coldata.tsv").exists()
    assert all((out / f"04_align/s{i}.done.json").exists() for i in range(3))
    assert _leftovers(env[0]) == [] and not (out / LOCK_NAME).exists()


@pytest.mark.parametrize("tool,step", [("bowtie2", "bowtie2"), ("samtools sort", "sort")])
def test_alignment_failure_reports_and_cleans(env, tool, step):
    with pytest.raises(R.SampleError):
        _run(env, FakeTools(fail={tool: 1}))
    f = _assert_clean_failure(env, "trim_align", step, "s0")
    assert "04_align/.s0.partial" in f["cleaned"]


def test_featurecounts_failure_removes_its_temp_files(env):
    def crash(cmd):
        tmp = Path(cmd[cmd.index("--tmpDir") + 1])
        (tmp / "temp-core-000000-a1b2.sam").write_text("x")
        (tmp / "temp-sort-000000-a1b2-").write_text("x")
        return Result(1, "", "featureCounts crashed")

    with pytest.raises(RuntimeError, match="featureCounts"):
        _run(env, FakeTools(hooks={"featureCounts": crash}))
    f = _assert_clean_failure(env, "featurecounts")
    assert "05_counts/.featurecounts.partial" in f["cleaned"]


def test_sigterm_mid_sample(env):
    tools = FakeTools()
    runner = CallableRunner(tools)

    def term(cmd):
        os.kill(os.getpid(), signal.SIGTERM)
        runner.cancelled.wait(5)            # the main thread has reacted

    tools.hooks["bowtie2"] = term
    before = signal.getsignal(signal.SIGTERM)
    with pytest.raises(Terminated):
        _run(env, tools, runner=runner)
    _assert_clean_failure(env, "trim_align")
    assert signal.getsignal(signal.SIGTERM) == before


def test_sigterm_during_featurecounts(env):
    def term(cmd):
        tmp = Path(cmd[cmd.index("--tmpDir") + 1])
        (tmp / "temp-core-000000-a1b2.sam").write_text("x")
        os.kill(os.getpid(), signal.SIGTERM)
        time.sleep(2)

    with pytest.raises(Terminated):
        _run(env, FakeTools(hooks={"featureCounts": term}))
    _assert_clean_failure(env, "featurecounts")


def test_first_failure_stops_queued_samples(env):
    tools = FakeTools(hooks={"bowtie2": lambda cmd: Result(1, "", "boom")
                             if any("s0_R1" in a for a in cmd) else None})
    with pytest.raises(R.SampleError):
        _run(env, tools)
    fastp_inputs = [c[c.index("--in1") + 1] for c in tools.called("fastp")]
    assert len(fastp_inputs) == 1 and "s0_R1" in fastp_inputs[0]


def test_live_lock_refuses_before_any_work(env):
    out = _out(env)
    out.mkdir(parents=True)
    tools = FakeTools()
    with RunLock(out):
        with pytest.raises(RunLocked):
            _run(env, tools)
    assert tools.calls == [] and not (out / "00_run_report.json").exists()


def test_preflight_error_keeps_the_previous_report(env):
    _run(env, FakeTools())
    before = (_out(env) / "00_run_report.json").read_text()
    bad = _cfg(contrasts={"explicit": [{"name": "x", "numerator": "C", "denominator": "A"}]})
    with pytest.raises(ValueError):
        _run(env, FakeTools(), cfg=bad)
    assert (_out(env) / "00_run_report.json").read_text() == before


def test_deseq2_failure_keeps_previous_results(env):
    _run(env, FakeTools())
    res = _out(env) / "06_deseq/results/B_vs_A.tsv"
    res.write_text("previous")
    with pytest.raises(RuntimeError, match="Rscript"):
        _run(env, FakeTools(fail={"Rscript": 1}))
    _assert_clean_failure(env, "deseq2")
    assert res.read_text() == "previous"


def test_deseq2_success_replaces_results(env):
    _run(env, FakeTools())
    res = _out(env) / "06_deseq/results/B_vs_A.tsv"
    res.write_text("previous")
    _run(env, FakeTools())
    assert res.read_text().startswith("Gene")


def test_stale_staging_and_lock_are_cleared_at_start(env):
    out = _out(env)
    for d in ("04_align/.s9.partial", "06_deseq/.partial", "05_counts/.featurecounts.partial"):
        (out / d).mkdir(parents=True)
        (out / d / "junk").write_text("x")
    (out / LOCK_NAME).write_text(json.dumps({"host": os.uname().nodename, "pid": 99999999,
                                             "started": "x", "heartbeat": time.time()}))
    rep = _run(env, FakeTools())
    assert rep["status"] == "ok" and _leftovers(env[0]) == []
    assert rep["provenance"]["stale_lock_cleared"]["pid"] == 99999999
    assert "06_deseq/.partial" in rep["provenance"]["stale_staging_removed"]


def test_folder_from_the_old_engine(env):
    out = _out(env)
    (out / "04_align").mkdir(parents=True)
    (out / "05_counts").mkdir(parents=True)
    (out / "04_align/s0.sam").write_text("old")
    (out / "04_align/s0.bam").write_text("old")
    (out / "05_counts/temp-core-1.tmp").write_text("old")
    tools = FakeTools()
    rep = _run(env, tools)
    assert rep["status"] == "ok"
    assert not any(r["resumed"] for r in rep["provenance"]["reads"].values())
    assert not (out / "04_align/s0.sam").exists()


def test_work_dir_with_a_space(env, tmp_path):
    work = tmp_path / "my work"
    rep = R.run_pipeline(_cfg(), work, None, env[1], runner=FakeTools())
    assert rep["status"] == "ok" and _leftovers(work) == []
