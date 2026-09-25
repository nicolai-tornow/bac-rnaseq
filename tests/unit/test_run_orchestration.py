import pytest
from engine.python.config import load_config
from engine.python import run as R

PE = ("sample_id\tfastq_r1\tfastq_r2\tcondition\n"
      "s1\tA_R1.fq.gz\tA_R2.fq.gz\t7H9\ns2\tB_R1.fq.gz\tB_R2.fq.gz\tSCFM2\n")


def _setup(tmp_path, monkeypatch, qc=None, sheet=PE):
    """Replace every stage that needs real files; keep the orchestration logic."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))   # no real site.yaml
    ss = tmp_path / "s.tsv"
    ss.write_text(sheet)
    calls, deseq_called, processed = [], [], []

    def fake_runner(cmd, **kw):
        calls.append(cmd)

        class R0:
            returncode = 0
            stdout = ""
            stderr = "96.50% overall alignment rate\n"
        return R0()

    monkeypatch.setattr(R, "build_bundle", lambda *a, **k: {
        "saf": str(tmp_path / "labels.saf"), "index_prefix": str(tmp_path / "ref"),
        "n_features": 4970, "seqids": ["NC_010397.1"], "fasta": "x", "gene_ids": [],
        "extra_ids": [], "structural": {}})
    monkeypatch.setattr(R, "_process_sample",
                        lambda s, paired, out, bundle, threads, runner, **kw:
                        processed.append((s, paired, threads)) or {})
    monkeypatch.setattr(R, "_sample_integrity", lambda *a, **k: [])   # no real outputs here
    monkeypatch.setattr(R, "_strand_check", lambda *a, **k: {})
    monkeypatch.setattr(R, "_ncrna_fracs", lambda *a, **k: {})
    monkeypatch.setattr(R, "_reshape_counts", lambda *a, **k: None)
    monkeypatch.setattr(R, "_collect_qc", lambda *a, **k: qc or {"s1": {"verdict": "PASS", "reasons": []}})
    monkeypatch.setattr(R, "_run_deseq2", lambda *a, **k: deseq_called.append(True))
    return ss, fake_runner, calls, deseq_called, processed


def _cfg(**extra):
    return load_config({"run_name": "t", "reference": {"species": "mabs"},
                        "contrasts": {"explicit": [{"name": "S_vs_7", "numerator": "SCFM2",
                                                    "denominator": "7H9"}]}, **extra})


def test_orchestrator_uses_reverse_strand(tmp_path, monkeypatch):
    ss, runner, calls, _, _ = _setup(tmp_path, monkeypatch)
    rep = R.run_pipeline(_cfg(), tmp_path, tmp_path, str(ss), runner=runner, check_files=False)
    fc = [c for c in calls if c and c[0] == "featureCounts"][0]
    assert fc[fc.index("-s") + 1] == "2"
    assert rep["params"]["strandedness"] == "reverse"
    assert (tmp_path / "out/t/00_inputs/config.yaml").exists()
    assert (tmp_path / "out/t/00_inputs/samplesheet.tsv").exists()
    assert any(c[0] == "multiqc" for c in calls) and any(c[0] == "fastqc" for c in calls)


def test_qc_fail_halts_before_deseq2(tmp_path, monkeypatch):
    qc = {"s1": {"verdict": "FAIL", "reasons": ["low alignment"]},
          "s2": {"verdict": "PASS", "reasons": []}}
    ss, runner, _, deseq_called, _ = _setup(tmp_path, monkeypatch, qc)
    rep = R.run_pipeline(_cfg(), tmp_path, tmp_path, str(ss), runner=runner, check_files=False)
    assert rep["status"] == "qc_fail"
    assert not deseq_called          # DESeq2 must NOT run on QC-failed data


def test_allow_qc_fail_proceeds(tmp_path, monkeypatch):
    qc = {"s1": {"verdict": "FAIL", "reasons": ["x"]}, "s2": {"verdict": "PASS", "reasons": []}}
    ss, runner, _, deseq_called, _ = _setup(tmp_path, monkeypatch, qc)
    rep = R.run_pipeline(_cfg(), tmp_path, tmp_path, str(ss), runner=runner,
                         allow_qc_fail=True, check_files=False)
    assert rep["status"] == "ok"
    assert deseq_called


def test_threads_fall_back_to_site_yaml(tmp_path, monkeypatch):
    ss, runner, _, _, _ = _setup(tmp_path, monkeypatch)
    from engine.python import siteconfig
    siteconfig.write_site({"threads": 12})
    rep = R.run_pipeline(_cfg(), tmp_path, tmp_path, str(ss), runner=runner, check_files=False)
    assert rep["params"]["threads"] == 12
    rep = R.run_pipeline(_cfg(resources={"threads": 3}), tmp_path, tmp_path, str(ss),
                         runner=runner, check_files=False)
    assert rep["params"]["threads"] == 3


def test_mate1_only_runs_mixed_sheet_single_end(tmp_path, monkeypatch):
    mixed = ("sample_id\tfastq_r1\tfastq_r2\tcondition\n"
             "s1\tA.fq.gz\t\t7H9\ns2\tB_R1.fq.gz\tB_R2.fq.gz\tSCFM2\n")
    ss, runner, calls, _, processed = _setup(tmp_path, monkeypatch, sheet=mixed)
    with pytest.raises(ValueError, match="mate1_only"):
        R.run_pipeline(_cfg(), tmp_path, tmp_path, str(ss), runner=runner, check_files=False)
    R.run_pipeline(_cfg(reads={"layout": "mate1_only"}), tmp_path, tmp_path, str(ss),
                   runner=runner, check_files=False)
    assert all(paired is False and s.fastq_r2 is None for s, paired, _ in processed)
    fc = [c for c in calls if c and c[0] == "featureCounts"][0]
    assert "-p" not in fc


def test_contrast_typo_fails_before_any_work(tmp_path, monkeypatch):
    ss, runner, calls, _, _ = _setup(tmp_path, monkeypatch)
    cfg = load_config({"run_name": "t", "reference": {"species": "mabs"},
                       "contrasts": {"explicit": [{"name": "x", "numerator": "SCFM",
                                                   "denominator": "7H9"}]}})
    with pytest.raises(ValueError, match="SCFM"):
        R.run_pipeline(cfg, tmp_path, tmp_path, str(ss), runner=runner, check_files=False)
    assert calls == []


SIX = "sample_id\tfastq_r1\tcondition\n" + "".join(
    f"s{i}\tA{i}.fq.gz\t{'7H9' if i < 3 else 'SCFM2'}\n" for i in range(6))


def test_default_runs_one_sample_at_a_time(tmp_path, monkeypatch):
    ss, runner, _, _, processed = _setup(tmp_path, monkeypatch, sheet=SIX)
    rep = R.run_pipeline(_cfg(resources={"threads": 32}), tmp_path, tmp_path, str(ss),
                         runner=runner, check_files=False)
    assert (rep["params"]["parallel_samples"], rep["params"]["threads_per_sample"]) == (1, 32)
    assert {t for *_, t in processed} == {32} and len(processed) == 6


@pytest.mark.parametrize("fs,par,warned", [("nfs", 4, True), ("nfs", 1, False),
                                           ("ext2/ext3", 4, False)])
def test_parallel_samples_on_nfs_warns(tmp_path, monkeypatch, fs, par, warned):
    ss, runner, _, _, _ = _setup(tmp_path, monkeypatch, sheet=SIX)
    monkeypatch.setattr(R, "_fs_type", lambda path: fs)
    rep = R.run_pipeline(_cfg(resources={"threads": 32, "parallel_samples": par}), tmp_path,
                         tmp_path, str(ss), runner=runner, check_files=False)
    assert bool([w for w in rep["warnings"] if "NFS" in w]) is warned


def test_fs_type_reports_the_filesystem(tmp_path):
    assert R._fs_type(tmp_path)
    assert R._fs_type(tmp_path / "missing") is None


def test_parallel_samples_override_and_small_budget(tmp_path, monkeypatch):
    ss, runner, _, _, _ = _setup(tmp_path, monkeypatch, sheet=SIX)
    rep = R.run_pipeline(_cfg(resources={"threads": 32, "parallel_samples": 2}), tmp_path,
                         tmp_path, str(ss), runner=runner, check_files=False)
    assert (rep["params"]["parallel_samples"], rep["params"]["threads_per_sample"]) == (2, 16)
    rep = R.run_pipeline(_cfg(resources={"threads": 6}), tmp_path, tmp_path, str(ss),
                         runner=runner, check_files=False)
    assert (rep["params"]["parallel_samples"], rep["params"]["threads_per_sample"]) == (1, 6)


def test_count_reads_uses_a_checked_pipe(tmp_path):
    from engine.python.procs import CallableRunner
    from tests.unit.fake_tools import FakeTools, write_fastq_gz
    fq = write_fastq_gz(tmp_path / "a.fq.gz", 7)
    tools = FakeTools()
    assert R._count_reads(CallableRunner(tools), fq) == 7
    assert [c[0] for c in tools.calls] == ["gzip", "wc"]
    (tmp_path / "bad.fq.gz").write_bytes(b"not gzip")
    with pytest.raises(RuntimeError, match="bad.fq.gz"):
        R._count_reads(CallableRunner(FakeTools()), str(tmp_path / "bad.fq.gz"))


def test_multiqc_config_names_raw_fastqc_by_sample(tmp_path):
    import yaml
    from engine.python.samplesheet import Sample
    samples = [Sample("ctl1", "/d/lib7_S1_L001_R1_001.fastq.gz", "c", "/d/lib7_S1_L001_R2_001.fastq.gz"),
               Sample("trt1", "/d/shared.fq.gz", "t"), Sample("trt2", "/d/shared.fq.gz", "t")]
    cfg, names = R._multiqc_config(tmp_path, samples)
    c = yaml.safe_load(cfg.read_text())
    assert c["use_filename_as_sample_name"] == ["fastp"] and "strand_check" in c["fn_ignore_dirs"]
    assert c["table_sample_merge"] == {"R1": "_R1", "R2": "_R2"}
    rows = dict(l.split("\t") for l in names.read_text().splitlines())
    assert rows == {"lib7_S1_L001_R1_001": "ctl1_R1", "lib7_S1_L001_R2_001": "ctl1_R2"}  # shared file: ambiguous
