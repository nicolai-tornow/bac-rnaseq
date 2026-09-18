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
                        lambda s, paired, *a, **k: processed.append((s, paired)) or {})
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
    assert all(paired is False and s.fastq_r2 is None for s, paired in processed)
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
