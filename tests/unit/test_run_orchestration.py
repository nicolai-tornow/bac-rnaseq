from engine.python.config import load_config
from engine.python import run as R


def test_orchestrator_uses_reverse_strand(tmp_path, monkeypatch):
    ss = tmp_path / "s.tsv"
    ss.write_text("sample_id\tfastq_r1\tfastq_r2\tcondition\n"
                  "s1\tA_R1.fq.gz\tA_R2.fq.gz\t7H9\ns2\tB_R1.fq.gz\tB_R2.fq.gz\tSCFM2\n")
    calls = []

    def fake_runner(cmd, **kw):
        calls.append(cmd)

        class R0:
            returncode = 0
            stdout = ""
            stderr = "96.50% overall alignment rate\n"
        return R0()

    monkeypatch.setattr(R, "build_bundle", lambda *a, **k: {
        "saf": str(tmp_path / "labels.saf"), "index_prefix": str(tmp_path / "ref"),
        "n_features": 4970, "seqids": ["NC_010397.1"], "fasta": "x"})
    monkeypatch.setattr(R, "_run_deseq2", lambda *a, **k: None)
    monkeypatch.setattr(R, "_collect_qc", lambda *a, **k: {"s1": {"verdict": "PASS", "reasons": []}})
    monkeypatch.setattr(R, "_reshape_counts", lambda fc, out_tsv: out_tsv)
    (tmp_path / "labels.saf").write_text("GeneID\tChr\tStart\tEnd\tStrand\n")

    cfg = load_config({"run_name": "t", "reference": {"species": "mabs"},
                       "contrasts": {"explicit": [{"name": "S_vs_7", "numerator": "SCFM2", "denominator": "7H9"}]}})
    rep = R.run_pipeline(cfg, tmp_path, tmp_path, str(ss), runner=fake_runner)

    fc = [c for c in calls if c and c[0] == "featureCounts"][0]
    assert fc[fc.index("-s") + 1] == "2"
    assert rep["params"]["strandedness"] == "reverse"


def _fail_runner_setup(tmp_path, monkeypatch, qc):
    ss = tmp_path / "s.tsv"
    ss.write_text("sample_id\tfastq_r1\tfastq_r2\tcondition\n"
                  "s1\tA_R1.fq.gz\tA_R2.fq.gz\t7H9\ns2\tB_R1.fq.gz\tB_R2.fq.gz\tSCFM2\n")

    def fake_runner(cmd, **kw):
        class R0:
            returncode = 0
            stdout = ""
            stderr = "80.00% overall alignment rate\n"
        return R0()

    monkeypatch.setattr(R, "build_bundle", lambda *a, **k: {
        "saf": str(tmp_path / "labels.saf"), "index_prefix": str(tmp_path / "ref"),
        "n_features": 4970, "seqids": ["NC_010397.1"], "fasta": "x"})
    monkeypatch.setattr(R, "_reshape_counts", lambda fc, out_tsv: out_tsv)
    monkeypatch.setattr(R, "_collect_qc", lambda *a, **k: qc)
    deseq_called = []
    monkeypatch.setattr(R, "_run_deseq2", lambda *a, **k: deseq_called.append(True))
    (tmp_path / "labels.saf").write_text("GeneID\tChr\tStart\tEnd\tStrand\n")
    cfg = load_config({"run_name": "t", "reference": {"species": "mabs"},
                       "contrasts": {"explicit": [{"name": "a", "numerator": "SCFM2", "denominator": "7H9"}]}})
    return ss, fake_runner, deseq_called, cfg


def test_qc_fail_halts_before_deseq2(tmp_path, monkeypatch):
    qc = {"s1": {"verdict": "FAIL", "reasons": ["low alignment"]},
          "s2": {"verdict": "PASS", "reasons": []}}
    ss, runner, deseq_called, cfg = _fail_runner_setup(tmp_path, monkeypatch, qc)
    rep = R.run_pipeline(cfg, tmp_path, tmp_path, str(ss), runner=runner)
    assert rep["status"] == "qc_fail"
    assert not deseq_called          # DESeq2 must NOT run on QC-failed data


def test_allow_qc_fail_proceeds(tmp_path, monkeypatch):
    qc = {"s1": {"verdict": "FAIL", "reasons": ["x"]}, "s2": {"verdict": "PASS", "reasons": []}}
    ss, runner, deseq_called, cfg = _fail_runner_setup(tmp_path, monkeypatch, qc)
    rep = R.run_pipeline(cfg, tmp_path, tmp_path, str(ss), runner=runner, allow_qc_fail=True)
    assert rep["status"] == "ok"
    assert deseq_called
