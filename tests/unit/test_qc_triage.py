import pandas as pd
from engine.python.qc_triage import (parse_bowtie2_log, parse_featurecounts_summary,
                                     infer_strandedness, ncrna_fractions, triage_sample)

BT2 = """1000 reads; of these:
  1000 (100.00%) were paired; of these:
    50 (5.00%) aligned concordantly 0 times
    900 (90.00%) aligned concordantly exactly 1 time
    50 (5.00%) aligned concordantly >1 times
96.50% overall alignment rate
"""

# Jackson 1239 (boulder, 2026-09-18): correct reverse library, 79% of reads on
# structural ncRNAs the GFF did not annotate.
JACKSON = {"reverse": 0.169, "forward": 0.016, "unstranded": 0.179}


def test_parse_bowtie2_rate():
    assert parse_bowtie2_log(BT2) == 96.5


def test_featurecounts_summary_is_per_sample(tmp_path):
    p = tmp_path / "fc.txt.summary"
    p.write_text("Status\t/w/04_align/good.bam\t/w/04_align/bad.bam\n"
                 "Assigned\t800\t100\n"
                 "Unassigned_NoFeatures\t200\t900\n")
    s = parse_featurecounts_summary(p)
    assert s["good"]["assigned_frac"] == 0.8
    assert s["bad"]["assigned_frac"] == 0.1        # not the pooled 0.45
    assert s["bad"]["Unassigned_NoFeatures"] == 900


def test_strand_confirmed_on_jackson_numbers():
    s = infer_strandedness(JACKSON, declared="reverse")
    assert s["inferred"] == "reverse" and s["verdict"] == "PASS"


def test_strand_wrong_is_fail():
    s = infer_strandedness(JACKSON, declared="forward")
    assert s["verdict"] == "FAIL" and "reverse" in s["message"]


def test_unstranded_library_detected():
    s = infer_strandedness({"reverse": 0.40, "forward": 0.38, "unstranded": 0.78}, "reverse")
    assert s["inferred"] == "unstranded" and s["verdict"] == "FAIL"


def test_declared_unstranded_on_stranded_library_warns():
    s = infer_strandedness(JACKSON, declared="unstranded")
    assert s["verdict"] == "WARN"


def test_ncrna_fractions_per_sample_and_class():
    # Same total, different cause: poor rRNA depletion (a) vs abundant Ms1 (b).
    counts = pd.DataFrame({"a": [90, 10, 0], "b": [10, 10, 80]},
                          index=["MABr5051", "MAB0001", "MABnc_ms1"])
    f = ncrna_fractions(counts, {"MABr5051": "rRNA", "MABnc_ms1": "Ms1_RNA"})
    assert f["a"]["total"] == f["b"]["total"] == 0.9
    assert f["a"]["by_class"] == {"Ms1_RNA": 0.0, "rRNA": 0.9}
    assert f["b"]["by_class"] == {"Ms1_RNA": 0.8, "rRNA": 0.1}


# Boulder, full Jackson 1239 run (2026-09-18), assigned fraction before the
# structural RNAs were counted: -s 2 vs -s 1 was 51-55% vs 3.4-4.0% (7H9/SCFM2)
# and 17-20% vs 1.6-1.9% (CF sputum). Both are correctly reverse-stranded.
def test_full_run_numbers_confirm_reverse():
    for rev, fwd, uns in ((0.51, 0.040, 0.53), (0.17, 0.019, 0.18)):
        s = infer_strandedness({"reverse": rev, "forward": fwd, "unstranded": uns}, "reverse")
        assert s["inferred"] == "reverse" and s["verdict"] == "PASS"
        v = triage_sample(96.0, rev, 0.05, strand=s)
        assert v["verdict"] == "WARN"                       # not FAIL, no strand advice
        assert not any("strandedness:" in r for r in v["reasons"])


def test_triage_pass():
    v = triage_sample(96.5, 0.72, 0.12, strand=infer_strandedness(JACKSON, "reverse"))
    assert v["verdict"] == "PASS"


def test_triage_fail_low_alignment():
    v = triage_sample(85.0, 0.72, 0.12)
    assert v["verdict"] == "FAIL"


def test_low_assigned_with_confirmed_strand_warns_not_strand_advice():
    v = triage_sample(97.0, 0.17, 0.10, strand=infer_strandedness(JACKSON, "reverse"),
                      nofeature_frac=0.80)
    assert v["verdict"] == "WARN"
    text = " ".join(v["reasons"]).lower()
    assert "confirmed" in text and "80%" in text
    assert "wrong strand" not in text


def test_wrong_strand_fails_with_advice():
    v = triage_sample(97.0, 0.02, 0.10, strand=infer_strandedness(JACKSON, "forward"))
    assert v["verdict"] == "FAIL"
    assert any("strandedness: reverse" in r for r in v["reasons"])


def test_high_ncrna_warns():
    v = triage_sample(97.0, 0.9, 0.9, strand=infer_strandedness(JACKSON, "reverse"))
    assert v["verdict"] == "WARN" and any("ncRNA" in r for r in v["reasons"])
