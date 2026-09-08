from engine.python.qc_triage import parse_bowtie2_log, triage_sample

BT2 = """1000 reads; of these:
  1000 (100.00%) were paired; of these:
    50 (5.00%) aligned concordantly 0 times
    900 (90.00%) aligned concordantly exactly 1 time
    50 (5.00%) aligned concordantly >1 times
96.50% overall alignment rate
"""


def test_parse_bowtie2_rate():
    assert parse_bowtie2_log(BT2) == 96.5


def test_triage_pass():
    v = triage_sample(96.5, 0.72, 0.12)
    assert v["verdict"] == "PASS"


def test_triage_fail_low_alignment():
    v = triage_sample(85.0, 0.72, 0.12)
    assert v["verdict"] == "FAIL"


def test_triage_fail_wrong_strand():
    v = triage_sample(97.0, 0.20, 0.10, strandedness="reverse")
    assert v["verdict"] == "FAIL"
    assert any("strand" in r.lower() for r in v["reasons"])
