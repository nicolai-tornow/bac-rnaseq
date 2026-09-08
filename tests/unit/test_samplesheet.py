import pytest
from engine.python.samplesheet import read_samplesheet, is_paired

PE = "sample_id\tfastq_r1\tfastq_r2\tcondition\n7H9_rep1\tA_R1.fq.gz\tA_R2.fq.gz\t7H9\nS_rep1\tB_R1.fq.gz\tB_R2.fq.gz\tSCFM2\n"
MISSING = "sample_id\tfastq_r1\n7H9_rep1\tA_R1.fq.gz\n"
DUP = "sample_id\tfastq_r1\tcondition\nx\tA.fq.gz\t7H9\nx\tB.fq.gz\tSCFM2\n"


def _w(tmp_path, txt, name="s.tsv"):
    p = tmp_path / name
    p.write_text(txt)
    return str(p)


def test_reads_paired(tmp_path):
    s = read_samplesheet(_w(tmp_path, PE))
    assert len(s) == 2 and s[0].sample_id == "7H9_rep1" and s[0].condition == "7H9"
    assert is_paired(s) is True


def test_missing_condition_raises(tmp_path):
    with pytest.raises(ValueError):
        read_samplesheet(_w(tmp_path, MISSING))


def test_duplicate_sample_id_raises(tmp_path):
    with pytest.raises(ValueError):
        read_samplesheet(_w(tmp_path, DUP))
