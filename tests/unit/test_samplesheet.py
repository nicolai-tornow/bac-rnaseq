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


MIXED = ("sample_id\tfastq_r1\tfastq_r2\tcondition\n"
         "a\tA.fq.gz\t\t7H9\nb\tB_R1.fq.gz\tB_R2.fq.gz\tSCFM2\n")


def test_mixed_layout_error_points_to_mate1_only(tmp_path):
    s = read_samplesheet(_w(tmp_path, MIXED))
    with pytest.raises(ValueError, match="mate1_only"):
        is_paired(s)


def test_mate1_only_makes_all_single_end(tmp_path):
    from engine.python.samplesheet import apply_layout
    s = apply_layout(read_samplesheet(_w(tmp_path, MIXED)), "mate1_only")
    assert is_paired(s) is False
    assert [x.fastq_r1 for x in s] == ["A.fq.gz", "B_R1.fq.gz"]


def test_check_samplesheet_catches_late_failures(tmp_path):
    from engine.python.samplesheet import check_samplesheet
    from engine.python.config import Contrast
    r1 = tmp_path / "A.fq.gz"
    r1.write_text("")
    txt = ("sample_id\tfastq_r1\tcondition\n"
           f"ok_1\t{r1}\t7H9\nbad id\t{tmp_path}/missing.fq.gz\tSCFM2\n")
    s = read_samplesheet(_w(tmp_path, txt))
    probs = check_samplesheet(s, [Contrast(name="S_vs_7", numerator="SCFM", denominator="7H9")])
    joined = "\n".join(probs)
    assert "bad id" in joined                      # unsafe sample_id
    assert "missing.fq.gz" in joined               # FASTQ does not exist
    assert "condition 'SCFM'" in joined            # contrast level typo
    assert "ok_1" not in joined
