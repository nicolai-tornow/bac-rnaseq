from engine.python import commands as c


def test_strand_flag():
    assert c.strand_flag("reverse") == "2"
    assert c.strand_flag("forward") == "1"
    assert c.strand_flag("unstranded") == "0"


def test_fastp_has_validated_params():
    cmd = c.fastp_cmd("a_R1.fq.gz", "o1.fq.gz", 4, r2="a_R2.fq.gz", out2="o2.fq.gz",
                      json="a.json", html="a.html")
    assert "--detect_adapter_for_pe" in cmd
    assert cmd[cmd.index("--qualified_quality_phred") + 1] == "20"
    assert cmd[cmd.index("--length_required") + 1] == "36"
    assert "--in2" in cmd and "--out2" in cmd


def test_fastp_single_end_has_no_r2():
    cmd = c.fastp_cmd("a_R1.fq.gz", "o1.fq.gz", 4)
    assert "--in2" not in cmd


def test_bowtie2_sensitive_no_unal_paired():
    cmd = c.bowtie2_cmd("ref", "r1.fq.gz", 8, r2="r2.fq.gz")
    assert "--sensitive" in cmd and "--no-unal" in cmd
    assert "-1" in cmd and "-2" in cmd


def test_bowtie2_single_end_uses_U():
    cmd = c.bowtie2_cmd("ref", "r1.fq.gz", 8)
    assert "-U" in cmd and "-1" not in cmd


def test_featurecounts_reverse_paired():
    cmd = c.featurecounts_cmd("labels.saf", "fc.txt", ["a.bam", "b.bam"], 16,
                              strandedness="reverse", paired=True)
    assert cmd[cmd.index("-s") + 1] == "2"
    assert "-p" in cmd and "--countReadPairs" in cmd
    assert cmd[cmd.index("-F") + 1] == "SAF"
    assert cmd[-2:] == ["a.bam", "b.bam"]


def test_featurecounts_single_end_omits_paired_flags():
    cmd = c.featurecounts_cmd("labels.saf", "fc.txt", ["a.bam"], 4,
                              strandedness="unstranded", paired=False)
    assert "-p" not in cmd and "--countReadPairs" not in cmd
    assert cmd[cmd.index("-s") + 1] == "0"
