import shutil
import subprocess
import pytest
from engine.python.build_refs import build_bundle
from engine.python import commands as C
from engine.python.qc_triage import parse_featurecounts_summary, infer_strandedness
from tests.testdata import REFS, READS_R1, READS_R2

have_tools = all(shutil.which(t) for t in ["bowtie2", "samtools", "featureCounts"])
pytestmark = pytest.mark.skipif(not have_tools, reason="need bowtie2, samtools, featureCounts")


def test_reverse_strand_wins_on_real_reads(tmp_path):
    # 20k read pairs of a public reverse-stranded library (SRR10958838, KAPA HyperPrep).
    b = build_bundle("mabs", refs_root=REFS, out_dir=tmp_path / "ref", threads=4)
    sam, bam = tmp_path / "a.sam", tmp_path / "a.bam"
    with open(sam, "w") as out:
        subprocess.run(C.bowtie2_cmd(b["index_prefix"], str(READS_R1), 4, r2=str(READS_R2)),
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)

    fracs = {}
    for strand in ("unstranded", "forward", "reverse"):
        outp = tmp_path / f"fc_{strand}.txt"
        subprocess.run(C.featurecounts_cmd(b["saf"], str(outp), [str(bam)], 4,
                                           strandedness=strand, paired=True),
                       capture_output=True, check=True)
        fracs[strand] = parse_featurecounts_summary(str(outp) + ".summary")["a"]["assigned_frac"]

    # A reverse library assigns far more reads at -s 2 than at -s 1. Unstranded
    # counting also collects antisense reads, so it is usually slightly HIGHER than
    # -s 2 on a correct library: require only that -s 2 is not far below it.
    assert fracs["reverse"] >= 5 * fracs["forward"]
    assert fracs["reverse"] >= 0.9 * fracs["unstranded"]
    assert infer_strandedness(fracs, "reverse")["verdict"] == "PASS"
    assert infer_strandedness(fracs, "forward")["verdict"] == "FAIL"
