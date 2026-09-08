import glob
import gzip
import shutil
import subprocess
import pytest
from pathlib import Path
from engine.python.build_refs import build_bundle
from engine.python import commands as C

REPO = Path(__file__).resolve().parents[2]
READS = sorted(glob.glob("/path/to/projects/rnaseq_mabs_media/reads/*_R1_001.fastq.gz"))
have_tools = all(shutil.which(t) for t in ["bowtie2", "samtools", "featureCounts"])
pytestmark = pytest.mark.skipif(not (READS and have_tools),
                                reason="need real Mabs reads + alignment tools")


def _subsample(src, dst, n_reads=50000):
    with gzip.open(src, "rt") as fh, gzip.open(dst, "wt") as out:
        for i, line in enumerate(fh):
            if i >= n_reads * 4:
                break
            out.write(line)


def _assigned(summary):
    tot = asg = 0
    for line in Path(summary).read_text().splitlines()[1:]:
        s, *v = line.split("\t")
        n = sum(map(int, v))
        tot += n
        if s == "Assigned":
            asg = n
    return asg / tot if tot else 0.0


def test_reverse_strand_wins_on_real_reads(tmp_path):
    r1 = READS[0]
    r2 = r1.replace("_R1_001", "_R2_001")
    assert Path(r2).exists(), "expected paired R2"
    s1, s2 = tmp_path / "r1.fq.gz", tmp_path / "r2.fq.gz"
    _subsample(r1, s1)
    _subsample(r2, s2)

    b = build_bundle("mabs", refs_root=REPO / "refs", out_dir=tmp_path / "ref", threads=4)
    sam = tmp_path / "a.sam"
    with open(sam, "w") as out:
        subprocess.run(C.bowtie2_cmd(b["index_prefix"], str(s1), 4, r2=str(s2)),
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    bam = tmp_path / "a.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)

    fracs = {}
    for strand in ("unstranded", "forward", "reverse"):
        outp = tmp_path / f"fc_{strand}.txt"
        subprocess.run(C.featurecounts_cmd(b["saf"], str(outp), [str(bam)], 4,
                                           strandedness=strand, paired=True),
                       capture_output=True, check=True)
        fracs[strand] = _assigned(str(outp) + ".summary")

    # Reverse-stranded libraries assign the most reads at -s 2 (proves strand-correct counting)
    assert fracs["reverse"] > fracs["forward"]
    assert fracs["reverse"] > fracs["unstranded"]
