"""rRNA-like fraction: unaligned reads aligned locally to the species' own rRNA genes.
Reads from M. abscessus 16S/23S with 8% substitutions stand in for another species'
rRNA; random reads for everything else."""
import gzip
import os
import random
import shutil
import pytest
from engine.python import unaligned as U
from engine.python.procs import ProcRunner
from tests.testdata import REFS

pytestmark = pytest.mark.skipif(not all(shutil.which(t) for t in ("bowtie2", "samtools")),
                                reason="need bowtie2 and samtools")
RRNA = [l.split("\t") for l in (REFS / "mabs/structural_rna.tsv").read_text().splitlines()[1:]
        if l.split("\t")[5] == "rRNA"]


def _genome():
    seq = []
    for line in open(REFS / "mabs/NC_010397.1.fasta"):
        if line.startswith(">"):
            if seq:
                break
            continue
        seq.append(line.strip())
    return "".join(seq)


def test_rrna_like_fraction_on_mutated_rrna(tmp_path):
    tmp_path = tmp_path / "a b,c"          # space and comma: -U takes a comma list
    bundle_dir = tmp_path / "bundle"
    bundle_dir.mkdir(parents=True)
    os.symlink(REFS / "mabs/NC_010397.1.fasta", bundle_dir / "genome.fasta")
    saf = bundle_dir / "labels.saf"
    saf.write_text("GeneID\tChr\tStart\tEnd\tStrand\n" +
                   "".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\n" for r in RRNA))
    bundle = {"saf": str(saf), "structural": {r[0]: "rRNA" for r in RRNA}}

    rng, g = random.Random(7), _genome()
    rrna = "".join(g[int(r[2]) - 1:int(r[3])] for r in RRNA[:2])          # 16S + 23S
    reads = []
    for _ in range(200):
        i = rng.randrange(len(rrna) - 60)
        reads.append("".join(b if rng.random() > 0.08 else rng.choice("ACGT".replace(b, ""))
                             for b in rrna[i:i + 60]))
    reads += ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(200)]
    fq = tmp_path / "unaligned.fq.gz"
    with gzip.open(fq, "wt") as fh:
        for n, r in enumerate(reads):
            fh.write(f"@r{n}\n{r}\n+\n{'I' * 60}\n")

    runner = ProcRunner()
    index = U.rrna_index(bundle, tmp_path / "idx", runner)
    frac = U.rrna_like_fraction(index, [str(fq)], runner, 2)
    assert 0.45 <= frac <= 0.55, frac                 # half the reads are rRNA
    n, gc = U.gc_and_count([str(fq)])
    assert n == 400 and 45 < gc < 65
