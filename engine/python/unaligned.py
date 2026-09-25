"""What the reads of a poorly aligning sample are: count, GC% and the fraction that
looks like rRNA. Reads of the species' own rRNA already aligned, so the rRNA-like ones
are rRNA of something else (contamination); they are found by local alignment to the
reference's own 16S/23S/5S genes, which are conserved enough across bacteria."""
from __future__ import annotations
import gzip
import os
from pathlib import Path
from . import qc_triage


def gc_and_count(paths) -> tuple[int, float]:
    n = gc = acgt = 0
    for p in paths:
        with gzip.open(p, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    s = line.strip().upper()
                    n += 1
                    g = s.count("G") + s.count("C")
                    gc += g
                    acgt += g + s.count("A") + s.count("T")
    return n, (100.0 * gc / acgt if acgt else 0.0)


def rrna_index(bundle, dst, runner) -> str | None:
    """bowtie2 index of the bundle's rRNA genes, or None if it annotates none."""
    ids = {g for g, c in (bundle.get("structural") or {}).items() if c == "rRNA"}
    if not ids:
        return None
    regions = []
    for line in Path(bundle["saf"]).read_text().splitlines()[1:]:
        f = line.split("\t")
        if f[0] in ids:
            regions.append(f"{f[1]}:{f[2]}-{f[3]}")
    if not regions:
        return None
    dst = Path(dst)
    dst.mkdir(parents=True, exist_ok=True)
    genome = Path(bundle["saf"]).parent / "genome.fasta"
    res = runner.run(["samtools", "faidx", str(genome), *regions])
    if res.returncode != 0:
        raise RuntimeError(f"samtools faidx of the rRNA genes failed: {res.stderr.strip()}")
    (dst / "rrna.fa").write_text(res.stdout)
    prefix = dst / "rrna"
    res = runner.run(["bowtie2-build", "--quiet", str(dst / "rrna.fa"), str(prefix)])
    if res.returncode != 0:
        raise RuntimeError(f"bowtie2-build of the rRNA genes failed: {res.stderr.strip()}")
    return str(prefix)


def rrna_like_fraction(index, paths, runner, threads) -> float:
    # A 15-nt seed with one mismatch finds ~100% of reads at 8-12% divergence from the
    # reference rRNA (--very-sensitive-local alone: 81% / 51%) and still 0% of random
    # sequence.
    res = runner.run(["bowtie2", "-x", index, "--very-sensitive-local", "-N", "1", "-L", "15",
                      "--no-unal",
                      "-p", str(threads), "-U", ",".join(map(str, paths)), "-S", os.devnull])
    if res.returncode != 0:
        raise RuntimeError(f"bowtie2 against the rRNA genes failed: {res.stderr.strip()}")
    return round(qc_triage.parse_bowtie2_log(res.stderr) / 100, 4)


def diagnose(sid, out, bundle, runner, threads, cache) -> dict | None:
    """Diagnostics from qc/unaligned/<sid>/, or None when nothing was captured."""
    d = Path(out) / "qc" / "unaligned" / sid
    files = sorted(str(p) for p in d.glob("*.fq.gz")) if d.is_dir() else []
    if not files:
        return None
    n, gc = gc_and_count(files)
    if "index" not in cache:
        cache["index"] = rrna_index(bundle, Path(out) / "qc" / "unaligned" / "rrna_index", runner)
    res = {"reads": n, "gc_pct": round(gc, 2), "rrna_like_frac": None,
           "files": [os.path.relpath(f, out) for f in files]}
    if cache["index"] is None:
        res["note"] = ("the reference annotates no rRNA genes, so the rRNA-like fraction "
                       "was not measured")
    elif n:
        res["rrna_like_frac"] = rrna_like_fraction(cache["index"], files, runner, threads)
    return res
