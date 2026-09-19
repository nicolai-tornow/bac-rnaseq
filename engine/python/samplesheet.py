from __future__ import annotations
import csv
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

REQUIRED = {"sample_id", "fastq_r1", "condition"}
# Sample IDs become file names, featureCounts column names and R column names.
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


@dataclass
class Sample:
    sample_id: str
    fastq_r1: str
    condition: str
    fastq_r2: Optional[str] = None
    replicate: Optional[str] = None
    batch: Optional[str] = None


def read_samplesheet(path) -> list[Sample]:
    p = Path(path)
    delim = "," if p.suffix.lower() == ".csv" else "\t"
    rows = list(csv.DictReader(p.read_text().splitlines(), delimiter=delim))
    if not rows:
        raise ValueError("empty sample sheet")
    missing = REQUIRED - set(rows[0].keys())
    if missing:
        raise ValueError(f"sample sheet missing columns: {sorted(missing)}")
    out, seen = [], set()
    for r in rows:
        sid = r["sample_id"].strip()
        if sid in seen:
            raise ValueError(f"duplicate sample_id: {sid}")
        seen.add(sid)
        r2 = (r.get("fastq_r2") or "").strip() or None
        out.append(Sample(sid, r["fastq_r1"].strip(), r["condition"].strip(),
                           r2, (r.get("replicate") or "").strip() or None,
                           (r.get("batch") or "").strip() or None))
    return out


def apply_layout(samples, layout: str = "auto"):
    """mate1_only: drop every fastq_r2 so all samples run single-end from read 1."""
    if layout == "mate1_only":
        return [replace(s, fastq_r2=None) for s in samples]
    return samples


def is_paired(samples) -> bool:
    flags = {s.fastq_r2 is not None for s in samples}
    if len(flags) > 1:
        se = [s.sample_id for s in samples if s.fastq_r2 is None]
        raise ValueError(
            "mixed single-end and paired-end samples are not supported "
            f"(single-end: {', '.join(se)}). Set `reads: {{layout: mate1_only}}` in the "
            "config to run every sample single-end from fastq_r1, so layout is not "
            "confounded with condition.")
    return flags == {True}


def check_samplesheet(samples, contrasts=(), check_files=True) -> list[str]:
    """Problems that would otherwise surface late (in R, or mid-run). Empty = ok."""
    problems = []
    for s in samples:
        if not SAFE_ID.match(s.sample_id):
            problems.append(f"sample_id '{s.sample_id}' must match {SAFE_ID.pattern} "
                            "(used as a file name and a column name)")
        if check_files:
            for f in (s.fastq_r1, s.fastq_r2):
                if f and not Path(f).is_file():
                    problems.append(f"{s.sample_id}: FASTQ not found: {f}")
    levels = {s.condition for s in samples}
    for c in contrasts:
        for lvl in (c.numerator, c.denominator):
            if lvl not in levels:
                problems.append(f"contrast '{c.name}': condition '{lvl}' is not in the "
                                f"sample sheet (have: {', '.join(sorted(levels))})")
        if not SAFE_ID.match(c.name):
            problems.append(f"contrast name '{c.name}' must match {SAFE_ID.pattern} "
                            "(used as a file name)")
    return problems
