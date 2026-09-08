from __future__ import annotations
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

REQUIRED = {"sample_id", "fastq_r1", "condition"}


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


def is_paired(samples) -> bool:
    flags = {s.fastq_r2 is not None for s in samples}
    if len(flags) > 1:
        raise ValueError("mixed single-end and paired-end samples are not supported")
    return flags.pop()
