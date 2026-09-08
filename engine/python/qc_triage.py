from __future__ import annotations
import re
from pathlib import Path


def parse_bowtie2_log(text: str) -> float:
    m = re.search(r"([\d.]+)%\s+overall alignment rate", text)
    if not m:
        raise ValueError("no overall alignment rate in bowtie2 log")
    return float(m.group(1))


def parse_featurecounts_summary(path) -> dict:
    total = assigned = 0
    for line in Path(path).read_text().splitlines()[1:]:
        status, *vals = line.split("\t")
        n = sum(int(v) for v in vals)
        total += n
        if status == "Assigned":
            assigned = n
    return {"assigned": assigned, "total": total,
            "assigned_frac": (assigned / total) if total else 0.0}


def ncrna_fraction(counts_df, ncrna_ids) -> float:
    total = counts_df.sum().sum()
    if total == 0:
        return 0.0
    nc = counts_df.loc[counts_df.index.intersection(ncrna_ids)].sum().sum()
    return float(nc / total)


def triage_sample(align_pct, assigned_frac, ncrna_frac, dup_pct=None,
                  strandedness="reverse") -> dict:
    reasons, verdict = [], "PASS"

    def worse(v):
        order = {"PASS": 0, "WARN": 1, "FAIL": 2}
        return v if order[v] > order[verdict] else verdict

    if align_pct < 90:
        verdict = worse("FAIL"); reasons.append(f"alignment {align_pct:.1f}% < 90%")
    elif align_pct < 95:
        verdict = worse("WARN"); reasons.append(f"alignment {align_pct:.1f}% < 95%")
    if assigned_frac < 0.40:
        verdict = worse("FAIL")
        reasons.append(f"assigned {assigned_frac:.0%} < 40% — possible wrong strandedness (declared {strandedness})")
    elif assigned_frac < 0.60:
        verdict = worse("WARN"); reasons.append(f"assigned {assigned_frac:.0%} < 60%")
    if ncrna_frac > 0.85:
        verdict = worse("WARN"); reasons.append(f"ncRNA {ncrna_frac:.0%} > 85% — low usable mRNA depth")
    return {"verdict": verdict, "reasons": reasons}
