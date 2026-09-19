from __future__ import annotations
import re
from pathlib import Path

_ORDER = {"PASS": 0, "WARN": 1, "FAIL": 2}
_FLAG = {"reverse": "-s 2", "forward": "-s 1", "unstranded": "-s 0"}


def parse_bowtie2_log(text: str) -> float:
    m = re.search(r"([\d.]+)%\s+overall alignment rate", text)
    if not m:
        raise ValueError("no overall alignment rate in bowtie2 log")
    return float(m.group(1))


def sample_name(column: str) -> str:
    """featureCounts column (a BAM path) -> sample_id."""
    name = Path(column).name
    return name[:-4] if name.endswith(".bam") else name


def parse_featurecounts_summary(path) -> dict:
    """Per sample: every status count, plus total and assigned_frac."""
    lines = Path(path).read_text().splitlines()
    samples = [sample_name(c) for c in lines[0].split("\t")[1:]]
    out = {s: {} for s in samples}
    for line in lines[1:]:
        status, *vals = line.split("\t")
        for s, v in zip(samples, vals):
            out[s][status] = int(v)
    for s, d in out.items():
        total = sum(d.values())
        d["total"] = total
        d["assigned_frac"] = d.get("Assigned", 0) / total if total else 0.0
        d["nofeature_frac"] = d.get("Unassigned_NoFeatures", 0) / total if total else 0.0
    return out


def infer_strandedness(fracs: dict, declared: str) -> dict:
    """Assigned fractions counted at all three -s settings -> library strandedness.

    A stranded library assigns far more reads on its own strand than on the other
    (ratio >= 5). Unstranded counting also collects antisense reads, so it is
    usually slightly HIGHER than the correct stranded setting; that is expected and
    is not evidence against the stranded setting.
    """
    rev, fwd, uns = fracs["reverse"], fracs["forward"], fracs["unstranded"]
    lo, hi = min(rev, fwd), max(rev, fwd)
    if hi < 0.01:
        inferred = "ambiguous"
    elif rev >= 5 * max(fwd, 1e-9):
        inferred = "reverse"
    elif fwd >= 5 * max(rev, 1e-9):
        inferred = "forward"
    elif hi <= 1.5 * lo:
        inferred = "unstranded"
    else:
        inferred = "ambiguous"

    shown = ", ".join(f"{k} ({_FLAG[k]}) {v:.1%}" for k, v in
                      (("reverse", rev), ("forward", fwd), ("unstranded", uns)))
    if inferred == declared:
        verdict, msg = "PASS", f"strandedness confirmed as {declared}: {shown}"
    elif inferred == "ambiguous":
        verdict, msg = "WARN", f"strandedness could not be inferred: {shown}"
    elif declared == "unstranded":
        verdict = "WARN"
        msg = (f"library looks {inferred}-stranded but is counted unstranded "
               f"(antisense reads are counted too); consider strandedness: {inferred}. {shown}")
    else:
        verdict = "FAIL"
        msg = (f"declared strandedness {declared} does not match the library, which "
               f"looks {inferred}: {shown}. Re-run with strandedness: {inferred}")
    return {"declared": declared, "inferred": inferred, "fractions": fracs,
            "verdict": verdict, "message": msg}


def ncrna_fractions(counts_df, ncrna_classes) -> dict:
    """Per sample: share of assigned reads on structural RNAs, in total and by class.

    ncrna_classes maps gene ID -> class (rRNA, tRNA, tmRNA, Ms1_RNA, ...). Keeping the
    classes apart matters: rRNA+tRNA measures depletion efficiency, whereas Ms1,
    tmRNA and RNase P RNA are abundant by biology and can dominate a library.
    """
    if not isinstance(ncrna_classes, dict):
        ncrna_classes = {i: "ncRNA" for i in ncrna_classes}
    ids = counts_df.index.intersection(list(ncrna_classes))
    out = {}
    for col in counts_df.columns:
        tot = counts_df[col].sum()
        by = {}
        for gid in ids:
            c = ncrna_classes[gid]
            by[c] = by.get(c, 0) + counts_df.at[gid, col]
        out[col] = {"total": float(sum(by.values()) / tot) if tot else 0.0,
                    "by_class": {c: float(n / tot) if tot else 0.0 for c, n in sorted(by.items())}}
    return out


def triage_sample(align_pct, assigned_frac, ncrna_frac, strand=None,
                  nofeature_frac=None) -> dict:
    reasons, verdict = [], "PASS"

    def worse(v):
        return v if _ORDER[v] > _ORDER[verdict] else verdict

    if align_pct < 90:
        verdict = worse("FAIL"); reasons.append(f"alignment {align_pct:.1f}% < 90%")
    elif align_pct < 95:
        verdict = worse("WARN"); reasons.append(f"alignment {align_pct:.1f}% < 95%")

    strand_ok = strand is not None and strand["verdict"] == "PASS"
    if strand is not None and strand["verdict"] != "PASS":
        verdict = worse(strand["verdict"]); reasons.append(strand["message"])

    if assigned_frac < 0.60:
        level = "FAIL" if (assigned_frac < 0.40 and strand is None) else "WARN"
        nf = f"; {nofeature_frac:.0%} of reads Unassigned_NoFeatures" if nofeature_frac is not None else ""
        if strand_ok:
            why = ("strandedness confirmed, so these reads fall outside the annotation "
                   "(unannotated RNA, contamination, or an incomplete GFF)")
        elif strand is None:
            why = "strandedness not assessed"
        else:
            why = "see the strandedness message"
        verdict = worse(level)
        reasons.append(f"assigned {assigned_frac:.0%} < {'40' if assigned_frac < 0.40 else '60'}%{nf} — {why}")

    if ncrna_frac > 0.85:
        verdict = worse("WARN")
        reasons.append(f"ncRNA {ncrna_frac:.0%} of assigned reads > 85% — low usable mRNA depth")
    return {"verdict": verdict, "reasons": reasons}
