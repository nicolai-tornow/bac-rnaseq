from __future__ import annotations
import pandas as pd
from ..python.ids import normalize_id, build_canonical_map


def load_genesets(xlsx, sheet, id_col, cat_col):
    df = pd.read_excel(xlsx, sheet_name=sheet)
    sets = {}
    for gid, cat in zip(df[id_col], df[cat_col]):
        if pd.isna(gid) or pd.isna(cat):
            continue
        sets.setdefault(str(cat), set()).add(str(gid))
    return sets


def align_hits(de_ids, geneset_ids, max_unmatched=0.05):
    canon = build_canonical_map(geneset_ids)
    mapping, unmatched = {}, 0
    for d in de_ids:
        ref = canon.get(normalize_id(d))
        if ref is None:
            unmatched += 1
        else:
            mapping[d] = ref
    frac = unmatched / len(de_ids) if de_ids else 0.0
    if frac > max_unmatched:
        raise ValueError(f"{frac:.1%} of DE IDs unmatched to category sheet (>{max_unmatched:.0%}) "
                         "— check ID case/spelling before trusting enrichment")
    return mapping
