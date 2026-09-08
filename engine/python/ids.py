from __future__ import annotations


def normalize_id(s: str) -> str:
    return s.strip().replace("MAB_", "MAB").replace("mab_", "MAB").upper()


def build_canonical_map(ids) -> dict:
    return {normalize_id(x): x for x in ids}


def match_ids(query_ids, ref_ids):
    canon = build_canonical_map(ref_ids)
    mapping, unmatched = {}, 0
    for q in query_ids:
        ref = canon.get(normalize_id(q))
        if ref is None:
            unmatched += 1
        else:
            mapping[q] = ref
    frac = unmatched / len(query_ids) if query_ids else 0.0
    return mapping, frac
