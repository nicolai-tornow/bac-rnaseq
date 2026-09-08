from __future__ import annotations
from ..python.ids import normalize_id


def resolve(df, selection, tpm=None):
    mode = selection["mode"]
    if mode == "genes":
        want = {normalize_id(g) for g in selection["genes"]}
        return [g for g in df.index if normalize_id(g) in want]
    n = int(selection.get("n", 10))
    by = selection.get("by", "l2fc")
    if by == "tpm":
        if tpm is None:
            raise ValueError("selection by tpm requires a tpm series")
        ranked = tpm.reindex(df.index).dropna().sort_values(ascending=False)
    else:
        ranked = df["log2FoldChange"].sort_values(ascending=False)
    return list(ranked.head(n).index) if mode == "top" else list(ranked.tail(n).index[::-1])
