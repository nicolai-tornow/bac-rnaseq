import pandas as pd
from engine.viz.selection import resolve


def _df():
    return pd.DataFrame({"log2FoldChange": [3.0, -2.0, 0.1, 5.0], "padj": [0.001, 0.01, 0.9, 1e-6]},
                        index=["MAB0001", "MAB0002c", "MAB0003", "MAB0004"])


def test_explicit_genes_case_insensitive():
    assert set(resolve(_df(), {"mode": "genes", "genes": ["mab0001", "MAB0004"]})) == {"MAB0001", "MAB0004"}


def test_top_by_l2fc():
    assert resolve(_df(), {"mode": "top", "n": 2, "by": "l2fc"}) == ["MAB0004", "MAB0001"]


def test_bottom_by_l2fc():
    assert resolve(_df(), {"mode": "bottom", "n": 1, "by": "l2fc"}) == ["MAB0002c"]
