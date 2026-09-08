from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def run_ora(hits, gene_sets, universe, alternative="two-sided", min_size=2, max_size=1000):
    universe = set(universe)
    hits = set(hits) & universe
    N, K = len(universe), len(hits)
    rows = []
    for name, s in gene_sets.items():
        s = set(s) & universe
        m = len(s)
        if m < min_size or m > max_size:
            continue
        a = len(s & hits)
        table = [[a, K - a], [m - a, N - K - (m - a)]]
        odds, p = fisher_exact(table, alternative=alternative)
        if odds == np.inf:
            log2_or = np.inf
        elif odds == 0.0:
            log2_or = -np.inf
        else:
            log2_or = float(np.log2(odds))
        rows.append({"set": name, "size": m, "n_hits": a, "odds_ratio": odds,
                     "log2_or": log2_or, "p": p, "effect": "over" if odds > 1 else "under"})
    df = pd.DataFrame(rows)
    if not df.empty:
        df["fdr"] = multipletests(df["p"].values, method="fdr_bh")[1]
        df = df.sort_values("p").reset_index(drop=True)
    return df
