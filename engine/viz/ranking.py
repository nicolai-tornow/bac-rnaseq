from __future__ import annotations
import matplotlib.pyplot as plt
from . import config as V


def top_bottom(tpm_df, n):
    mean = tpm_df.mean(axis=1).sort_values(ascending=False)
    return list(mean.head(n).index), list(mean.tail(n).index)


def ranking(tpm_df, n, out_stem):
    V.apply_style()
    top, bottom = top_bottom(tpm_df, n)
    mean = tpm_df.mean(axis=1)
    genes = bottom[::-1] + top[::-1]
    vals = [mean[g] for g in genes]
    colors = [V.SIG_DOWN] * len(bottom) + [V.SIG_UP] * len(top)
    fig, ax = plt.subplots(figsize=(V.FIGSIZE[0], 0.2 * len(genes) + 1))
    ax.barh(range(len(genes)), vals, color=colors)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes)
    ax.set_xlabel("mean TPM")
    return V.save(fig, out_stem)
