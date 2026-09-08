from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from . import config as V
from .selection import resolve


def load_results(path):
    df = pd.read_csv(path, sep="\t")
    df = df.dropna(subset=[V.PADJ_COL]).set_index(V.GENE_COL)
    return df


def _cloud_color(row):
    if row[V.PADJ_COL] <= V.PADJ and row[V.LOG2FC_COL] >= V.LOG2FC:
        return V.SIG_UP_MUTED
    if row[V.PADJ_COL] <= V.PADJ and row[V.LOG2FC_COL] <= -V.LOG2FC:
        return V.SIG_DOWN_MUTED
    return V.NONSIG


def _neglog(padj):
    return -np.log10(padj.replace(0, np.nextafter(0, 1)))


def volcano(df, out_stem, selection=None, tpm=None, label=True):
    V.apply_style()
    x = df[V.LOG2FC_COL].values
    y = _neglog(df[V.PADJ_COL]).values
    colors = df.apply(_cloud_color, axis=1).values
    fig, ax = plt.subplots(figsize=V.FIGSIZE)
    ax.scatter(x, y, s=5, c=colors, linewidths=0, rasterized=True)
    lim = float(np.nanmax(np.abs(x))) + 0.5
    ax.set_xlim(-lim, lim)
    ax.axhline(-np.log10(V.PADJ), ls="--", lw=0.5, color="k")
    ax.axvline(V.LOG2FC, ls="--", lw=0.5, color="k")
    ax.axvline(-V.LOG2FC, ls="--", lw=0.5, color="k")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(padj)")
    ax.set_box_aspect(1)
    if selection:
        genes = resolve(df, selection, tpm=tpm)
        sub = df.loc[genes]
        ax.scatter(sub[V.LOG2FC_COL], _neglog(sub[V.PADJ_COL]), s=22, facecolors="none",
                   edgecolors=[V.SIG_UP if v >= 0 else V.SIG_DOWN for v in sub[V.LOG2FC_COL]],
                   linewidths=0.6)
        if label:
            for g, r in sub.iterrows():
                ax.annotate(str(g), (r[V.LOG2FC_COL], _neglog(pd.Series([r[V.PADJ_COL]])).iloc[0]),
                            fontsize=5, ha="center", va="bottom")
    return V.save(fig, out_stem)
