from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from . import config as V


def compute_tpm(counts_df, lengths, ncrna_ids=()):
    drop = set(ncrna_ids)
    genes = [g for g in counts_df.index if g not in drop and g in lengths]
    cnt = counts_df.loc[genes]
    kb = pd.Series({g: lengths[g] / 1000.0 for g in genes})
    rpk = cnt.div(kb, axis=0)
    return rpk.div(rpk.sum(axis=0) / 1e6, axis=1)


def load_counts(counts_tsv):
    return pd.read_csv(counts_tsv, sep="\t", index_col=0)


def lengths_from_saf(saf_path):
    df = pd.read_csv(saf_path, sep="\t")
    return {r.GeneID: int(r.End) - int(r.Start) + 1 for r in df.itertuples(index=False)}


def heatmap(tpm_df, genes, out_stem):
    V.apply_style()
    keep = [g for g in genes if g in tpm_df.index]
    sub = np.log10(tpm_df.loc[keep] + 1)
    fig, ax = plt.subplots(figsize=(0.5 * sub.shape[1] + 2, 0.25 * sub.shape[0] + 1))
    im = ax.imshow(sub.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(sub.shape[1])); ax.set_xticklabels(sub.columns, rotation=90)
    ax.set_yticks(range(sub.shape[0])); ax.set_yticklabels(sub.index)
    fig.colorbar(im, ax=ax, label="log10(TPM+1)")
    return V.save(fig, out_stem)


def bar(tpm_df, genes, out_stem):
    V.apply_style()
    keep = [g for g in genes if g in tpm_df.index]
    sub = tpm_df.loc[keep].mean(axis=1).sort_values()
    fig, ax = plt.subplots(figsize=V.FIGSIZE)
    ax.barh(range(len(sub)), sub.values, color=V.SIG_DOWN)
    ax.set_yticks(range(len(sub))); ax.set_yticklabels(sub.index)
    ax.set_xlabel("mean TPM")
    return V.save(fig, out_stem)
