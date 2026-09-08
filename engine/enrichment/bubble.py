from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from ..viz import config as V


def bubble(ora_df, out_stem, top=25, p_col="p", p_thresh=0.05):
    V.apply_style()
    df = ora_df[(ora_df["effect"] == "over") & (ora_df[p_col] < p_thresh)].copy()
    df = df.sort_values(p_col).head(top)
    if df.empty:
        fig, ax = plt.subplots(figsize=V.FIGSIZE)
        ax.text(0.5, 0.5, "no over-represented sets", ha="center", va="center", fontsize=7)
        ax.axis("off")
        return V.save(fig, out_stem)
    finite = df["log2_or"].replace([np.inf, -np.inf], np.nan)
    mx = float(np.nanmax(np.abs(finite))) or 1.0
    df = df.iloc[::-1]  # smallest p at top
    mag = df["log2_or"].replace([np.inf], mx).replace([-np.inf], mx).abs().clip(upper=mx)
    sizes = 20 + 120 * (mag / mx)
    fig, ax = plt.subplots(figsize=(4.6, 0.3 * len(df) + 1.2))
    ax.scatter([0] * len(df), range(len(df)), s=sizes, color=V.SIG_UP,
               alpha=0.85, edgecolors="k", linewidths=0.3)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["set"])
    ax.set_xticks([])
    ax.set_xlim(-1, 1)
    ax.set_title("over-represented sets (bubble area ∝ |log2 OR|)", fontsize=7)
    return V.save(fig, out_stem)
