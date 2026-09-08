from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENE_COL, LOG2FC_COL, PADJ_COL = "Gene", "log2FoldChange", "padj"
PADJ, LOG2FC = 0.05, 1.0
SIG_UP, SIG_DOWN = "#d52928", "#1f77b4"
SIG_UP_MUTED, SIG_DOWN_MUTED, NONSIG = "#e8a3a3", "#a3c4db", "#cccccc"
FIGSIZE = (3.7, 3.5)
DPI = 300


def apply_style():
    matplotlib.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial",
        "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6,
    })


def save(fig, out_stem):
    pdf, png = f"{out_stem}.pdf", f"{out_stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return pdf, png
