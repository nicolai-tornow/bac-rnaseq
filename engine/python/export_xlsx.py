from __future__ import annotations
from pathlib import Path
import pandas as pd
from .ids import normalize_id


def export_workbook(results_dir, out_xlsx, gene_names=None, normalized_tsv=None,
                    vst_tsv=None, tpm_df=None):
    results_dir = Path(results_dir)
    gene_names = gene_names or {}
    # Match gene IDs case-insensitively (results may use mab0001, the GFF MAB0001).
    names = {normalize_id(str(k)): v for k, v in gene_names.items()}
    result_files = sorted(results_dir.glob("*.tsv"))
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        summary = []
        for f in result_files:
            df = pd.read_csv(f, sep="\t")
            sig = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 1)]
            summary.append({"contrast": f.stem, "n_genes": len(df), "n_sig": len(sig)})
        pd.DataFrame(summary).to_excel(xw, sheet_name="summary", index=False)
        for f in result_files:
            df = pd.read_csv(f, sep="\t")
            df.insert(1, "gene_name", [names.get(normalize_id(str(g)), g) for g in df["Gene"]])
            df.to_excel(xw, sheet_name=f.stem[:31], index=False)
        for name, path in [("normalized", normalized_tsv), ("vst", vst_tsv)]:
            if path and Path(path).exists():
                pd.read_csv(path, sep="\t").to_excel(xw, sheet_name=name, index=False)
        if tpm_df is not None:
            tpm_df.reset_index().rename(columns={"index": "Gene"}).to_excel(
                xw, sheet_name="tpm", index=False)
    return out_xlsx
