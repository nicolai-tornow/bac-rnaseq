---
name: batch-integration
description: Use to check and visualize batch effects across two or more RNAseq runs (ComBat on VST for PCA/correlation only). Differential expression stays per-batch; this does NOT merge counts for DE.
---

# Batch Integration

Guardrails (enforced): inputs must share a reference (same gene universe) — the
tool aborts otherwise; ComBat is for **visualization only** (PCA/correlation),
never for DE calls; DE comparisons stay per-batch.

1. Collect the VST matrices (`06_deseq/vst_counts.tsv`) from each run and a
   metadata TSV with columns `sample, batch, condition`.
2. Run: `${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq batch --vst run1/vst.tsv run2/vst.tsv --meta meta.tsv --out-dir <dir>`.
3. Inspect PCA before/after ComBat. **If there was no strong batch effect, ComBat
   OVER-CORRECTS — trust the "before" plot.** For DE across batches, re-run DESeq2
   with `design = ~ batch + condition` rather than merging the ComBat output.
