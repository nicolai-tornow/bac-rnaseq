---
name: visualize-results
description: Use to make figures from DESeq2 result tables — volcano plots, TPM bar/heatmap, and expression rankings. Requires a gene selection input, one of a pathway/operon locus-tag list, explicit locus tags, or top-N/bottom-N genes (by fold change or TPM).
---

# Visualize Results

Consumes any DESeq2 result TSV (`Gene, baseMean, log2FoldChange, lfcSE, stat,
pvalue, padj`) — a full run is NOT required. Tool names per
`skills/_shared/references/<harness>-tools.md`. Figures are written as PDF + PNG
in the figure-1 house style (Arial, muted red/blue/grey, dashed thresholds at
padj 0.05 and |log2FC| 1, square axes).

## Ask the user which genes to feature (required)

One of:
- **pathway / operon**: a list of locus tags for the pathway or operon of interest.
- **specific locus tags**: an explicit list (matched case-insensitively).
- **top-N / bottom-N**: the N most up-/down-regulated genes, ranked by log2 fold
  change or by TPM.

## Produce the figures

- **Volcano** (`${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq visualize <results.tsv> --out <stem> --top 10`
  or `--genes MAB0001,MAB0002,...`): the significance cloud plus the selected
  genes highlighted and labeled.
- **TPM bar / heatmap / ranking** (call `engine.viz.tpm` / `engine.viz.ranking`
  with counts + the SAF for gene lengths) for absolute-expression views.

Report the output file paths.
