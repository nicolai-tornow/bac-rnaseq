---
name: pathway-enrichment
description: Use to find functional categories / pathways over-represented among differentially expressed genes (Fisher over-representation, BH-FDR). For M. abscessus the category sheet is bundled; for M. tuberculosis or a custom genome, the user must supply their own category sheet.
---

# Pathway Enrichment

Consumes a DESeq2 result TSV plus a gene-category sheet. Tool names per
`skills/_shared/references/<harness>-tools.md`.

1. Significant-gene set: default `padj < 0.05` and `|log2FC| >= 1`.
2. Category sheet:
   - **M. abscessus**: bundled at `${CLAUDE_PLUGIN_ROOT}/refs/mabs/categories.xlsx`
     (sheet `2 Gene Annotations`, id column `Gene ID`, category column `Module`).
   - **M. tuberculosis / custom**: ask the user for an `.xlsx` and its columns.
     If no category sheet exists, say so plainly and stop — there is no bundled one.
3. Run: `${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq enrich <results.tsv> --categories <xlsx> --out <stem>`
   (optionally `--sheet --id-col --cat-col --padj --log2fc --alternative --min-size`).
4. Gene IDs are matched **case-insensitively**; the run **aborts if >5%** of DE IDs
   fail to match (a sign of wrong ID case/spelling — e.g. upper-casing `*c` genes).
5. Report the top over-represented sets + the bubble plot. NOTE: default is
   **two-sided** Fisher; the plot shows the over-represented side. FDRs are not
   comparable between one- and two-sided runs.
