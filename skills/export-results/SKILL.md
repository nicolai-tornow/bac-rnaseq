---
name: export-results
description: Use to package a run's DESeq2 results into one tidy multi-sheet Excel workbook — a summary sheet, one sheet per contrast with human-readable gene names, plus normalized/VST/TPM matrices.
---

# Export Results

Consumes a run's `06_deseq/results/` directory (and optionally normalized/VST/TPM).
Tool names per `skills/_shared/references/<harness>-tools.md`.

1. Gene names are derived from the reference GFF (`product`/`Name`/`gene`), so pass
   the same GFF used for the run (e.g. `${CLAUDE_PLUGIN_ROOT}/refs/mabs/NC_010397.1.gff3`).
2. Run: `python -m engine.python.cli export --results-dir out/<run>/06_deseq/results --gff <gff> --out results.xlsx [--normalized <tsv> --vst <tsv>]`.
3. The workbook has: `summary` (sig-gene counts per contrast), one sheet per contrast
   (DE table + a `gene_name` column), and `normalized`/`vst` sheets when provided.
   Report the workbook path.
