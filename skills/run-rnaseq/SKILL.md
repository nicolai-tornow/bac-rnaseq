---
name: run-rnaseq
description: Use to run the RNAseq pipeline from raw Illumina FASTQ to DESeq2 results. Needs a sample sheet, a reference selection (mabs/mtb/custom), and one or more contrasts. Produces counts, DESeq2 result tables, and a run report; auto-checks QC.
---

# Run RNAseq

Prereq: `setup-environment` has run (env + reference bundle ready). Tool names
per `skills/_shared/references/<harness>-tools.md`.

1. Collect inputs: a **sample sheet** (`sample_id, fastq_r1[, fastq_r2], condition[, replicate, batch]`),
   the **reference** (`mabs`/`mtb`/`custom` + FASTA/GFF for custom), and the
   **contrasts** (explicit `numerator vs denominator`, and/or `all_vs_all`).
   Confirm strandedness (default `reverse`) and the thread budget.
2. Write a config YAML (spec §7.2) and validate it:
   `python -m engine.python.cli validate <config.yaml>`.
3. Run: `python -m engine.python.cli run <config.yaml> --work-dir <dir> --samplesheet <tsv> --refs-root <refs>`.
4. Read `out/<run_name>/00_run_report.json`. If any sample's QC verdict is FAIL,
   STOP and surface it (do not present the DE as trustworthy) — invoke `qc-triage`
   to explain. Otherwise report the outputs (counts, per-contrast result tables).
