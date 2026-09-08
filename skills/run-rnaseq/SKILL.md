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
4. Read `out/<run_name>/00_run_report.json`. **A FAIL halts the run before DESeq2**
   (`status: "qc_fail"`, CLI exit 1) so bad data cannot silently produce a DE table —
   only counts + QC are written. Surface the failing sample(s) and reasons (invoke
   `qc-triage`); the user should fix the cause (e.g. correct `strandedness`, or drop /
   re-sequence the bad library) and re-run, or pass `--allow-qc-fail` to proceed
   deliberately. WARN samples proceed but are flagged. On `status: "ok"`, report the
   outputs (counts, per-contrast result tables).
