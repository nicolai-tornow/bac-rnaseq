---
name: run-rnaseq
description: Use to run the RNAseq pipeline from raw Illumina FASTQ to DESeq2 results. Needs a sample sheet, a reference selection (mabs/mtb/custom), and one or more contrasts. Produces counts, DESeq2 result tables, and a run report; auto-checks QC.
---

# Run RNAseq

Prereq: `setup-environment` has run (env + reference bundle ready). Tool names
per `skills/_shared/references/<harness>-tools.md`. The CLI is
`${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq` (run it with the env's Python).

1. Collect inputs: a **sample sheet** (`sample_id, fastq_r1[, fastq_r2], condition[, replicate, batch]`),
   the **reference** (`mabs`/`mtb`/`custom` + FASTA/GFF for custom), and the
   **contrasts** (explicit `numerator vs denominator`, and/or `all_vs_all`).
   Confirm strandedness (default `reverse`; the run checks it) and the thread budget.
   Sample IDs may contain letters, digits, `_`, `.` and `-`.
   If the sheet mixes single- and paired-end samples, set `reads: {layout: mate1_only}`
   (every sample then runs single-end from `fastq_r1`).
2. Write a config YAML and validate it together with the sample sheet:
   `bac-rnaseq validate <config.yaml> --samplesheet <tsv>`. Unknown keys are errors.
3. Run: `bac-rnaseq run <config.yaml> --work-dir <dir> --samplesheet <tsv>`.
4. Read `out/<run_name>/00_run_report.json`. **A FAIL halts the run before DESeq2**
   (`status: "qc_fail"`, exit code 1); counts and QC are still written. Explain it
   with `qc-triage`. The user fixes the cause and re-runs, or passes
   `--allow-qc-fail` to proceed deliberately. A re-run skips samples whose completion
   marker still matches their inputs. WARN samples proceed but are flagged.
   `status: "failed"` (exit 1): `failure` in the report names the stage, sample and
   error; half-written files were removed, so fix the cause and re-run. Exit 2 with
   "in use by another run": another run holds the folder; do not work around it.
   On `status: "ok"`, report the outputs and `de_summary` (DE gene counts per contrast).

Outputs in `out/<run_name>/`: `05_counts/counts.tsv` (GFF genes, the DESeq2 input),
`05_counts/ncrna_counts.tsv` (structural RNAs added by the bundle),
`05_counts/strand_check/`, `06_deseq/results/<contrast>.tsv`, `qc/multiqc/`,
`00_inputs/` (config and sample sheet as run) and `00_run_report.json`.
