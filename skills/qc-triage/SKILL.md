---
name: qc-triage
description: Use to interpret RNAseq QC and flag bad libraries before trusting differential expression. Aggregates FastQC/MultiQC and checks per-sample alignment rate, assigned-read fraction (a wrong-strandedness detector), and ncRNA fraction.
---

# QC Triage

Runs automatically at the end of `run-rnaseq`; also invocable standalone on a
finished run directory.

1. Aggregate reports: `multiqc -f -o out/<run>/qc out/<run>` (FastQC + fastp + bowtie2 + featureCounts).
2. Read the per-sample verdicts in `00_run_report.json` (`samples_qc`). Rules:
   alignment <90% = FAIL, <95% = WARN; assigned fraction <40% = FAIL
   (**likely wrong strandedness** — offer to re-run with the correct `strandedness`),
   <60% = WARN; ncRNA >85% = WARN (low usable mRNA depth — the batch_1 SCFM2 lesson).
3. Report a PASS/WARN/FAIL table with reasons. On any FAIL, recommend the fix
   (re-run with corrected strandedness, or drop/re-sequence the contaminated library)
   before using the DE results.
