---
name: qc-triage
description: Use to interpret RNAseq QC and flag bad libraries before trusting differential expression. Reads the per-sample verdicts of a run (alignment rate, three-way strandedness check, assigned fraction, ncRNA fraction) and the MultiQC report, and recommends a fix for each problem.
---

# QC Triage

`run-rnaseq` computes QC automatically. Use this skill to explain the result, either
right after a run or later on a finished run directory.

1. Read `out/<run>/00_run_report.json` → `samples_qc`. For each sample it holds
   `verdict`, `reasons`, `alignment_pct`, `assigned_frac`, `nofeature_frac`,
   `ncrna_frac` and `strandedness` (assigned fraction at reverse `-s 2`,
   forward `-s 1`, unstranded `-s 0`, plus the inferred setting).
   The MultiQC report is at `out/<run>/qc/multiqc/multiqc_report.html`.
2. Rules:

   | Check | FAIL | WARN |
   |---|---|---|
   | Alignment rate | < 90% | < 95% |
   | Strandedness | declared setting contradicted (the other strand assigns ≥ 5x more, or the library looks unstranded) | declared unstranded on a stranded library, or not inferable |
   | Assigned fraction | — | < 60% (or < 40%) with strandedness confirmed |
   | ncRNA (rRNA, tRNA, tmRNA, RNase P, Ms1, SRP) | — | > 85% of assigned reads |

3. Report a PASS/WARN/FAIL table with the reasons, then recommend:
   - **Strandedness FAIL:** re-run with the `strandedness` the check inferred.
   - **Low assigned, strandedness confirmed:** do NOT change strandedness. The reads
     fall outside the annotation: unannotated RNA, contamination, or an incomplete
     GFF (check `ncrna_counts.tsv` and the NoFeatures share).
   - **Low alignment:** contamination, or a strain that differs from the reference.
     The user decides whether to re-run with `--allow-qc-fail`.
   - **High ncRNA:** poor rRNA depletion; the sample has little usable mRNA depth.
4. A re-run reuses trimmed reads and BAMs, so fixing a setting only repeats counting
   and DESeq2.
