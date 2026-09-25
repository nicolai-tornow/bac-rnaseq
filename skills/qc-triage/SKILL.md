---
name: qc-triage
description: Use to interpret RNAseq QC and flag bad libraries before trusting differential expression. Reads the per-sample verdicts of a run (alignment rate, three-way strandedness check, assigned fraction, ncRNA fraction) and the MultiQC report, and recommends a fix for each problem.
---

# QC Triage

`run-rnaseq` computes QC automatically. Use this skill to explain the result, either
right after a run or later on a finished run directory.

1. Read `out/<run>/00_run_report.json` → `samples_qc`. For each sample it holds
   `verdict`, `reasons`, `alignment_pct`, `assigned_frac`, `nofeature_frac`,
   `ncrna_frac`, `ncrna_by_class` and `strandedness` (assigned fraction at reverse `-s 2`,
   forward `-s 1`, unstranded `-s 0`, plus the inferred setting). A sample below 95%
   alignment also has `unaligned`: the `reads` that did not align, their `gc_pct` and
   `rrna_like_frac` (share aligning locally to the reference's own rRNA genes), with
   the reads in `out/<run>/qc/unaligned/<sample>/`.
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
     Use `unaligned`: a high `rrna_like_frac` is rRNA of another organism
     (contamination, or a contaminant's rRNA that depletion missed); a `gc_pct` far from
     the reference genome's (*M. abscessus* 64.1%, *M. tuberculosis* 65.6%) also points
     to contamination; unaligned reads at the reference's GC% with little rRNA point to
     a strain or annotation difference. The user decides whether to re-run with
     `--allow-qc-fail --reason "<their reason>"`; ask them for the reason, never write
     one yourself.
   - **High ncRNA:** little usable mRNA depth. Use `ncrna_by_class` to tell why:
     rRNA+tRNA measures depletion efficiency, while Ms1, tmRNA and RNase P RNA are
     abundant by biology (in CF sputum they can exceed 75% of reads).
   - Always compare `ncrna_by_class` across replicates: one replicate with a much
     higher rRNA+tRNA share (e.g. 20% vs < 1%) was depleted less well, even when no
     rule fires.
4. A re-run reuses each sample's BAM while its completion marker matches (same FASTQs,
   layout, reference and trimming/alignment settings), so fixing a counting setting
   such as strandedness only repeats counting and DESeq2.
