# bac-rnaseq

Run a bacterial RNAseq analysis end to end — raw Illumina reads to DESeq2 results
and figures — in **Codex and Claude Code**. Validated for *Mycobacterium abscessus*
and *M. tuberculosis*; runs any bacterial genome (FASTA + GFF3).

Developed by [Nicolai Tornow](https://github.com/nicolaitornow). I'm sharing this
plugin privately with colleagues. **Ask me for repository access and send me your
GitHub username.**

## Included skills

| Skill | What it does |
|---|---|
| **setup-environment** | Builds or locates the analysis environment, confirms a CPU budget, and builds reference bundles (*M. abscessus*, *M. tuberculosis*, or a custom genome). |
| **run-rnaseq** | Raw FASTQ → QC → trim → align → count → DESeq2, with validated parameters and a machine-readable run report. |
| **qc-triage** | Per-sample QC gate: alignment rate, assigned-read fraction (a wrong-strandedness detector) and ncRNA fraction. |
| **visualize-results** | Volcano, TPM bar/heatmap and expression-ranking figures, selected by pathway/operon, locus tags, or top-N/bottom-N genes. |
| **pathway-enrichment** | Fisher over-representation of functional categories (BH-FDR) with a bubble plot. |
| **batch-integration** | ComBat/PCA to visualize batch effects across runs (differential expression stays per-batch). |
| **export-results** | One tidy Excel workbook: summary, per-contrast tables with gene names, plus normalized/VST/TPM. |

## Installation

You need access to this private repository, Git, GitHub CLI, `micromamba` (or conda),
and Codex or Claude Code. Linux/WSL is tested. Sign in with your own GitHub account
using `gh auth login`, then run `gh auth setup-git`.

### Codex

```bash
codex plugin marketplace add https://github.com/nicolai-tornow/plugins.git
codex plugin add bac-rnaseq@nicolai-tornow
```

Browse with `/plugins`. Start a new session after installation.

### Claude Code

```bash
claude plugin marketplace add https://github.com/nicolai-tornow/plugins.git
claude plugin install bac-rnaseq@nicolai-tornow
```

Browse with `/plugin`. Start a new session after installation.

## Run an analysis

1. **Set up once** — `setup-environment` builds the environment and the reference
   bundle, and asks you to confirm a CPU budget.
2. **Run** — give `run-rnaseq` a sample sheet
   (`sample_id, fastq_r1[, fastq_r2], condition[, replicate, batch]`), a reference
   (`mabs` / `mtb` / `custom`), and one or more contrasts. It produces counts,
   per-contrast DESeq2 tables, and `00_run_report.json`.
3. **Figures and tables** — `visualize-results`, `pathway-enrichment` and
   `export-results` run on any DESeq2 result table; a full run is not required.

Keep data and outputs in their own workspace, outside the plugin installation.

## Quality control has teeth

`run-rnaseq` runs `qc-triage` automatically. **If any sample FAILs QC** (alignment
below 90%, or an assigned-read fraction that implies the wrong strandedness), the run
**stops before DESeq2**: it writes the counts and the QC report (`status: "qc_fail"`)
and exits non-zero, so a contaminated or mis-stranded library can never silently
produce a differential-expression table. Fix the cause (correct `strandedness`, or
drop / re-sequence the library) and re-run, or pass `--allow-qc-fail` to proceed
deliberately. WARN samples proceed but are flagged.

## Correctness

Pipeline parameters come from the lab's validated pipelines; load-bearing invariants
(reverse-stranded counting, plasmid exclusion, feature counts, seqid reconciliation,
case-insensitive gene IDs) are enforced at runtime, not just in tests.

## Verification status

Built and unit/integration-tested (Linux/WSL, micromamba). **Before broad use, run
once end to end on real raw reads on boulder** — see `docs/BOULDER_RUNBOOK.md`.
Still to verify there:

- [ ] Full raw-FASTQ → counts on real reads, including the strand-correctness gate
      (`tests/integration/test_correctness_gate.py`, which self-skips without raw FASTQ).
- [ ] Pinned environment creation from `env/environment.yml` resolves on the target host.
- [ ] MultiQC aggregation in `qc-triage` on a real run.
- [ ] Codex install path.

Please ask me before redistributing this private package.
