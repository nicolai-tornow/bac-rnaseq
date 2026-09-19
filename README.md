# bac-rnaseq

A plugin for **Claude Code** and **Codex** that runs a bacterial RNAseq analysis from
raw Illumina reads to DESeq2 results, with checks that stop a mis-stranded or
contaminated library from silently producing a differential-expression table.
Validated for *Mycobacterium abscessus* and *M. tuberculosis*; runs any bacterial
genome with a FASTA and a GFF3.

Developed by [Nicolai Tornow](https://github.com/nicolaitornow). Source-available,
not open source: see [License](#license).

## What it does

| Skill | What it does |
|---|---|
| **setup-environment** | Creates or finds the analysis environment, confirms a CPU budget, builds reference bundles (*M. abscessus*, *M. tuberculosis*, or a custom genome). |
| **run-rnaseq** | FASTQ → FastQC → fastp → bowtie2 → featureCounts → QC gate → DESeq2, with a machine-readable run report. |
| **qc-triage** | Explains the per-sample QC verdicts (alignment, strandedness, assigned fraction, ncRNA fraction) and recommends fixes. |
| **visualize-results** | Volcano, TPM bar/heatmap and expression-ranking figures, selected by pathway, locus tags, or top/bottom N genes. |
| **pathway-enrichment** | Fisher over-representation of functional categories (BH-FDR) with a bubble plot. |
| **batch-integration** | ComBat/PCA to look at batch effects across runs (differential expression stays per batch). |
| **export-results** | One Excel workbook: summary, per-contrast tables with gene names, normalized/VST/TPM values. |
| **report-feedback** | Files a bug report or suggestion as a GitHub issue on this repository. |

You talk to the agent ("run RNAseq on these FASTQs, SCFM2 vs 7H9"); the skills tell
it which commands to run and what to check.

## Install

No GitHub account or credentials are needed. You need `git`, `micromamba` (or
conda) and Claude Code or Codex. Tested on Linux and WSL.

**Claude Code**

```bash
claude plugin marketplace add https://github.com/nicolai-tornow/bac-rnaseq.git
claude plugin install bac-rnaseq@bac-rnaseq
```

Start a new session afterwards. Alternatively, without the marketplace:

```bash
git clone https://github.com/nicolai-tornow/bac-rnaseq.git
claude --plugin-dir ./bac-rnaseq
```

**Codex** (not yet tested)

```bash
codex plugin marketplace add https://github.com/nicolai-tornow/bac-rnaseq.git
codex plugin add bac-rnaseq@bac-rnaseq
```

**The analysis environment** is created by the `setup-environment` skill, or by hand:

```bash
micromamba create -y -n bac-rnaseq -f env/environment.yml
```

On a shared server, set up one environment and one set of reference bundles for
everyone: see [docs/BOULDER_RUNBOOK.md](docs/BOULDER_RUNBOOK.md).

## Run an analysis

Ask the agent to set up and run, or use the command line directly. The command is
`bin/bac-rnaseq` in the plugin folder; it works from any directory inside the
`bac-rnaseq` environment. To type just `bac-rnaseq`:

```bash
micromamba activate bac-rnaseq
export PATH="/path/to/bac-rnaseq/bin:$PATH"
```

**1. Sample sheet** (tab-separated; `fastq_r2`, `replicate` and `batch` are optional):

```
sample_id	fastq_r1	fastq_r2	condition
7H9_rep1	/data/7H9_1_R1.fastq.gz	/data/7H9_1_R2.fastq.gz	7H9
SCFM2_rep1	/data/SCFM2_1_R1.fastq.gz	/data/SCFM2_1_R2.fastq.gz	SCFM2
...
```

**2. Config** (`config.yaml`). Unknown keys are rejected, so a typo cannot silently
fall back to a default.

```yaml
run_name: media
reference:
  species: mabs            # mabs | mtb | custom (custom also needs fasta: and gff:)
  strandedness: reverse    # default; checked on the data during the run
contrasts:
  explicit:
    - {name: SCFM2_vs_7H9, numerator: SCFM2, denominator: 7H9}
# optional:
# reads: {layout: mate1_only}      # mixed single/paired-end sheet: run all from read 1
# design: {batch_variable: batch}  # ~ batch + condition
# resources: {threads: 8}
# thresholds: {padj: 0.05, log2fc: 1}
```

**3. Validate, then run:**

```bash
bac-rnaseq validate config.yaml --samplesheet samples.tsv
bac-rnaseq run config.yaml --work-dir . --samplesheet samples.tsv
```

`validate` catches missing FASTQs, unsafe sample IDs, mixed read layouts and contrast
levels that are not in the sample sheet, before any work starts.

**4. Outputs** in `out/<run_name>/`:

| Path | Contents |
|---|---|
| `00_run_report.json` | Status, per-sample QC, DE gene counts per contrast, tool versions, plugin commit, reference checksums |
| `00_inputs/` | The config and sample sheet exactly as run |
| `05_counts/counts.tsv` | Raw counts, GFF genes × samples (the DESeq2 input) |
| `05_counts/ncrna_counts.tsv` | Counts for structural RNAs the GFF lacks (see below) |
| `05_counts/strand_check/` | featureCounts at the two other strand settings |
| `06_deseq/results/<contrast>.tsv` | DESeq2 results (`lfcShrink` normal, alpha 0.05) |
| `06_deseq/` | Normalized counts, VST, size factors |
| `qc/multiqc/` | MultiQC over FastQC, fastp, bowtie2 and featureCounts |

**5. Figures and tables:** `visualize-results`, `pathway-enrichment` and
`export-results` work on any DESeq2 result table, not only on runs made here.

## Quality control

Every run checks each sample and **stops before DESeq2 if any sample FAILs**
(`status: "qc_fail"`, exit code 1). Counts and QC are still written.

| Check | FAIL | WARN |
|---|---|---|
| Alignment rate | < 90% | < 95% |
| Strandedness | the declared setting is contradicted by the data | declared unstranded on a stranded library |
| Assigned fraction | — | < 60% while strandedness is confirmed |
| Structural ncRNA share of assigned reads | — | > 85% (reported per class: rRNA, tRNA, tmRNA, RNase P, Ms1, SRP) |

**Strandedness is measured, not assumed.** Reads are counted at all three
featureCounts settings (`-s 2` reverse, `-s 1` forward, `-s 0` unstranded). A
stranded library assigns at least 5x more reads on its own strand. Unstranded
counting is usually a little *higher* than the correct stranded setting, because
it also counts antisense reads; that is expected. A low assigned fraction with
confirmed strandedness means reads fall outside the annotation, which is a WARN, not
a reason to change the strand setting.

**After a FAIL:** fix the cause and run again, or accept it deliberately with
`--allow-qc-fail`. Samples whose BAM is complete are not re-trimmed or re-aligned,
so a re-run only repeats counting and DESeq2.

## Reference bundles and structural RNAs

The bundled genomes use the Rock-lab annotation: *M. abscessus* ATCC 19977
(`NC_010397.1`, plasmid excluded, 4,970 features) and *M. tuberculosis* H37Rv.

Total-RNA libraries contain abundant structural RNAs even after rRNA depletion.
The *M. abscessus* GFF lacks four of them: Ms1 RNA, tmRNA, RNase P RNA and 4.5S SRP
RNA. In one public dataset they took most of the reads, which then looked like a
strandedness problem. The bundles add these RNAs (called with Rfam/Infernal and
Aragorn; coordinates and evidence in [refs/README.md](refs/README.md)), so their reads
are counted and reported in `ncrna_counts.tsv`. They are **kept out of the DESeq2
matrix**, which stays the GFF gene set, so new count tables line up row for row with
tables from earlier runs of this pipeline. Reads overlapping both a structural RNA
and a gene are counted for neither.

For a custom genome, pass `reference.structural_rna: <table>` in the same format;
otherwise GFF genes with an rRNA, tRNA, tmRNA or ncRNA `gene_biotype` are used for
the ncRNA share.

Bundles are built once and reused while their inputs are unchanged. Save a shared
location with `bac-rnaseq doctor --save-refs-root <dir>`.

## Things to know

- **Gene IDs keep the GFF's case** (`MAB0001`, `MAB3648`). Some earlier tables from
  the Rock lab use lowercase tags (`mab0001`). Gene matching in enrichment is
  case-insensitive; when you join tables yourself, normalize the case first.
- **Strains other than the reference.** Reads from *M. abscessus* subsp.
  *massiliense* or *bolletii* can be run against the ATCC 19977 bundle, but genes
  the strain has and ATCC 19977 lacks are invisible, and the alignment rate is lower
  (it may FAIL the 90% gate for that reason alone). Read the alignment rate with
  that in mind, and treat genes with zero counts in every sample as possibly absent.
- **Mixed read layouts** are rejected by default. `reads: {layout: mate1_only}` runs
  every sample single-end from read 1, so layout is not confounded with condition.
- **Threads:** `resources.threads` in the config, else the budget saved with
  `bac-rnaseq doctor --save-threads N`, else 4.

## Tests

```bash
python -m pytest tests -q
```

The suite includes a full run on real reads: 20,000 read pairs of a public
reverse-stranded *M. abscessus* library bundled in `tests/fixtures/`. Tests that need
unpublished lab data skip unless `BAC_RNASEQ_TESTDATA` points at it.

## Verification status

- [x] Pinned environment from `env/environment.yml` resolves on boulder
      (DESeq2 1.42.0, R 4.3.3).
- [x] Raw FASTQ → counts on real reads, including the strand-correctness gate, run
      on the bundled public fixture (Linux/WSL).
- [x] MultiQC over a real run (bundled fixture).
- [x] A full-size real dataset on boulder (0.1.0 plus manual QC workarounds): the
      public *M. abscessus* subsp. *massiliense* 1239 RNAseq (PRJNA602697; 9 libraries,
      7H9 / SCFM2 / CF sputum). It reproduces the published DE tables (Pearson r
      0.935-0.960, 92-98% same direction, all five qRT-PCR genes agree), and a
      mate-1-only run agrees with the paired-end run at r = 0.989. The QC problems it
      exposed are what 0.2.0 fixes.
- [ ] A full-size run on boulder with 0.2.0.
- [ ] Codex install path.

## Feedback

Open an issue on this repository, or ask the agent to use the `report-feedback`
skill.

## License

Copyright (c) 2026 Nicolai Tornow. Source-available, **not** open source. You may
download, install, run and modify this software for your own research or internal
use, including private copies and GitHub forks. Redistributing it (or a modified
version) outside GitHub's fork feature, or using it commercially, needs written
permission. The full terms are in [LICENSE](LICENSE).
