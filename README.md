# bac-rnaseq

Bacterial RNAseq pipeline plugin for Claude Code and Codex. Feed raw Illumina
reads and a reference strain; get standard RNAseq outputs (counts, DESeq2
results) and publication-style figures. Validated for *Mycobacterium abscessus*
and *M. tuberculosis*; runs any bacterial genome (FASTA + GFF3). Lab-internal.

## Install (Claude Code)

```
/plugin marketplace add nicolai-tornow/bac-rnaseq
/plugin install bac-rnaseq@bac-rnaseq
```

## Install (Codex)

Copy or symlink `skills/` into `~/.codex/skills/` (see docs).

## Skills

- `setup-environment` — build the env + reference bundles, confirm the CPU budget.
- (later phases: `run-rnaseq`, `qc-triage`, `visualize-results`,
  `pathway-enrichment`, `batch-integration`, `export-results`)

## Correctness

Pipeline parameters come from the lab's validated pipelines. Load-bearing
invariants (strandedness, plasmid exclusion, feature counts, seqid
reconciliation) are enforced at runtime, not just in tests.

## Verification status

Built and unit/integration-tested on a dev machine (micromamba `deseq` env).
**Before publishing to labmates, the pipeline must be run once end-to-end on real
raw reads (planned on boulder).**

**Verified**
- Reference-bundle building on the real *M. abscessus* + *M. tuberculosis* genomes
  (Bowtie2 index, SAF, Mtb `H37RvBD`→`NC_018143.2` seqid remap; Mabs = 4970 features).
- DESeq2 recipe on the real `counts.tsv` (prefilter ≥10, `lfcShrink` normal, α=0.05;
  produces real differential-expression signal).
- Stage command parameters via unit tests (reverse-stranded `-s 2`, `--sensitive`,
  fastp Q20/len36, paired `--countReadPairs`).
- Orchestration wiring (stage order; `-s 2` reaches featureCounts).

**Still to verify — before publishing (boulder end-to-end run)**
- [ ] Full raw-FASTQ → trim → align → counts on real reads, including the
      **strand-correctness gate** (`tests/integration/test_correctness_gate.py`),
      which self-skips wherever no raw FASTQ are present.
- [ ] Creation of the pinned environment from `env/environment.yml` on the target
      host (resolves without version-pin conflicts).
- [ ] MultiQC aggregation in the `qc-triage` skill on a real run.
- [ ] Codex install path (`skills/` discoverable under `~/.codex/skills/`).
