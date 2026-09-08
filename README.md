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
