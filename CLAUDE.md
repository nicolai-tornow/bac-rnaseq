# bac-rnaseq

Bacterial RNAseq pipeline. Skills live in `skills/`; real logic in `engine/`.
When a skill body says "run the CLI", it means `bin/bac-rnaseq <command>`
(`${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq`), run with the Python of the plugin's
micromamba env `bac-rnaseq`, or a shared/Lmod env. It works from any directory.

Harness tool mappings: `skills/_shared/references/<harness>-tools.md`.

Correctness invariants (strandedness, plasmid exclusion, feature counts, seqid
reconciliation, structural RNAs kept out of the DESeq2 matrix) are enforced by the
engine — never bypass them.

Tests: `python -m pytest tests -q` from the plugin root. Lab data for the extra
real-data tests is found under `$BAC_RNASEQ_TESTDATA`.
