# bac-rnaseq

Bacterial RNAseq pipeline. Skills live in `skills/`; real logic in `engine/`.
When a skill body says "run the CLI", it means `python -m engine.python.cli`
(inside the plugin's micromamba env `bac-rnaseq`, or a shared/Lmod env).

Harness tool mappings: `skills/_shared/references/<harness>-tools.md`.

Correctness invariants (strandedness, plasmid exclusion, feature counts, seqid
reconciliation) are enforced by the engine — never bypass them.
