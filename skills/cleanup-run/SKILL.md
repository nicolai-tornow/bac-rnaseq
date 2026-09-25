---
name: cleanup-run
description: Use when the user wants disk space back from a finished RNAseq run. Removes regenerable intermediates (trimmed FASTQs, leftover SAMs, featureCounts temp files, optionally BAMs) from out/<run_name>/ after showing a dry run and getting explicit confirmation. Never runs on its own.
---

# Clean Up a Run

Only when the user asks for it (or accepts an offer). Never as part of a run.
Tool names per `skills/_shared/references/<harness>-tools.md`. The CLI is
`${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq` (run it with the env's Python).

1. Dry run: `bac-rnaseq cleanup out/<run_name>`. It deletes nothing. Show its table
   to the user as printed: files and GB per tier, marked `delete` or `keep`.
   - `trimmed` (trimmed FASTQs), `sam` (SAMs left by engines before 0.3.0) and
     `fc_temp` (featureCounts temp files from crashed runs) are deleted by default.
   - `bam` is kept unless `--include-bams`. Logs, QC, counts, DESeq2 results, the run
     report and the completion markers are always kept. Raw FASTQs are never touched.
2. If it says REFUSED (exit code 2), explain the reason and stop. Never delete files by
   hand to get around a refusal. A refusal only because files changed recently (for
   example right after a run) says when to retry; use `--idle-minutes 0` only if the
   user explicitly asks for it.
3. Ask whether to delete the files marked `delete`. Ask separately about BAMs: they are
   needed for IGV and for re-counting, and without them the next run re-trims and
   re-aligns those samples. Mention any `note:` lines.
4. Only after an explicit yes: `bac-rnaseq cleanup out/<run_name> --yes`, adding
   `--include-bams` only if the user explicitly agreed to delete BAMs. Never pass
   `--yes` without that answer.
5. Report the freed GB and the manifest, `out/<run_name>/cleanup_manifest.tsv` (one row
   per deleted file). The run report's `cleanup` list records each clean-up.
