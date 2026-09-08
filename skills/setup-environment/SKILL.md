---
name: setup-environment
description: Use FIRST, before running any bac-rnaseq analysis. Creates or locates the analysis environment, confirms a CPU budget, builds and verifies reference bundles (M. abscessus, M. tuberculosis, or a custom bacterial genome), and persists site config.
---

# Setup Environment

Follow these steps. Tool names per `skills/_shared/references/<harness>-tools.md`.

## 1. Locate or create the environment

- If `~/.config/bac-rnaseq/site.yaml` has `env_prefix`, or `$RNASEQ_ENV_PREFIX`
  is set, or Lmod modules are configured: USE that shared environment (no install).
- Else, on a shared host (e.g. boulder), offer to create a shared env at a
  group-readable prefix the user names; otherwise create a personal env:
  `micromamba create -y -n bac-rnaseq -f "${CLAUDE_PLUGIN_ROOT}/env/environment.yml"`
  (shared: `micromamba create -y -p <prefix> -f .../env/environment.yml`).

## 2. Confirm the CPU budget

- Run `doctor` (`python -m engine.python.cli doctor`) to detect cores, get a
  suggested budget (leaves 1-2 free), and check tool presence.
- ASK the user to confirm or override the suggested thread count; save it to
  `site.yaml`.

## 3. Build + verify the reference bundle

- For `mabs`/`mtb`: `... build-refs --species <sp> --refs-root ${CLAUDE_PLUGIN_ROOT}/refs --out <refs_root>/<sp> --threads <T>`.
- For a custom bacterial genome: pass `--species custom --fasta <FASTA> --gff <GFF3>`.
- VERIFY: the reported feature count is > 0 and every SAF `Chr` is a FASTA
  sequence ID (the builder raises if not). For `mabs`, expect 4970 features.
- On a shared host, build once into the shared `refs_root`; other users reuse it.

## 4. Report

- Summarize: env location, confirmed threads, per-species bundle (feature count,
  seqids, index path). The user is now ready for `run-rnaseq`.
