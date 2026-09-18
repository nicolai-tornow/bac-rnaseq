---
name: setup-environment
description: Use FIRST, before running any bac-rnaseq analysis. Creates or locates the analysis environment, confirms a CPU budget, builds and verifies reference bundles (M. abscessus, M. tuberculosis, or a custom bacterial genome), and persists site config.
---

# Setup Environment

Follow these steps. Tool names per `skills/_shared/references/<harness>-tools.md`.
The CLI is `${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq`; run it with the env's Python.

## 1. Locate or create the environment

- If `~/.config/bac-rnaseq/site.yaml` has `env_prefix`, or `$RNASEQ_ENV_PREFIX`
  is set, or Lmod modules are configured: USE that shared environment (no install).
- Else, on a shared host (e.g. boulder), offer to create a shared env at a
  group-readable prefix the user names; otherwise create a personal env:
  `micromamba create -y -n bac-rnaseq -f "${CLAUDE_PLUGIN_ROOT}/env/environment.yml"`
  (shared: `micromamba create -y -p <prefix> -f .../env/environment.yml`).

## 2. Confirm the CPU budget

- Run `bac-rnaseq doctor`: it detects cores, suggests a budget (leaves 1-2 free),
  shows the saved budget and checks the tools. It does not change anything.
- ASK the user to confirm or override the suggestion, then save the confirmed value:
  `bac-rnaseq doctor --save-threads <N>`. Runs use it unless a config sets
  `resources.threads`.

## 3. Reference bundles

- Pick where bundles live and save it: `bac-rnaseq doctor --save-refs-root <dir>`
  (on a shared host, the shared folder the lab uses). Without it, each work
  directory builds its own bundle.
- For `mabs`/`mtb`: `bac-rnaseq build-refs --species <sp> --out <dir>/<sp> --threads <N>`
  (source files come from the plugin's `refs/`). A later run reuses the bundle as
  long as its inputs are unchanged.
- For a custom bacterial genome, the bundle is built on the first run from the
  config's `reference: {species: custom, fasta: ..., gff: ...}`.
- VERIFY: the reported feature count is > 0 and every SAF `Chr` is a FASTA
  sequence ID (the builder raises if not). For `mabs`, expect 4970 features.

## 4. Report

- Summarize: env location, confirmed threads, bundle folder, and per species the
  feature count, seqids and index path. The user is now ready for `run-rnaseq`.
