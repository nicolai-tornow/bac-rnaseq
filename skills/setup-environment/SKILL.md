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
- Samples are trimmed and aligned several at a time. The default is one sample per
  8 threads (a 32-thread budget runs 4 samples x 8 threads; 8 threads runs one).
- ASK the user to confirm or override both numbers, then save them:
  `bac-rnaseq doctor --save-threads <N> --save-parallel <P>`. Runs use them unless a
  config sets `resources.threads` / `resources.parallel_samples`. Each parallel
  sample's sort uses about 2 GB of memory at 8 threads.

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
