# Plan 4 — Boulder Integration & Publish Runbook

This is the deployment runbook (not a code plan): get `bac-rnaseq` onto boulder,
stand up a **shared** environment and reference bundles with the sysadmin, run the
pipeline **end-to-end once** on real reads, run the deferred strand-correctness
gate, verify Codex, then publish to labmates. Tick the checklist in the top-level
`README.md` ("Still to verify") as you go.

> Status when this was written: Plans 1–3 built and locally tested (53 pass / 1
> skip). The one skip is the strand gate below — it needs raw FASTQ, which live on
> boulder. **This runbook closes it.**

---

## 0. Prerequisites

- SSH access to boulder; ability to start Claude Code or Codex there.
- A shared, group-readable filesystem path on boulder and a common Unix group for
  the lab (ask the sysadmin). Call the path `$SHARED` below (e.g. `/data/lab/bac-rnaseq`).
- `micromamba` available on boulder (or Lmod modules for the tools — see B2).
- One real Illumina dataset with raw FASTQ (e.g. the 7H9/SCFM2/ALI libraries).

---

## A. Get the plugin onto boulder

```bash
git clone https://github.com/nicolai-tornow/bac-rnaseq.git
cd bac-rnaseq
```
(Or, inside Claude Code on boulder: `/plugin marketplace add nicolai-tornow/bac-rnaseq`
then `/plugin install bac-rnaseq@bac-rnaseq` — but for the first end-to-end run a
plain clone is simplest.)

---

## B. Shared environment (with the sysadmin)

**B1 — Shared micromamba env (recommended).** Build once at a group path:
```bash
micromamba create -y -p "$SHARED/envs/rnaseq" -f env/environment.yml
chgrp -R <lab-group> "$SHARED/envs" && chmod -R g+rX "$SHARED/envs" && chmod g+s "$SHARED/envs"
export PYTHONDONTWRITEBYTECODE=1   # keep the read-only env immutable at runtime
```
Every user thereafter: `micromamba run -p "$SHARED/envs/rnaseq" <cmd>` — no reinstall.
Record it so `setup-environment` finds it:
```bash
mkdir -p ~/.config/bac-rnaseq
printf 'env_prefix: %s\n' "$SHARED/envs/rnaseq" >> ~/.config/bac-rnaseq/site.yaml
```

**B2 — Lmod alternative.** If boulder already exposes these tools as modules, load
them instead of building an env, and record `use_lmod: true` + the module list in
`site.yaml`. Confirm all six are present:
```bash
micromamba run -p "$SHARED/envs/rnaseq" bash -lc \
  'for t in fastp bowtie2 samtools featureCounts multiqc Rscript; do command -v $t || echo "MISSING $t"; done'
```

**⚠ Verify the env resolves** (a README checklist item): if `micromamba create`
errors on a version-pin conflict, that must be fixed in `env/environment.yml` and
re-pushed before labmates install.

---

## C. Shared reference bundles (build once)

```bash
E="micromamba run -p $SHARED/envs/rnaseq"
for sp in mabs mtb; do
  $E python -m engine.python.cli build-refs --species $sp \
     --refs-root ./refs --out "$SHARED/refs/$sp" --threads 8
done
chgrp -R <lab-group> "$SHARED/refs" && chmod -R g+rX "$SHARED/refs"
printf 'refs_root: %s\n' "$SHARED/refs" >> ~/.config/bac-rnaseq/site.yaml
```
Expect `mabs: 4970 features on ['NC_010397.1']` and `mtb: … on ['NC_018143.2']`.

---

## D. End-to-end run on real reads

Write a sample sheet `samples.tsv` (tab-separated) pointing at the raw FASTQ:
```
sample_id   fastq_r1                 fastq_r2                 condition   replicate   batch
7H9_rep1    /path/7H9_1_R1.fastq.gz  /path/7H9_1_R2.fastq.gz  7H9         1           Batch2
...
```
Write `config.yaml`:
```yaml
run_name: mabs_media
reference: {species: mabs}      # strandedness defaults to reverse (-s 2)
design: {variable: condition, batch_variable: batch}
contrasts:
  explicit:
    - {name: SCFM2_vs_7H9, numerator: SCFM2, denominator: 7H9}
    - {name: ALI_vs_7H9,   numerator: ALI,   denominator: 7H9}
resources: {threads: 8}
```
Run:
```bash
$E python -m engine.python.cli validate config.yaml
$E python -m engine.python.cli run config.yaml --work-dir . --samplesheet samples.tsv --refs-root "$SHARED/refs"
```
Inspect `out/mabs_media/00_run_report.json` — every sample's `samples_qc.verdict`
should be PASS (a FAIL on assigned-fraction means wrong strandedness).

**Cross-check against the known-good counts** (the real correctness proof):
```bash
$E python - <<'PY'
import pandas as pd
a = pd.read_csv("out/mabs_media/05_counts/counts.tsv", sep="\t", index_col=0)
b = pd.read_csv("/path/to/rnaseq_mabs_media/05_counts/counts.tsv", sep="\t", index_col=0)
common = a.index.intersection(b.index)
print("gene overlap:", len(common), "of", len(b))
# same libraries should match closely; identical libraries should match exactly
PY
```

---

## E. Run the deferred strand-correctness gate  ← the key verification

With raw reads present, the skipped test now runs:
```bash
$E bash -c 'PYTHONPATH=$PWD python -m pytest tests/integration/test_correctness_gate.py -v'
```
**Must PASS** ("reverse assigns the most reads"). If it skips, the reads glob in
the test doesn't match your layout — point it at a real `*_R1_001.fastq.gz` pair.
Then run the whole suite: `pytest tests -q` → expect all pass, 0 skipped.

---

## F. Downstream skills on the real run

```bash
$E python -m engine.python.cli visualize out/mabs_media/06_deseq/results/SCFM2_vs_7H9.tsv --out fig_volcano --top 10
$E python -m engine.python.cli enrich    out/mabs_media/06_deseq/results/SCFM2_vs_7H9.tsv --categories refs/mabs/categories.xlsx --out fig_enrich
$E python -m engine.python.cli export    --results-dir out/mabs_media/06_deseq/results --gff refs/mabs/NC_010397.1.gff3 --out results.xlsx
```

---

## G. Codex install verification

```bash
mkdir -p ~/.codex/skills
cp -r skills/* ~/.codex/skills/     # or symlink
```
Start Codex on boulder and confirm it discovers the seven skills (setup-environment,
run-rnaseq, qc-triage, visualize-results, pathway-enrichment, batch-integration,
export-results) and that a skill body's `python -m engine.python.cli …` command runs.

---

## H. Patch, re-push, publish

1. Tick every box in `README.md` → "Still to verify". Fix anything that broke on
   boulder (paths, env pins, strandedness edge cases) and commit.
2. Push the updated version:
   ```bash
   git commit -am "boulder: end-to-end validated; <patches>"
   git push origin main
   ```
   (Do **not** add a `Co-Authored-By: Claude` trailer.)
3. Announce to labmates: install with
   `/plugin marketplace add nicolai-tornow/bac-rnaseq` then
   `/plugin install bac-rnaseq@bac-rnaseq` (they must be in the `nicolai-tornow`
   org). Codex users: copy `skills/` into `~/.codex/skills/`.

---

## Troubleshooting

- **All samples FAIL on assigned fraction** → strandedness is wrong for this
  library prep; set `reference.strandedness: forward` or `unstranded` and re-run.
- **`build-refs` raises "SAF Chr not in FASTA sequence IDs"** → a custom genome's
  GFF seqids don't match its FASTA headers; pass a `seqid_map` (see `build_refs.py`).
- **ComBat over-corrects** → there was no real batch effect; trust the "before" PCA.
- **Enrichment aborts ">5% unmatched"** → DE gene IDs don't match the category
  sheet's spelling; check case and the `--id-col`.
