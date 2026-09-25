# Shared-server runbook (boulder)

How to set up `bac-rnaseq` once on a shared server so everyone in the lab uses one
environment and one set of reference bundles. Individual users then only need the
"Run an analysis" steps in the README.

`$SHARED` below is a group-readable folder (for example `/data/lab/bac-rnaseq`) and
`<lab-group>` the lab's Unix group. Ask the sysadmin for both.

## 1. Get the plugin

```bash
git clone https://github.com/nicolai-tornow/bac-rnaseq.git "$SHARED/bac-rnaseq"
```

The repository is public, so no GitHub account is needed. To update later:
`git -C "$SHARED/bac-rnaseq" pull`.

## 2. Shared environment

```bash
micromamba create -y -p "$SHARED/envs/bac-rnaseq" -f "$SHARED/bac-rnaseq/env/environment.yml"
chgrp -R <lab-group> "$SHARED/envs" && chmod -R g+rX "$SHARED/envs"
```

Each user then activates it with `micromamba activate "$SHARED/envs/bac-rnaseq"`
and sets `export PYTHONDONTWRITEBYTECODE=1` (keeps the shared install unchanged).
On boulder the pinned environment resolves (DESeq2 1.42.0, R 4.3.3).

## 3. Shared reference bundles

Bundles (bowtie2 index + SAF) are built once into a shared folder and reused by every
run whose reference inputs are unchanged.

```bash
B="$SHARED/bac-rnaseq/bin/bac-rnaseq"
mkdir -p "$SHARED/refs"
for sp in mabs mtb; do $B build-refs --species $sp --out "$SHARED/refs/$sp" --threads 8; done
chgrp -R <lab-group> "$SHARED/refs" && chmod -R g+rX "$SHARED/refs"
```

Expect `mabs: 4970 features on ['NC_010397.1']` and `mtb: ... on ['NC_018143.2']`.

Each user records the shared locations once:

```bash
$B doctor --save-refs-root "$SHARED/refs"
$B doctor                      # checks tools, suggests a thread budget
$B doctor --save-threads 8     # the budget you agreed with the sysadmin
```

After a plugin update that changes a reference input (FASTA, GFF or
`structural_rna.tsv`), rebuild the shared bundles as above. A user's run stops with
a clear message if the shared bundle is out of date and read-only for them.

## 4. Check the install

```bash
cd "$SHARED/bac-rnaseq" && python -m pytest tests -q
```

Expect every test to pass except those that need lab data: 79 passed, 5 skipped.
The strand-correctness gate and a full end-to-end run use the bundled read fixture
(`tests/fixtures/`), so they run here. To include the lab-data tests, set
`BAC_RNASEQ_TESTDATA` to a folder that contains `rnaseq_mabs_media/`.

## 5. Troubleshooting

- **A sample FAILs on strandedness.** The report shows the assigned fraction at all
  three `-s` settings and which one the library supports. Set
  `reference.strandedness` to that value and re-run. BAMs are reused.
- **Assigned fraction is low but strandedness is confirmed (WARN).** Reads fall
  outside the annotation: check `ncrna_counts.tsv`, the `Unassigned_NoFeatures`
  share and the alignment rate. Do not change strandedness.
- **Alignment below 90% (FAIL).** Common for a strain other than the reference
  (for example *M. abscessus* subsp. *massiliense* against ATCC 19977) or for host
  contamination. Decide whether to proceed; if so, re-run with `--allow-qc-fail`
  (alignment is reused).
- **`build-refs` raises "SAF Chr not in FASTA sequence IDs".** A custom genome's GFF
  seqids differ from its FASTA headers; add a `seqid_map` to `reference:`.
- **ComBat over-corrects.** There was no real batch effect; trust the "before" PCA.
- **Enrichment aborts ">5% unmatched".** DE gene IDs do not match the category
  sheet; check case and `--id-col`.
