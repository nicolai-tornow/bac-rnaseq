# Reference bundles

Committed inputs (Rock-lab annotation):

- `mabs/NC_010397.1.{fasta,gff3}`: *Mycobacterium abscessus* ATCC 19977 chromosome.
  The plasmid `NC_010394.1` is excluded.
- `mtb/H37RvBD.{fasta,gff3}`: *M. tuberculosis* H37Rv. The GFF seqid is `H37RvBD`
  but the FASTA header is `NC_018143.2`; `build-refs` remaps the GFF seqid to the
  FASTA one so the SAF `Chr` matches the alignment references.
- `mabs/categories.xlsx`: functional categories for `pathway-enrichment`.
- `<species>/structural_rna.tsv`: structural RNAs (see below).

Built by `build-refs` or `run`, outside this folder: `labels.saf`, `genome.fasta`
(+ `.fai`), `index/ref.*`, `bundle.json`. A bundle whose inputs are unchanged is
reused.

## Structural RNAs (`structural_rna.tsv`)

Columns: `GeneID Chr Start End Strand class in_gff source`.

These RNAs are abundant in total-RNA libraries, even after rRNA depletion. They are
used for two things:

1. **ncRNA % in QC:** the share of assigned reads that fall on any row.
2. **Counting:** rows with `in_gff = no` are added to the SAF, so their reads are
   assigned to them. Without them those reads are reported as `Unassigned_NoFeatures`
   (which looks like a strandedness problem), or are counted to an overlapping gene.
   These extra rows are written to `ncrna_counts.tsv` and kept out of the DESeq2
   matrix, which stays the GFF gene set.

`in_gff = yes` rows are GFF features: rRNA and tRNA by `gene_biotype` (mabs), or
by description (mtb, whose GFF has no biotypes).

`in_gff = no` rows were called with Infernal 1.1 `cmscan --cut_ga --nohmmonly`
against Rfam families RF00010, RF00011, RF00023, RF00169, RF00013 and RF02566,
and tmRNA was confirmed with Aragorn `-m`:

| Species | GeneID | Class | Coordinates | Evidence |
|---|---|---|---|---|
| mabs | MABnc_rnpB | RNase P RNA | 1,909,668-1,910,072 + | RF00010, score 298.4. Overlaps the first 79 bp of MAB1913 (same strand) |
| mabs | MABnc_ms1 | Ms1 RNA | 431,917-432,215 - | RF02566, score 230.4 |
| mabs | MABnc_ssrA | tmRNA | 3,513,404-3,513,772 - | RF00023, score 142.4; Aragorn 3,513,406-3,513,772 |
| mabs | MABnc_ffs | 4.5S SRP RNA | 306,563-306,657 + | RF00169, score 65.3. Overlaps the last 30 bp of MAB0304 (same strand) |
| mtb | RVBDnc_ms1 | Ms1 RNA (MTS2823) | 4,100,862-4,101,163 + | RF02566, score 227.8 |

For mtb, the same scan reproduces the GFF's `rnpB`, `ssr` and `4.5S` to the base,
so those are listed as `in_gff = yes`.

**Correction to the bundled mtb GFF (2026-09-18):** `ssr` (tmRNA, 3,467,964-3,468,331)
was annotated on the `+` strand. Rfam (RF00023, score 133.6) and Aragorn (`-m`,
c[3467965,3468331]) both place it on `-`, and with reverse-stranded counting the
`+` feature missed the tmRNA reads. The strand was changed to `-`; this is the only
edit to the original file.

A read that overlaps both a structural RNA and a gene is counted as ambiguous
(featureCounts default), so it no longer inflates MAB1913 or MAB0304.

**Custom genomes:** pass `reference.structural_rna: <table>` in the config. Without
a table, GFF genes whose `gene_biotype` is rRNA, tRNA, tmRNA, ncRNA, RNase_P_RNA or
SRP_RNA are used.
