from __future__ import annotations
import hashlib
import json
import os
import subprocess
from pathlib import Path
from .gff_to_saf import gff_to_saf, write_saf

PLUGIN_REFS = Path(__file__).resolve().parents[2] / "refs"

BUNDLES = {
    "mabs": {"fasta": "mabs/NC_010397.1.fasta", "gff": "mabs/NC_010397.1.gff3",
             "structural_rna": "mabs/structural_rna.tsv",
             "exclude_seqids": ["NC_010394.1"], "seqid_map": {}},
    "mtb":  {"fasta": "mtb/H37RvBD.fasta", "gff": "mtb/H37RvBD.gff3",
             "structural_rna": "mtb/structural_rna.tsv",
             "exclude_seqids": [], "seqid_map": {"H37RvBD": "NC_018143.2"}},
}
STRUCTURAL_BIOTYPES = {"rRNA", "tRNA", "tmRNA", "RNase_P_RNA", "SRP_RNA", "ncRNA"}


def md5(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fasta_seqids(fasta) -> set[str]:
    ids = set()
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def read_structural_rna(path) -> list[dict]:
    """Structural RNA table: GeneID Chr Start End Strand class in_gff source."""
    lines = Path(path).read_text().splitlines()
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l.strip()]


def _structural_from_biotype(gff, id_attribute):
    """Custom genome without a table: GFF genes whose gene_biotype is structural."""
    rows = []
    for line in Path(gff).read_text().splitlines():
        f = line.split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "gene":
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if a.get("gene_biotype") in STRUCTURAL_BIOTYPES and id_attribute in a:
            rows.append({"GeneID": a[id_attribute], "class": a["gene_biotype"], "in_gff": "yes"})
    return rows


def build_bundle(species, refs_root=None, out_dir=None, threads=4,
                 fasta=None, gff=None, feature_types=("gene",),
                 id_attribute="locus_tag", exclude_seqids=None, seqid_map=None,
                 structural_rna=None):
    """Build (or reuse) a reference bundle: bowtie2 index + SAF + structural-RNA IDs.

    The SAF holds the GFF features (the DESeq2 gene set) plus any structural RNAs
    the GFF lacks, so featureCounts sees both: reads from an unannotated ncRNA are
    then assigned to it, not left as NoFeatures or given to an overlapping gene.
    A bundle whose inputs are unchanged is reused (bundle.json stamp).
    """
    refs_root = Path(refs_root) if refs_root else PLUGIN_REFS
    out_dir = Path(out_dir)
    if species in BUNDLES:
        b = BUNDLES[species]
        fasta = refs_root / b["fasta"]
        gff = refs_root / b["gff"]
        structural_rna = structural_rna or refs_root / b["structural_rna"]
        exclude_seqids = b["exclude_seqids"] if exclude_seqids is None else exclude_seqids
        seqid_map = b["seqid_map"] if seqid_map is None else seqid_map
    else:  # custom
        fasta, gff = Path(fasta), Path(gff)
        exclude_seqids = exclude_seqids or []
        seqid_map = seqid_map or {}

    stamp = {"species": species, "fasta_md5": md5(fasta), "gff_md5": md5(gff),
             "structural_rna_md5": md5(structural_rna) if structural_rna else None,
             "feature_types": list(feature_types), "id_attribute": id_attribute,
             "exclude_seqids": list(exclude_seqids), "seqid_map": dict(seqid_map)}
    meta_path = out_dir / "bundle.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if meta.get("stamp") == stamp and Path(meta["index_prefix"] + ".1.bt2").exists():
            meta["reused"] = True
            return meta

    seqids = fasta_seqids(fasta)
    gene_rows = gff_to_saf(gff, seqids, list(feature_types), id_attribute,
                           list(exclude_seqids), dict(seqid_map))
    gene_ids = {r[0] for r in gene_rows}
    structural = (read_structural_rna(structural_rna) if structural_rna
                  else _structural_from_biotype(gff, id_attribute))
    extra_rows = []
    for s in structural:
        if s.get("in_gff") == "no":
            if s["Chr"] not in seqids:
                raise ValueError(f"structural RNA {s['GeneID']}: Chr '{s['Chr']}' "
                                 f"not in FASTA sequence IDs {sorted(seqids)}")
            if s["GeneID"] in gene_ids:
                raise ValueError(f"structural RNA ID {s['GeneID']} duplicates a GFF feature")
            extra_rows.append((s["GeneID"], s["Chr"], int(s["Start"]), int(s["End"]), s["Strand"]))
        elif s["GeneID"] not in gene_ids:
            raise ValueError(f"structural RNA {s['GeneID']} is marked in_gff=yes but "
                             "is not a feature in the SAF")

    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        saf = write_saf(gene_rows + extra_rows, out_dir / "labels.saf")
    except OSError as e:
        raise RuntimeError(f"reference bundle at {out_dir} is missing or out of date and "
                           f"cannot be rebuilt there ({e}); ask the owner of the shared "
                           "bundle to rebuild it, or build a personal one") from e

    # Index a link to the FASTA inside the bundle, so the .fai never lands in the
    # (possibly read-only) plugin install.
    link = out_dir / "genome.fasta"
    if link.is_symlink() or link.exists():
        link.unlink()
    os.symlink(Path(fasta).resolve(), link)
    subprocess.run(["samtools", "faidx", str(link)], check=True)
    idx_dir = out_dir / "index"
    idx_dir.mkdir(exist_ok=True)
    prefix = idx_dir / "ref"
    subprocess.run(["bowtie2-build", "--quiet", "--threads", str(threads),
                    str(fasta), str(prefix)], check=True)

    meta = {"stamp": stamp, "saf": str(saf), "saf_md5": md5(saf),
            "index_prefix": str(prefix), "fasta": str(fasta), "gff": str(gff),
            "n_features": len(gene_rows), "seqids": sorted({r[1] for r in gene_rows}),
            "gene_ids": sorted(gene_ids),
            "extra_ids": [r[0] for r in extra_rows],
            "structural": {s["GeneID"]: s["class"] for s in structural},
            "reused": False}
    meta_path.write_text(json.dumps(meta, indent=1))
    return meta
