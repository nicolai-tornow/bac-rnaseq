from __future__ import annotations
import subprocess
from pathlib import Path
from .gff_to_saf import gff_to_saf, write_saf

BUNDLES = {
    "mabs": {"fasta": "mabs/NC_010397.1.fasta", "gff": "mabs/NC_010397.1.gff3",
             "exclude_seqids": ["NC_010394.1"], "seqid_map": {}},
    "mtb":  {"fasta": "mtb/H37RvBD.fasta", "gff": "mtb/H37RvBD.gff3",
             "exclude_seqids": [], "seqid_map": {"H37RvBD": "NC_018143.2"}},
}


def fasta_seqids(fasta) -> set[str]:
    ids = set()
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            ids.add(line[1:].split()[0])
    return ids


def build_bundle(species, refs_root, out_dir, threads=4,
                 fasta=None, gff=None, feature_types=("gene",),
                 id_attribute="locus_tag", exclude_seqids=None, seqid_map=None):
    refs_root, out_dir = Path(refs_root), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if species in BUNDLES:
        b = BUNDLES[species]
        fasta = refs_root / b["fasta"]
        gff = refs_root / b["gff"]
        exclude_seqids = b["exclude_seqids"] if exclude_seqids is None else exclude_seqids
        seqid_map = b["seqid_map"] if seqid_map is None else seqid_map
    else:  # custom
        fasta, gff = Path(fasta), Path(gff)
        exclude_seqids = exclude_seqids or []
        seqid_map = seqid_map or {}

    seqids = fasta_seqids(fasta)
    rows = gff_to_saf(gff, seqids, list(feature_types), id_attribute,
                      list(exclude_seqids), dict(seqid_map))
    saf = write_saf(rows, out_dir / "labels.saf")

    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    idx_dir = out_dir / "index"
    idx_dir.mkdir(exist_ok=True)
    prefix = idx_dir / "ref"
    subprocess.run(["bowtie2-build", "--quiet", "--threads", str(threads),
                    str(fasta), str(prefix)], check=True)

    used = sorted({r[1] for r in rows})
    return {"saf": str(saf), "index_prefix": str(prefix),
            "n_features": len(rows), "seqids": used, "fasta": str(fasta)}
