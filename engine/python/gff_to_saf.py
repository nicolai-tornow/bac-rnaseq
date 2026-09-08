from __future__ import annotations
import re
from pathlib import Path


def _attr(col9: str, key: str):
    m = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", col9)
    return m.group(1) if m else None


def gff_to_saf(gff_path, fasta_seqids, feature_types, id_attribute,
               exclude_seqids, seqid_map):
    feature_types = set(feature_types)
    exclude = set(exclude_seqids)
    rows = []
    for line in Path(gff_path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in feature_types:
            continue
        seqid = seqid_map.get(f[0], f[0])
        if seqid in exclude:
            continue
        gid = _attr(f[8], id_attribute)
        if gid is None:
            continue
        if seqid not in fasta_seqids:
            raise ValueError(
                f"SAF Chr '{seqid}' (from GFF seqid '{f[0]}') not in FASTA "
                f"sequence IDs {sorted(fasta_seqids)} — counting would drop all reads")
        rows.append((gid, seqid, int(f[3]), int(f[4]), f[6]))
    return rows


def write_saf(rows, out_path):
    with open(out_path, "w") as fh:
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return out_path
