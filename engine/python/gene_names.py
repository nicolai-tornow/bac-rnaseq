from __future__ import annotations
import re
from pathlib import Path


def _attr(col9, key):
    m = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", col9)
    return m.group(1) if m else None


def gene_names_from_gff(gff, id_attr="locus_tag",
                        name_attrs=("gene", "Name", "gene_name", "product")):
    """Map locus tag -> a human-readable name, scanning every feature line
    (gene/CDS/…) and taking the first name attribute found."""
    out = {}
    for line in Path(gff).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        gid = _attr(f[8], id_attr)
        if not gid:
            continue
        for a in name_attrs:
            v = _attr(f[8], a)
            if v:
                out.setdefault(gid, v)
                break
    return out
