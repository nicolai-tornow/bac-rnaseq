import pytest
from engine.python.gff_to_saf import gff_to_saf

GFF = """##gff-version 3
chrA\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;locus_tag=MAB0001
chrA\tsrc\tgene\t200\t300\t.\t-\t.\tID=g2;locus_tag=MAB0002c
plasmid\tsrc\tgene\t1\t50\t.\t+\t.\tID=gp;locus_tag=MABp01
chrA\tsrc\tCDS\t1\t100\t.\t+\t0\tID=c1;locus_tag=MAB0001
"""


def _write(tmp_path):
    p = tmp_path / "x.gff3"
    p.write_text(GFF)
    return str(p)


def test_basic_gene_rows(tmp_path):
    rows = gff_to_saf(_write(tmp_path), {"chrA", "plasmid"}, ["gene"], "locus_tag", [], {})
    ids = {r[0] for r in rows}
    assert ids == {"MAB0001", "MAB0002c", "MABp01"}
    assert ("MAB0002c", "chrA", 200, 300, "-") in rows


def test_exclude_seqid(tmp_path):
    rows = gff_to_saf(_write(tmp_path), {"chrA", "plasmid"}, ["gene"], "locus_tag", ["plasmid"], {})
    assert {r[0] for r in rows} == {"MAB0001", "MAB0002c"}


def test_seqid_map_applied_and_validated(tmp_path):
    rows = gff_to_saf(_write(tmp_path), {"NC_1", "plasmid"}, ["gene"], "locus_tag", ["plasmid"], {"chrA": "NC_1"})
    assert all(r[1] == "NC_1" for r in rows)


def test_chr_not_in_fasta_raises(tmp_path):
    with pytest.raises(ValueError):
        gff_to_saf(_write(tmp_path), {"other"}, ["gene"], "locus_tag", [], {})
