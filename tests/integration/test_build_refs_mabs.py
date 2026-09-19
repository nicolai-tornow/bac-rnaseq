import shutil
import pytest
from pathlib import Path
from engine.python.build_refs import build_bundle, fasta_seqids

REPO = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(not shutil.which("bowtie2-build"),
                                reason="bowtie2-build not on PATH")


def test_build_mabs_bundle(tmp_path):
    res = build_bundle("mabs", refs_root=REPO / "refs", out_dir=tmp_path, threads=2)
    assert res["n_features"] == 4970
    assert res["seqids"] == ["NC_010397.1"]
    assert Path(res["saf"]).exists()
    assert Path(str(res["index_prefix"]) + ".1.bt2").exists()
    lines = Path(res["saf"]).read_text().splitlines()[1:]
    assert all(l.split("\t")[1] == "NC_010397.1" for l in lines)


def test_fasta_seqids_mabs():
    ids = fasta_seqids(REPO / "refs/mabs/NC_010397.1.fasta")
    assert ids == {"NC_010397.1"}


def test_build_mtb_bundle_remaps_seqid(tmp_path):
    res = build_bundle("mtb", refs_root=REPO / "refs", out_dir=tmp_path, threads=2)
    # GFF seqid 'H37RvBD' must be remapped to the FASTA header 'NC_018143.2'
    assert res["seqids"] == ["NC_018143.2"]
    assert res["n_features"] > 4000
    lines = Path(res["saf"]).read_text().splitlines()[1:]
    assert all(line.split("\t")[1] == "NC_018143.2" for line in lines)


def test_mabs_bundle_counts_structural_rnas(tmp_path):
    res = build_bundle("mabs", refs_root=REPO / "refs", out_dir=tmp_path, threads=2)
    # The GFF lacks Ms1, tmRNA, RNase P RNA and 4.5S SRP RNA; the bundle adds them
    # to the SAF but they are not part of the 4970-feature DESeq2 gene set.
    assert set(res["extra_ids"]) == {"MABnc_rnpB", "MABnc_ms1", "MABnc_ssrA", "MABnc_ffs"}
    assert not set(res["extra_ids"]) & set(res["gene_ids"])
    saf_ids = [l.split("\t")[0] for l in Path(res["saf"]).read_text().splitlines()[1:]]
    assert len(saf_ids) == 4970 + 4
    classes = res["structural"]
    assert classes["MABr5051"] == "rRNA" and classes["MABnc_ms1"] == "Ms1_RNA"
    assert sum(c == "tRNA" for c in classes.values()) == 47


def test_bundle_is_reused_and_leaves_plugin_refs_clean(tmp_path):
    first = build_bundle("mabs", refs_root=REPO / "refs", out_dir=tmp_path, threads=2)
    second = build_bundle("mabs", refs_root=REPO / "refs", out_dir=tmp_path, threads=2)
    assert first["reused"] is False and second["reused"] is True
    assert (tmp_path / "genome.fasta.fai").exists()
    assert not (REPO / "refs/mabs/NC_010397.1.fasta.fai").exists()


def test_custom_bundle_honours_exclude_and_feature_types(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr\n" + "A" * 400 + "\n>plasmid\n" + "C" * 200 + "\n")
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "chr\tx\tCDS\t1\t90\t.\t+\t.\tID=c1;locus_tag=G1\n"
        "chr\tx\tgene\t1\t90\t.\t+\t.\tID=g1;locus_tag=G1\n"
        "chr\tx\tCDS\t120\t300\t.\t-\t.\tID=c2;locus_tag=G2\n"
        "plasmid\tx\tCDS\t1\t90\t.\t+\t.\tID=c3;locus_tag=P1\n")
    res = build_bundle("custom", out_dir=tmp_path / "b", threads=1, fasta=fa, gff=gff,
                       feature_types=("CDS",), exclude_seqids=["plasmid"])
    assert res["gene_ids"] == ["G1", "G2"] and res["seqids"] == ["chr"]


def test_mtb_tmrna_on_minus_strand(tmp_path):
    # Rfam RF00023 and Aragorn both place H37Rv tmRNA (ssr) on the minus strand.
    res = build_bundle("mtb", refs_root=REPO / "refs", out_dir=tmp_path, threads=2)
    row = [l.split("\t") for l in Path(res["saf"]).read_text().splitlines()
           if l.startswith("ssr\t")]
    assert row and row[0][4] == "-"
    assert res["structural"]["RVBDnc_ms1"] == "Ms1_RNA"
