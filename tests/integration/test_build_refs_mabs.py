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
