import pytest
from engine.python.config import load_config, Config


def test_mabs_defaults():
    cfg = load_config({"run_name": "r1", "reference": {"species": "mabs"}})
    assert isinstance(cfg, Config)
    assert cfg.reference.species == "mabs"
    assert cfg.reference.strandedness == "reverse"
    assert cfg.reference.feature_types == ["gene"]
    assert cfg.reference.id_attribute == "locus_tag"
    assert cfg.reference.exclude_seqids == ["NC_010394.1"]
    assert cfg.resources.threads is None


def test_custom_requires_fasta_and_gff():
    with pytest.raises(ValueError):
        load_config({"run_name": "r", "reference": {"species": "custom"}})


def test_custom_ok_with_paths():
    cfg = load_config({"run_name": "r", "reference": {
        "species": "custom", "fasta": "g.fa", "gff": "g.gff3"}})
    assert cfg.reference.exclude_seqids == []


def test_bad_strandedness_rejected():
    with pytest.raises(ValueError):
        load_config({"run_name": "r", "reference": {"species": "mabs", "strandedness": "maybe"}})


def test_unknown_key_rejected():
    # A misspelt key must not silently fall back to the default strandedness.
    with pytest.raises(ValueError):
        load_config({"run_name": "r", "reference": {"species": "mabs", "strandness": "forward"}})


def test_unknown_top_level_key_rejected():
    with pytest.raises(ValueError):
        load_config({"run_name": "r", "reference": {"species": "mabs"}, "treshold": {}})


def test_reads_layout_default_and_mate1_only():
    assert load_config({"run_name": "r", "reference": {"species": "mabs"}}).reads.layout == "auto"
    cfg = load_config({"run_name": "r", "reference": {"species": "mabs"},
                       "reads": {"layout": "mate1_only"}})
    assert cfg.reads.layout == "mate1_only"
