import pandas as pd
import pytest
from engine.python.batch import read_common_vst, assert_same_reference


def test_common_genes(tmp_path):
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    pd.DataFrame({"s1": [1, 2, 3]}, index=["g1", "g2", "g3"]).to_csv(a, sep="\t")
    pd.DataFrame({"s2": [4, 5]}, index=["g2", "g3"]).to_csv(b, sep="\t")
    mat = read_common_vst([str(a), str(b)])
    assert list(mat.index) == ["g2", "g3"]
    assert list(mat.columns) == ["s1", "s2"]


def test_diff_reference_raises(tmp_path):
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    pd.DataFrame({"s1": [1, 2, 3]}, index=["g1", "g2", "g3"]).to_csv(a, sep="\t")
    pd.DataFrame({"s2": [4, 5]}, index=["x1", "x2"]).to_csv(b, sep="\t")
    with pytest.raises(ValueError):
        assert_same_reference([str(a), str(b)])


def test_common_genes_case_insensitive(tmp_path):
    # One run written with upper-case locus tags, another with lower-case.
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    pd.DataFrame({"s1": [1, 2, 3]}, index=["MAB0001", "MAB0002c", "MAB0003"]).to_csv(a, sep="\t")
    pd.DataFrame({"s2": [4, 5]}, index=["mab0002c", "mab0003"]).to_csv(b, sep="\t")
    assert_same_reference([str(a), str(a)])
    mat = read_common_vst([str(a), str(b)])
    assert list(mat.index) == ["MAB0002c", "MAB0003"]      # first input's spelling
    assert mat.loc["MAB0002c", "s2"] == 4
