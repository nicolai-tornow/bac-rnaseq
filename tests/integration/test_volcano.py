import glob
import os
import pytest
from engine.viz.volcano import load_results, volcano
from tests.testdata import lab

CANDS = sorted(glob.glob(str(lab("rnaseq_mabs_media/06_deseq/results")) + "/*.tsv"))
pytestmark = pytest.mark.skipif(not CANDS, reason="no real result TSVs")


def test_volcano_renders_real(tmp_path):
    df = load_results(CANDS[0])
    assert {"log2FoldChange", "padj"}.issubset(df.columns)
    pdf, png = volcano(df, str(tmp_path / "v"),
                       selection={"mode": "top", "n": 10, "by": "l2fc"})
    assert os.path.exists(pdf) and os.path.exists(png)
