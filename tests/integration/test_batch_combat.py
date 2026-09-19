import shutil
import pytest
import pandas as pd
from pathlib import Path
from engine.python.batch import run_combat
from tests.testdata import lab

VST = lab("rnaseq_mabs_media/06_deseq/vst_counts.tsv")
pytestmark = pytest.mark.skipif(not (VST.exists() and shutil.which("Rscript")),
                                reason="need real VST + Rscript")


def test_combat_runs_on_real_vst(tmp_path):
    df = pd.read_csv(VST, sep="\t", index_col=0)
    cols = [c for c in df.columns if c.startswith("7H9") or c.startswith("SCFM2")]
    df = df[cols]
    meta = pd.DataFrame(index=cols)
    meta["condition"] = ["7H9" if c.startswith("7H9") else "SCFM2" for c in cols]
    # two batches, each containing both conditions (not confounded)
    meta["batch"] = ["A" if c.endswith(("rep1", "rep2")) else "B" for c in cols]
    res = run_combat(df, meta, tmp_path)
    assert Path(res["corrected"]).exists()
    assert Path(res["pca_before"]).exists() and Path(res["pca_after"]).exists()
    corrected = pd.read_csv(res["corrected"], sep="\t", index_col=0)
    assert corrected.shape[0] == df.shape[0]   # gene count preserved
