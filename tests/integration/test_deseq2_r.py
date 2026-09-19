import csv
import re
import shutil
import subprocess
import pytest
from pathlib import Path
from tests.testdata import lab

REPO = Path(__file__).resolve().parents[2]
REAL_COUNTS = lab("rnaseq_mabs_media/05_counts/counts.tsv")
pytestmark = pytest.mark.skipif(
    not (shutil.which("Rscript") and REAL_COUNTS.exists()),
    reason="need Rscript + real Mabs counts.tsv")


def _condition(sample):
    return re.sub(r"_rep\d+$", "", sample)


def test_deseq2_on_real_mabs_counts(tmp_path):
    # Use REAL pipeline output — synthetic matrices are too uniform for DESeq2's
    # parametric dispersion fit (the lab's real-data-over-fixtures lesson).
    header = REAL_COUNTS.read_text().splitlines()[0].split("\t")
    samples = header[1:]
    coldata = tmp_path / "coldata.tsv"
    coldata.write_text("sample\tcondition\n" +
                       "".join(f"{s}\t{_condition(s)}\n" for s in samples))
    (tmp_path / "contrasts.tsv").write_text("SCFM2_vs_7H9\tSCFM2\t7H9\n")
    out = tmp_path / "out"
    r = subprocess.run(["Rscript", str(REPO / "engine/r/deseq2.R"),
                        str(REAL_COUNTS), str(coldata),
                        str(tmp_path / "contrasts.tsv"), str(out), "0"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    res = out / "results" / "SCFM2_vs_7H9.tsv"
    assert res.exists()
    lines = res.read_text().splitlines()
    assert lines[0].split("\t") == ["Gene", "baseMean", "log2FoldChange",
                                    "lfcSE", "stat", "pvalue", "padj"]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    sig = [x for x in rows if x["padj"] not in ("NA", "") and float(x["padj"]) < 0.05]
    assert len(rows) > 4000          # prefilter keeps most of the 4970 genes
    assert len(sig) > 50             # real SCFM2-vs-7H9 DE signal
    assert (out / "normalized_counts.tsv").exists()
    assert (out / "vst_counts.tsv").exists()
