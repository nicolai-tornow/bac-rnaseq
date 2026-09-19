"""Full run on real reads: trim -> align -> count -> strand check -> QC -> MultiQC.

DESeq2 is stubbed: four 5k-pair slices of one library are too small for DESeq2's
vst (it needs >= 1000 expressed genes). DESeq2 itself is covered by
test_deseq2_r.py on real count tables.
"""
import gzip
import json
import shutil
import pytest
from engine.python.config import load_config
from engine.python import run as R
from tests.testdata import READS_R1, READS_R2

TOOLS = ["fastp", "fastqc", "multiqc", "bowtie2", "samtools", "featureCounts"]
pytestmark = pytest.mark.skipif(not all(shutil.which(t) for t in TOOLS),
                                reason="need the full bac-rnaseq environment")


def _split_fixture(dst, n=4, pairs=5000):
    rows = ["sample_id\tfastq_r1\tfastq_r2\tcondition"]
    for m, src in ((1, READS_R1), (2, READS_R2)):
        lines = gzip.open(src, "rt").read().splitlines(True)
        for i in range(n):
            with gzip.open(dst / f"s{i}_R{m}.fq.gz", "wt") as fh:
                fh.writelines(lines[i * pairs * 4:(i + 1) * pairs * 4])
    for i in range(n):
        rows.append(f"s{i}\t{dst}/s{i}_R1.fq.gz\t{dst}/s{i}_R2.fq.gz\t{'A' if i < n // 2 else 'B'}")
    (dst / "ss.tsv").write_text("\n".join(rows) + "\n")
    return dst / "ss.tsv"


def test_full_run_on_real_reads(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    deseq = []
    monkeypatch.setattr(R, "_run_deseq2", lambda *a, **k: deseq.append(a))
    ss = _split_fixture(tmp_path)
    cfg = load_config({"run_name": "e2e", "reference": {"species": "mabs"},
                       "resources": {"threads": 2},
                       "contrasts": {"explicit": [{"name": "B_vs_A", "numerator": "B",
                                                   "denominator": "A"}]}})
    rep = R.run_pipeline(cfg, tmp_path, None, ss)
    out = tmp_path / "out" / "e2e"

    assert rep["status"] == "ok" and deseq
    for sid, v in rep["samples_qc"].items():
        assert v["strandedness"]["inferred"] == "reverse", sid
        assert v["verdict"] == "PASS", (sid, v["reasons"])
        by = v["ncrna_by_class"]                       # depletion vs biology, kept apart
        assert by["Ms1_RNA"] > 0.05 and by["tmRNA"] > 0.05 and "rRNA" in by
        assert rep["provenance"]["reads"][sid]["reads_in"] == 10000   # fastp read everything

    counts = (out / "05_counts/counts.tsv").read_text().splitlines()
    assert len(counts) - 1 == 4970                     # DESeq2 matrix = GFF gene set
    nc = {l.split("\t")[0]: l.split("\t")[1:] for l in
          (out / "05_counts/ncrna_counts.tsv").read_text().splitlines()[1:]}
    assert set(nc) == {"MABnc_rnpB", "MABnc_ms1", "MABnc_ssrA", "MABnc_ffs"}
    assert int(nc["MABnc_ms1"][0]) > 100               # Ms1 is abundant in this library

    assert not list((out / "04_align").glob("*.sam"))  # intermediate SAMs removed
    assert (out / "qc/multiqc/multiqc_report.html").exists()
    assert list((out / "01_qc_raw").glob("*_fastqc.html"))
    assert list((out / "03_qc_trimmed").glob("*_fastqc.html"))
    rpt = json.loads((out / "00_run_report.json").read_text())
    assert rpt["provenance"]["reference"]["saf_md5"]
    assert (out / "00_inputs/config.yaml").exists()

    # Re-run: every sample resumes from its BAM and the bundle is reused.
    rep2 = R.run_pipeline(cfg, tmp_path, None, ss)
    assert all(r["resumed"] for r in rep2["provenance"]["reads"].values())
    assert rep2["provenance"]["reference"]["reused"] is True
