"""cleanup: what a dry run lists, and every reason it refuses."""
import json
import os
import time
from pathlib import Path
import pytest
from engine.python import cleanup as CL
from engine.python.runlock import RunLock


def _age(root, secs=3600):
    t = time.time() - secs
    for p in [root, *Path(root).rglob("*")]:
        os.utime(p, (t, t), follow_symlinks=False)


def make_run(tmp_path, status="ok", sheet=None):
    """A finished run folder with every kind of file a real one holds."""
    raw = tmp_path / "raw"
    raw.mkdir()
    r1, r2 = raw / "A_R1.fq.gz", raw / "A_R2.fq.gz"
    r1.write_bytes(b"raw1")
    r2.write_bytes(b"raw2")
    out = tmp_path / "work" / "out" / "r"
    files = {
        "00_run_report.json": json.dumps({"run_name": "r", "status": status}),
        "00_inputs/samplesheet.tsv": sheet or
        f"sample_id\tfastq_r1\tfastq_r2\tcondition\ns1\t{r1}\t{r2}\tA\n",
        "00_inputs/config.yaml": "run_name: r\n",
        "01_qc_raw/A_R1_fastqc.html": "x",
        "02_trimmed/s1_R1.fq.gz": "t" * 1000, "02_trimmed/s1_R2.fq.gz": "t" * 1000,
        "02_trimmed/s1.json": "{}", "02_trimmed/s1.html": "h",
        "04_align/s1.bam": "b" * 3000, "04_align/s1.bam.bai": "i" * 10,
        "04_align/s1.bowtie2.log": "log", "04_align/s1.done.json": "{}",
        "04_align/s0.sam": "s" * 5000,
        "05_counts/counts.tsv": "c", "05_counts/temp-core-000001-ab.sam": "x" * 700,
        "05_counts/strand_check/temp-core-000002-cd.sam": "x" * 300,
        "05_counts/temp-sort-000001-ab-": "x" * 10,
        "06_deseq/results/B_vs_A.tsv": "d", "qc/multiqc/multiqc_report.html": "q"}
    for rel, text in files.items():
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        (out / rel).write_text(text)
    _age(tmp_path)
    return out, (r1, r2)


def _files(root):
    return sorted(str(p) for p in Path(root).rglob("*") if p.is_file() or p.is_symlink())


def test_dry_run_lists_tiers_and_deletes_nothing(tmp_path):
    out, _ = make_run(tmp_path)
    before = _files(tmp_path)
    p = CL.plan(out)
    assert p.refusals == []
    t = p.by_tier()
    assert t["trimmed"] == (2, 2000) and t["sam"] == (1, 5000)
    assert t["fc_temp"] == (3, 1010) and t["bam"] == (2, 3010)
    assert {c.tier for c in p.selected} == {"trimmed", "sam", "fc_temp"}
    assert p.kept_files == len(_files(out)) - 8   # all but the 8 tier files (bam has its own row)
    text = CL.format_plan(p)
    assert "delete" in text and "--include-bams" in text
    assert _files(tmp_path) == before


def test_include_bams_selects_bams(tmp_path):
    out, _ = make_run(tmp_path)
    p = CL.plan(out, include_bams=True)
    assert {c.tier for c in p.selected} == {"trimmed", "sam", "fc_temp", "bam"}
    assert any("re-aligns" in n for n in p.notes)


def test_run_without_markers_gets_a_note(tmp_path):
    out, _ = make_run(tmp_path)
    (out / "04_align/s1.done.json").unlink()
    _age(tmp_path)
    assert any("0.3.0" in n for n in CL.plan(out).notes)


def _refusal(tmp_path, case):
    out, (r1, r2) = make_run(tmp_path, status="qc_fail" if case == "qc_fail" else "ok")
    trimmed = out / "02_trimmed/s1_R1.fq.gz"
    cwd = None
    if case == "no_report":
        (out / "00_run_report.json").unlink()
    elif case == "partial":
        (out / "04_align/.s1.partial").mkdir()
    elif case == "no_sheet":
        (out / "00_inputs/samplesheet.tsv").unlink()
    elif case == "symlink":
        (tmp_path / "elsewhere").write_text("x")
        trimmed.unlink()
        trimmed.symlink_to(tmp_path / "elsewhere")
    elif case == "outside":
        ext = tmp_path / "scratch"
        (out / "02_trimmed").rename(ext)
        (out / "02_trimmed").symlink_to(ext)
    elif case == "raw_by_path":
        (out / "00_inputs/samplesheet.tsv").write_text(
            f"sample_id\tfastq_r1\tcondition\ns1\t{trimmed}\tA\n")
    elif case == "raw_by_hardlink":
        trimmed.unlink()
        os.link(r1, trimmed)
    elif case == "relative_path":       # relative to the work dir; cleanup runs elsewhere
        (out / "00_inputs/samplesheet.tsv").write_text(
            "sample_id\tfastq_r1\tcondition\ns1\tout/r/02_trimmed/s1_R1.fq.gz\tA\n")
        cwd = tmp_path
    elif case == "csv_sheet":           # run copies a CSV sheet verbatim to samplesheet.tsv
        (out / "00_inputs/samplesheet.tsv").write_text(
            f"sample_id,fastq_r1,condition\ns1,{trimmed},A\n")
    _age(tmp_path)
    if case == "recent":
        (out / "05_counts/counts.tsv").write_text("new")
    lk = RunLock(out) if case == "live_lock" else None
    if lk:
        lk.acquire()
    try:
        snap = _files(tmp_path)
        p = CL.plan(out, cwd=cwd)
        assert _files(tmp_path) == snap               # a plan never deletes
        return p
    finally:
        if lk:
            lk.release()


@pytest.mark.parametrize("case,words", [
    ("no_report", "00_run_report.json"), ("qc_fail", "qc_fail"), ("live_lock", "using this folder"),
    ("partial", "staging"), ("recent", "modified"), ("no_sheet", "samplesheet.tsv"),
    ("symlink", "symlink"), ("outside", "outside the run folder"),
    ("raw_by_path", "raw FASTQ"), ("raw_by_hardlink", "raw FASTQ"),
    ("relative_path", "raw FASTQ"), ("csv_sheet", "raw FASTQ")])
def test_refusals(tmp_path, case, words):
    p = _refusal(tmp_path, case)
    assert any(words in r for r in p.refusals), p.refusals
    assert "REFUSED" in CL.format_plan(p)


def test_idle_minutes_zero_skips_the_recent_check(tmp_path):
    out, _ = make_run(tmp_path)
    (out / "05_counts/counts.tsv").write_text("new")
    assert CL.plan(out, idle_minutes=0).refusals == []


def _manifest(out):
    return (out / "cleanup_manifest.tsv").read_text().splitlines()


def test_execute_deletes_exactly_what_the_plan_lists(tmp_path):
    out, _ = make_run(tmp_path)
    p = CL.plan(out)
    expected = {str(c.path) for c in p.selected}
    before = set(_files(tmp_path))
    entry = CL.execute(out)
    assert before - set(_files(tmp_path)) == expected
    assert entry["files"] == len(expected) == 6
    assert entry["bytes"] == sum(c.bytes for c in p.selected)
    rows = _manifest(out)
    assert rows[0] == "path\tbytes\ttier\ttimestamp\tplugin_commit" and len(rows) == 7
    assert {r.split("\t")[0] for r in rows[1:]} == {str(Path(e).relative_to(out)) for e in expected}
    rep = json.loads((out / "00_run_report.json").read_text())
    assert rep["status"] == "ok" and rep["cleanup"][0]["files"] == 6
    assert not (out / ".lock").exists()


def test_second_cleanup_appends_to_manifest_and_report(tmp_path):
    out, _ = make_run(tmp_path)
    CL.execute(out)
    _age(tmp_path)
    entry = CL.execute(out, include_bams=True)
    assert entry["files"] == 2 and "bam" in entry["tiers"]
    rows = _manifest(out)
    assert rows.count(rows[0]) == 1 and len(rows) == 1 + 6 + 2
    assert len(json.loads((out / "00_run_report.json").read_text())["cleanup"]) == 2


@pytest.mark.parametrize("case", ["qc_fail", "recent", "raw_by_path", "live_lock"])
def test_execute_refuses_and_deletes_nothing(tmp_path, case):
    out, _ = make_run(tmp_path, status="qc_fail" if case == "qc_fail" else "ok")
    if case == "raw_by_path":
        (out / "00_inputs/samplesheet.tsv").write_text(
            f"sample_id\tfastq_r1\tcondition\ns1\t{out / '02_trimmed/s1_R1.fq.gz'}\tA\n")
        _age(tmp_path)
    if case == "recent":
        (out / "05_counts/counts.tsv").write_text("new")
    lk = RunLock(out) if case == "live_lock" else None
    if lk:
        lk.acquire()
    try:
        before = _files(tmp_path)
        with pytest.raises(CL.CleanupRefused):
            CL.execute(out)
        assert _files(tmp_path) == before
    finally:
        if lk:
            lk.release()
