"""Raw FASTQ md5: from sample-sheet columns or an md5sum file; checked by validate and run."""
import hashlib
import pytest
from engine.python import checksums as K
from engine.python import run as R
from engine.python.procs import CallableRunner
from engine.python.samplesheet import Sample, read_samplesheet, check_samplesheet
from tests.unit.fake_tools import FakeTools, write_fastq_gz


def _md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def test_manifest_formats_and_relative_paths(tmp_path):
    (tmp_path / "d").mkdir()
    f = tmp_path / "d" / "a.fq.gz"
    f.write_bytes(b"x")
    m = tmp_path / "d" / "md5.txt"
    m.write_text(f"{'A' * 32}  a.fq.gz\n{'b' * 32} *sub/b.fq.gz\n# comment\n")
    man = K.read_md5_manifest(m)
    assert K.lookup(man, str(f)) == "a" * 32                  # by resolved path
    assert K.lookup(man, "/elsewhere/b.fq.gz") == "b" * 32    # by file name
    assert K.lookup(man, "/elsewhere/c.fq.gz") is None


def test_ambiguous_file_name_needs_a_full_path(tmp_path):
    m = tmp_path / "md5.txt"
    m.write_text(f"{'a' * 32}  L1/x.fq.gz\n{'b' * 32}  L2/x.fq.gz\n")
    man = K.read_md5_manifest(m)
    with pytest.raises(ValueError, match="x.fq.gz"):
        K.lookup(man, "/data/x.fq.gz")
    assert K.lookup(man, str(tmp_path / "L2/x.fq.gz")) == "b" * 32


def test_bad_manifest_line(tmp_path):
    (tmp_path / "m").write_text("not-an-md5  a.fq.gz\n")
    with pytest.raises(ValueError, match="md5sum"):
        K.read_md5_manifest(tmp_path / "m")


def test_expected_prefers_columns_and_skips_dropped_mates(tmp_path):
    man = {"a.fq.gz": "1" * 32, "b.fq.gz": "2" * 32}
    s1 = Sample("s1", "/d/a.fq.gz", "c", "/d/b.fq.gz", md5_r1="f" * 32)
    s2 = Sample("s2", "/d/b.fq.gz", "c", None, md5_r2="e" * 32)    # mate1_only
    assert K.expected_md5s([s1, s2], man) == {"/d/a.fq.gz": "f" * 32, "/d/b.fq.gz": "2" * 32}


def test_check_md5s(tmp_path):
    f = tmp_path / "a"
    f.write_bytes(b"data")
    assert K.check_md5s({str(f): _md5(f)}) == []
    assert "does not match" in K.check_md5s({str(f): "0" * 32})[0]


def test_samplesheet_md5_columns(tmp_path):
    ss = tmp_path / "s.tsv"
    ss.write_text("sample_id\tfastq_r1\tcondition\tmd5_r1\ns1\ta\tc\t" + "A" * 32 + "\n"
                  "s2\tb\tc\tnope\n")
    samples = read_samplesheet(ss)
    assert samples[0].md5_r1 == "a" * 32
    assert any("md5_r1" in p for p in check_samplesheet(samples, check_files=False))


def _sample(tmp_path, md5=None):
    raw = tmp_path / "raw"
    raw.mkdir()
    s = Sample("s1", write_fastq_gz(raw / "A_R1.fq.gz", 10), "c",
               write_fastq_gz(raw / "A_R2.fq.gz", 10))
    out = tmp_path / "out" / "t"
    for d in ("02_trimmed", "04_align"):
        (out / d).mkdir(parents=True)
    md5s = {s.fastq_r1: md5 or _md5(s.fastq_r1), s.fastq_r2: _md5(s.fastq_r2)}
    return s, out, md5s


def test_run_rejects_a_changed_raw_fastq(tmp_path):
    s, out, md5s = _sample(tmp_path, md5="0" * 32)
    with pytest.raises(R.SampleError) as e:
        R._process_sample(s, True, out, {"index_prefix": "i", "stamp": {}}, 2,
                          CallableRunner(FakeTools()), md5s=md5s)
    assert e.value.step == "md5" and "does not match" in str(e.value)
    assert not (out / "04_align/s1.done.json").exists()


def test_run_checks_md5_only_when_processing(tmp_path, monkeypatch):
    s, out, md5s = _sample(tmp_path)
    calls = []
    real = R.check_md5s
    monkeypatch.setattr(R, "check_md5s", lambda exp: calls.append(exp) or real(exp))
    b = {"index_prefix": "i", "stamp": {}}
    assert R._process_sample(s, True, out, b, 2, CallableRunner(FakeTools()), md5s=md5s)["resumed"] is False
    assert R._process_sample(s, True, out, b, 2, CallableRunner(FakeTools()), md5s=md5s)["resumed"] is True
    assert len(calls) == 1


def test_validate_cli_checks_md5(tmp_path):
    from tests.unit.test_cli import _bin
    import yaml
    fq = tmp_path / "a.fq.gz"
    fq.write_bytes(b"reads")
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "r", "reference": {"species": "mabs"},
                                   "contrasts": {"explicit": [{"name": "B_vs_A", "numerator": "B",
                                                               "denominator": "A"}]}}))
    ss = tmp_path / "s.tsv"
    ss.write_text(f"sample_id\tfastq_r1\tcondition\na1\t{fq}\tA\nb1\t{fq}\tB\n")
    good, bad = tmp_path / "good.md5", tmp_path / "bad.md5"
    good.write_text(f"{_md5(fq)}  a.fq.gz\n")
    bad.write_text(f"{'0' * 32}  a.fq.gz\n")
    ok = _bin(["validate", str(cfg), "--samplesheet", str(ss), "--md5", str(good)], tmp_path)
    assert ok.returncode == 0 and "md5 ok" in ok.stdout, ok.stderr
    no = _bin(["validate", str(cfg), "--samplesheet", str(ss), "--md5", str(bad)], tmp_path)
    assert no.returncode == 2 and "does not match" in no.stderr
