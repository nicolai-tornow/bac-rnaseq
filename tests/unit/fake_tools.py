"""Stand-ins for the pipeline's tools that write the files the real tools write, so
the orchestration (staging, promotion, markers, clean-up) runs for real in tests."""
import gzip
import json
import re
from pathlib import Path
from engine.python.procs import Result


def write_fastq(path, n, seq="GGCCAATT", tag="u"):
    with open(path, "w") as fh:
        for i in range(n):
            fh.write(f"@{tag}{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    return str(path)


def write_fastq_gz(path, n, seq="ACGTACGTACGTACGTACGTACGTACGTACGTACGT", tag="r"):
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            fh.write(f"@{tag}{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    return str(path)


def _arg(cmd, flag):
    return cmd[cmd.index(flag) + 1] if flag in cmd else None


class FakeTools:
    """Callable like subprocess.run. fail={key: returncode}; hooks={key: fn(cmd)} run
    first (a hook returning a Result replaces the tool). Keys: the executable, or
    "samtools <subcommand>"."""

    def __init__(self, align_pct=96.5, bam_records=100, fail=None, hooks=None):
        self.align_pct = align_pct
        self.bam_records = bam_records
        self.fail = dict(fail or {})
        self.hooks = dict(hooks or {})
        self.calls = []

    @staticmethod
    def key(cmd):
        exe = Path(cmd[0]).name
        return f"samtools {cmd[1]}" if exe == "samtools" else exe

    def called(self, key):
        return [c for c in self.calls if self.key(c) == key]

    def __call__(self, cmd, input=None, **kw):
        self.calls.append(list(cmd))
        k = self.key(cmd)
        if k in self.hooks:
            r = self.hooks[k](cmd)
            if isinstance(r, Result):
                return r
        if k in self.fail:
            return Result(self.fail[k], "", f"{k} failed")
        fn = getattr(self, "_" + re.sub(r"\W", "_", k), None)
        return fn(cmd, input) if fn else Result(0)

    def _fastp(self, cmd, input):
        n = 0
        for i, o in (("--in1", "--out1"), ("--in2", "--out2")):
            src, dst = _arg(cmd, i), _arg(cmd, o)
            if src:
                data = gzip.open(src, "rb").read()
                Path(dst).write_bytes(gzip.compress(data))
                n += data.count(b"\n") // 4
        if _arg(cmd, "--json"):
            Path(_arg(cmd, "--json")).write_text(json.dumps({"summary": {
                "before_filtering": {"total_reads": n}, "after_filtering": {"total_reads": n}}}))
        if _arg(cmd, "--html"):
            Path(_arg(cmd, "--html")).write_text("<html/>")
        return Result(0)

    def _bowtie2(self, cmd, input):
        if _arg(cmd, "--un-conc"):             # 3 unaligned pairs, uncompressed
            for m in ("1", "2"):
                write_fastq(_arg(cmd, "--un-conc").replace("%", m), 3)
        if _arg(cmd, "--un"):
            write_fastq(_arg(cmd, "--un"), 3)
        return Result(0, "", f"{self.align_pct:.2f}% overall alignment rate\n")

    def _samtools_sort(self, cmd, input):
        Path(_arg(cmd, "-o")).write_text(f"FAKEBAM {self.bam_records}\n")
        return Result(0)

    def _samtools_index(self, cmd, input):
        Path(cmd[-1] + ".bai").write_text("bai")
        return Result(0)

    def _samtools_quickcheck(self, cmd, input):
        return Result(0 if Path(cmd[-1]).exists() else 1)

    def _samtools_view(self, cmd, input):
        return Result(0, Path(cmd[-1]).read_text().split()[1] + "\n")

    def _gzip(self, cmd, input):
        try:
            return Result(0, gzip.open(cmd[-1], "rt").read())
        except (OSError, EOFError) as e:
            return Result(1, "", str(e))

    def _wc(self, cmd, input):
        return Result(0, f"{(input or '').count(chr(10))}\n")

    def _featureCounts(self, cmd, input):
        out = _arg(cmd, "-o")
        bams = [a for a in cmd if a.endswith(".bam")]
        cols = "\t".join(bams)
        Path(out).write_text("# fake\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t" + cols +
                             "\nG1\tc\t1\t10\t+\t10\t" + "\t".join("5" for _ in bams) + "\n")
        Path(out + ".summary").write_text("Status\t" + cols + "\nAssigned\t" +
                                          "\t".join("5" for _ in bams) + "\n")
        return Result(0)

    def _Rscript(self, cmd, input):
        out_dir = Path(cmd[5])
        (out_dir / "results").mkdir(parents=True, exist_ok=True)
        for line in Path(cmd[4]).read_text().splitlines():
            (out_dir / "results" / f"{line.split(chr(9))[0]}.tsv").write_text(
                "Gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\nG1\t1\t0\t1\t0\t1\t1\n")
        (out_dir / "normalized_counts.tsv").write_text("Gene\n")
        return Result(0)
