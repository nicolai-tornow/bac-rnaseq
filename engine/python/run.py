from __future__ import annotations
import hashlib
import json
import os
import shutil
import signal
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import yaml
from .samplesheet import read_samplesheet, is_paired, apply_layout, check_samplesheet
from .contrasts import expand_contrasts
from .build_refs import build_bundle, md5
from . import commands as C
from . import qc_triage, run_report, siteconfig
from .procs import Cancelled, ProcRunner, Terminated, as_runner
from .runlock import RunLock

STRANDS = ("reverse", "forward", "unstranded")
CORES_PER_SAMPLE = 8      # default split of the thread budget across parallel samples
SORT_MEM = "256M"         # samtools sort memory per thread (4 x 8 threads ~ 8 GB)


def _check(runner, cmd):
    res = runner.run(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"stage failed: {cmd[0]}\n{res.stderr}")
    return res


def _run_deseq2(runner, r_script, counts, coldata, contrasts_tsv, out_dir, use_batch):
    _check(runner, ["Rscript", str(r_script), str(counts), str(coldata),
                    str(contrasts_tsv), str(out_dir), str(use_batch)])


def _threads(config) -> int:
    return config.resources.threads or siteconfig.read_site().get("threads") or 4


def _parallelism(config, total, n_samples):
    """(samples at once, threads per sample). Default: one sample per 8 threads."""
    par = (config.resources.parallel_samples or siteconfig.read_site().get("parallel_samples")
           or max(1, total // CORES_PER_SAMPLE))
    par = max(1, min(par, n_samples))
    return par, max(1, total // par)


def _bundle_dir(config, work_dir) -> Path:
    """Built bundles live in site.yaml `refs_root` (shared, reused) when set."""
    ref = config.reference
    name = ref.species
    if ref.species == "custom":
        key = json.dumps([ref.fasta, ref.gff, ref.feature_types, ref.id_attribute,
                          ref.exclude_seqids, ref.seqid_map, ref.structural_rna])
        name = "custom-" + hashlib.md5(key.encode()).hexdigest()[:10]
    root = siteconfig.read_site().get("refs_root")
    return Path(root) / name if root else Path(work_dir) / "out" / "refs" / name


def _count_reads(runner, fastq):
    res = runner.pipe(["gzip", "-cdf", str(fastq)], ["wc", "-l"])
    if res.rc1 != 0 or res.rc2 != 0:
        raise RuntimeError(f"{fastq}: cannot be read as a (gzipped) FASTQ: "
                           f"{(res.stderr1 or res.stderr).strip()}")
    txt = res.stdout.strip()
    return int(txt) // 4 if txt else None


class SampleError(RuntimeError):
    """A sample failed at `step` (fastp, bowtie2, sort, index, promote, integrity)."""

    def __init__(self, sample_id, step, message):
        super().__init__(f"{sample_id}: {step} failed: {message}")
        self.sample_id, self.step, self.message = sample_id, step, message


_SIGPIPE = (-signal.SIGPIPE, 128 + signal.SIGPIPE)


@lru_cache(maxsize=1)
def _plugin_commit():
    return run_report.plugin_version().get("commit")


def _final_paths(out, sid) -> dict:
    t, a = Path(out) / "02_trimmed", Path(out) / "04_align"
    return {"r1": t / f"{sid}_R1.fq.gz", "r2": t / f"{sid}_R2.fq.gz", "json": t / f"{sid}.json",
            "html": t / f"{sid}.html", "bam": a / f"{sid}.bam", "bai": a / f"{sid}.bam.bai",
            "log": a / f"{sid}.bowtie2.log", "marker": a / f"{sid}.done.json",
            "sam": a / f"{sid}.sam"}


def _remove_outputs(fin):
    """Marker first: a sample without a marker is never trusted."""
    for k in ("marker", "bam", "bai", "log", "sam", "r1", "r2", "json", "html"):
        fin[k].unlink(missing_ok=True)


def _fastq_id(path):
    if not path:
        return None
    rp = os.path.realpath(path)
    try:
        st = os.stat(rp)
        return {"path": rp, "size": st.st_size, "mtime": st.st_mtime}
    except OSError:
        return {"path": rp, "size": None, "mtime": None}


def _sample_inputs(s, paired, bundle) -> dict:
    """What a sample's BAM was made from; any change makes the sample re-run."""
    return {"fastq_r1": _fastq_id(s.fastq_r1), "fastq_r2": _fastq_id(s.fastq_r2),
            "paired": paired,
            "reference_fasta_md5": (bundle.get("stamp") or {}).get("fasta_md5"),
            "fastp_args": list(C.FASTP_SETTINGS), "bowtie2_args": list(C.BOWTIE2_SETTINGS)}


def _marker_valid(fin, inputs) -> tuple[bool, str]:
    if not fin["marker"].exists():
        return False, "no completion marker"
    try:
        m = json.loads(fin["marker"].read_text())
    except ValueError:
        return False, "unreadable completion marker"
    if m.get("schema") != 1:
        return False, "unknown completion marker schema"
    if m.get("inputs") != inputs:
        return False, "inputs changed since the sample was processed"
    if not fin["bam"].exists():
        return False, "BAM missing"
    if _fastp_stats(fin["json"]).get("reads_passed") != m.get("reads_passed"):
        return False, "fastp report does not match the completion marker"
    if md5(fin["bam"]) != m.get("bam_md5"):
        return False, "BAM changed since the sample was processed"
    return True, ""


def _tail(path, n=15) -> str:
    p = Path(path)
    return "\n".join(p.read_text(errors="replace").splitlines()[-n:]) if p.exists() else ""


def _process_sample(s, paired, out, bundle, threads, runner, cleaned=None):
    """fastp -> bowtie2 | samtools sort into staging folders; outputs are moved into
    place only after every check passed, and the completion marker is written last.
    A sample with a valid marker is skipped (its trimmed reads need not exist)."""
    sid = s.sample_id
    fin = _final_paths(out, sid)
    inputs = _sample_inputs(s, paired, bundle)
    ok, why = _marker_valid(fin, inputs)
    if ok:
        if not fin["bai"].exists():
            _check(runner, ["samtools", "index", str(fin["bam"])])
        return {"resumed": True, **_fastp_stats(fin["json"])}

    _remove_outputs(fin)
    tstage = out / "02_trimmed" / f".{sid}.partial"
    astage = out / "04_align" / f".{sid}.partial"
    for d in (tstage, astage):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
    step, done = "fastp", False
    try:
        t1 = tstage / f"{sid}_R1.fq.gz"
        t2 = tstage / f"{sid}_R2.fq.gz" if paired else None
        fjson = tstage / f"{sid}.json"
        # fastp can stop early on a truncated/corrupt input and still exit 0, so
        # compare what it read with the raw read count.
        raw = _count_reads(runner, s.fastq_r1)
        _check(runner, C.fastp_cmd(s.fastq_r1, str(t1), threads,
                                   r2=s.fastq_r2, out2=str(t2) if t2 else None,
                                   json=str(fjson), html=str(fjson.with_suffix(".html"))))
        stats = _fastp_stats(fjson)
        expected = raw * (2 if paired else 1) if raw is not None else None
        if expected is not None and stats.get("reads_in") is not None and stats["reads_in"] != expected:
            raise RuntimeError(f"fastp read {stats['reads_in']} reads but the input has "
                               f"{expected}; the FASTQ may be truncated or corrupt")
        stats["raw_reads"] = expected

        step = "bowtie2"
        log, bam = astage / f"{sid}.bowtie2.log", astage / f"{sid}.bam"
        bt2 = C.bowtie2_cmd(bundle["index_prefix"], str(t1), threads, r2=str(t2) if t2 else None)
        sort = ["samtools", "sort", "-@", str(min(4, max(1, threads // 4))), "-m", SORT_MEM,
                "-T", str(astage / "sort"), "-o", str(bam), "-"]
        res = runner.pipe(bt2, sort, stderr1=str(log))
        if res.rc1 != 0 and not (res.rc2 != 0 and res.rc1 in _SIGPIPE):
            raise RuntimeError(f"bowtie2 exited {res.rc1}:\n{_tail(log)}")
        if res.rc2 != 0:
            step = "sort"
            raise RuntimeError(f"samtools sort exited {res.rc2}: {res.stderr.strip()[-800:]}")

        step = "index"
        _check(runner, ["samtools", "index", str(bam)])
        _check(runner, ["samtools", "quickcheck", str(bam)])
        records = int(_check(runner, ["samtools", "view", "-c", "-@", str(threads),
                                      str(bam)]).stdout.strip())

        step = "promote"
        for src in (t1, t2, fjson, fjson.with_suffix(".html")):
            if src is not None and src.exists():
                os.replace(src, out / "02_trimmed" / src.name)
        for src in (bam, Path(f"{bam}.bai"), log):
            os.replace(src, out / "04_align" / src.name)
        marker = {"schema": 1, "sample_id": sid,
                  "created": datetime.now(timezone.utc).isoformat(),
                  "plugin_commit": _plugin_commit(), "bam_md5": md5(fin["bam"]),
                  "bam_records": records, "reads_passed": stats.get("reads_passed"),
                  "inputs": inputs}
        tmp = astage / "done.json"
        tmp.write_text(json.dumps(marker, indent=1))
        os.replace(tmp, fin["marker"])
        done = True
        return {"resumed": False, "resume_skipped": why, **stats}
    except Cancelled:
        raise
    except Exception as e:
        raise SampleError(sid, step, str(e)) from e
    finally:
        for d in (tstage, astage):
            if d.exists():
                shutil.rmtree(d, ignore_errors=True)
                if not done and cleaned is not None:
                    cleaned.append(str(d.relative_to(out)))


def _sample_integrity(sid, paired, out, runner, threads) -> list[str]:
    """Trimmed FASTQs (when still present) decompress and hold what fastp passed;
    the BAM holds the records its completion marker recorded."""
    fin = _final_paths(out, sid)
    try:
        m = json.loads(fin["marker"].read_text())
    except (OSError, ValueError):
        return ["completion marker missing or unreadable"]
    problems = []
    rp = m.get("reads_passed")
    want = (rp // 2 if paired else rp) if isinstance(rp, int) else None
    for k in ("r1", "r2") if paired else ("r1",):
        f = fin[k]
        if not f.exists():          # cleaned up; the BAM check still applies
            continue
        res = runner.pipe(["gzip", "-cd", str(f)], ["wc", "-l"])
        if res.rc1 != 0 or res.rc2 != 0:
            problems.append(f"{f.name}: not a valid gzip stream")
            continue
        n = int(res.stdout.strip() or 0) // 4
        if want is not None and n != want:
            problems.append(f"{f.name}: {n} reads, but fastp passed {want}")
    res = runner.run(["samtools", "view", "-c", "-@", str(threads), str(fin["bam"])])
    n = res.stdout.strip()
    if res.returncode != 0 or not n.isdigit() or int(n) != m.get("bam_records"):
        problems.append(f"{fin['bam'].name}: {n or 'unreadable'} records, but the "
                        f"completion marker says {m.get('bam_records')}")
    return problems


def _fastp_stats(fjson) -> dict:
    p = Path(fjson)
    if not p.exists():
        return {}
    s = json.loads(p.read_text())["summary"]
    return {"reads_in": s["before_filtering"]["total_reads"],
            "reads_passed": s["after_filtering"]["total_reads"]}


def _featurecounts(runner, tmp, cleaned=None, **kw):
    """featureCounts with its temp files in a private folder that is removed
    whether it succeeds or not (a crash would otherwise leave GBs of temp-core-*)."""
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    ok = False
    try:
        _check(runner, C.featurecounts_cmd(tmp_dir=str(tmp), **kw))
        ok = True
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        if not ok and cleaned is not None:
            cleaned.append(f"{tmp.parent.name}/{tmp.name}")


def _strand_check(runner, saf, bams, threads, paired, out_dir, declared, main_summary,
                  tmp=None, cleaned=None):
    """Assigned fraction per sample at all three -s settings."""
    per = {declared: qc_triage.parse_featurecounts_summary(main_summary)}
    out_dir.mkdir(parents=True, exist_ok=True)
    for st in STRANDS:
        if st == declared:
            continue
        fc = out_dir / f"fc_{st}.txt"
        _featurecounts(runner, tmp or out_dir / ".featurecounts.partial", cleaned,
                       saf=saf, out=str(fc), bams=bams, threads=threads,
                       strandedness=st, paired=paired)
        per[st] = qc_triage.parse_featurecounts_summary(str(fc) + ".summary")
    samples = per[declared].keys()
    return {s: {st: per[st][s]["assigned_frac"] for st in STRANDS} for s in samples}


def _reshape_counts(fc_txt, out_tsv, keep_ids=None, other_tsv=None):
    """featureCounts output -> Gene x sample table. With keep_ids, rows not in it
    (the added structural RNAs) go to other_tsv instead."""
    lines = [l for l in Path(fc_txt).read_text().splitlines() if not l.startswith("#")]
    header = lines[0].split("\t")
    keep = list(range(6, len(header)))  # sample columns follow Geneid,Chr,Start,End,Strand,Length
    samp = [qc_triage.sample_name(header[i]) for i in keep]
    head = "Gene\t" + "\t".join(samp) + "\n"
    main, other = [head], [head]
    for l in lines[1:]:
        f = l.split("\t")
        row = f[0] + "\t" + "\t".join(f[i] for i in keep) + "\n"
        (main if keep_ids is None or f[0] in keep_ids else other).append(row)
    Path(out_tsv).write_text("".join(main))
    if other_tsv is not None:
        Path(other_tsv).write_text("".join(other))
    return out_tsv


def _collect_qc(align_logs, fc_summary, strand_fracs, ncrna, declared):
    fc = qc_triage.parse_featurecounts_summary(fc_summary)
    out = {}
    for sid, log in align_logs.items():
        pct = qc_triage.parse_bowtie2_log(Path(log).read_text())
        strand = qc_triage.infer_strandedness(strand_fracs[sid], declared)
        nc = ncrna.get(sid, {"total": 0.0, "by_class": {}})
        v = qc_triage.triage_sample(pct, fc[sid]["assigned_frac"], nc["total"],
                                    strand=strand, nofeature_frac=fc[sid]["nofeature_frac"])
        v.update({"alignment_pct": pct, "assigned_frac": fc[sid]["assigned_frac"],
                  "nofeature_frac": fc[sid]["nofeature_frac"],
                  "ncrna_frac": nc["total"], "ncrna_by_class": nc["by_class"],
                  "strandedness": strand})
        out[sid] = v
    return out


def _ncrna_fracs(fc_txt, structural_classes):
    import pandas as pd
    df = pd.read_csv(fc_txt, sep="\t", comment="#", index_col=0).iloc[:, 5:]
    df.columns = [qc_triage.sample_name(c) for c in df.columns]
    return qc_triage.ncrna_fractions(df, structural_classes)


def _de_summary(results_dir, contrasts, padj, log2fc):
    import pandas as pd
    out = {}
    for c in contrasts:
        p = Path(results_dir) / f"{c.name}.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t")
        sig = df[df["padj"] < padj]
        big = sig[sig["log2FoldChange"].abs() >= log2fc]
        out[c.name] = {"tested": int(df["padj"].notna().sum()),
                       f"padj<{padj}": int(len(sig)),
                       "up": int((big["log2FoldChange"] > 0).sum()),
                       "down": int((big["log2FoldChange"] < 0).sum()),
                       "thresholds": {"padj": padj, "log2fc": log2fc}}
    return out


_STAGING = ("02_trimmed/.*.partial", "04_align/.*.partial",
            "05_counts/.featurecounts.partial", "06_deseq/.partial")
_NO_HANDLER = object()


def _sweep_staging(out) -> list[str]:
    """Remove staging folders a killed run left behind."""
    removed = []
    for pat in _STAGING:
        for p in sorted(Path(out).glob(pat)):
            shutil.rmtree(p, ignore_errors=True)
            removed.append(str(p.relative_to(out)))
    return removed


def _install_sigterm():
    """SIGTERM raises Terminated in the main thread, so the run unwinds like on
    Ctrl-C. Only the main thread can install a handler."""
    if threading.current_thread() is not threading.main_thread():
        return _NO_HANDLER
    fired = []

    def handler(signum, frame):
        if not fired:
            fired.append(signum)
            raise Terminated(f"received signal {signum}")

    return signal.signal(signal.SIGTERM, handler)


def _restore_sigterm(prev):
    if prev is not _NO_HANDLER:
        signal.signal(signal.SIGTERM, signal.SIG_DFL if prev is None else prev)


def _promote_dir(src, dst):
    """Move every entry of src into dst, replacing what is there."""
    for p in list(Path(src).iterdir()):
        target = Path(dst) / p.name
        if target.is_dir() and not target.is_symlink():
            shutil.rmtree(target)
        elif target.exists() or target.is_symlink():
            target.unlink()
        os.replace(p, target)


def _trim_align(samples, paired, out, bundle, threads, par, runner, cleaned):
    """Samples in parallel; the first failure cancels the rest (queued ones never
    start, running ones are killed and remove their staging)."""
    pool = ThreadPoolExecutor(max_workers=par)
    futures = {pool.submit(_process_sample, s, paired, out, bundle, threads, runner,
                           cleaned=cleaned): s.sample_id for s in samples}
    reads = {}
    try:
        for f in as_completed(futures):
            reads[futures[f]] = f.result()
    except BaseException:
        runner.cancel()
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    return {s.sample_id: reads[s.sample_id] for s in samples}


def _write_failed_report(out, config, state, exc):
    failure = {"stage": state["stage"],
               "step": getattr(exc, "step", None),
               "sample": getattr(exc, "sample_id", None),
               "error": f"{type(exc).__name__}: {exc}"[:4000],
               "cleaned": sorted(set(state["cleaned"]))}
    try:
        rep = run_report.build_report(config.run_name, state["params"], {}, {}, [], {},
                                      "failed", provenance=state["provenance"],
                                      warnings=state["warnings"], failure=failure)
        run_report.write_report(rep, out / "00_run_report.json")
    except Exception as e:           # never hide the original error
        exc.add_note(f"could not write the failed run report: {e}")


def run_pipeline(config, work_dir, refs_root, samplesheet_path, runner=None,
                 allow_qc_fail=False, check_files=True):
    """Pre-flight errors and a live lock raise without touching the run folder's
    report. Once the run holds the lock, any failure (incl. Ctrl-C and SIGTERM)
    kills running tools, removes staging, writes a `failed` report and re-raises."""
    work_dir = Path(work_dir)
    out = work_dir / "out" / config.run_name
    samples = apply_layout(read_samplesheet(samplesheet_path), config.reads.layout)
    paired = is_paired(samples)
    conds = [s.condition for s in samples]
    cons = expand_contrasts(config.contrasts.explicit, config.contrasts.all_vs_all, conds)
    problems = check_samplesheet(samples, cons, check_files=check_files)
    if problems:
        raise ValueError("sample sheet problems:\n  " + "\n  ".join(problems))

    runner = ProcRunner() if runner is None else as_runner(runner)
    out.mkdir(parents=True, exist_ok=True)
    lock = RunLock(out)
    stale = lock.acquire()
    state = {"stage": "setup", "cleaned": [], "params": {}, "warnings": [],
             "provenance": {"plugin": run_report.plugin_version()}}
    if stale is not None:
        state["provenance"]["stale_lock_cleared"] = stale
    prev = _install_sigterm()
    try:
        swept = _sweep_staging(out)
        if swept:
            state["provenance"]["stale_staging_removed"] = swept
        return _run_stages(config, work_dir, out, refs_root, samplesheet_path, samples,
                           paired, cons, runner, allow_qc_fail, state)
    except BaseException as e:
        runner.cancel()
        state["cleaned"] += _sweep_staging(out)
        _write_failed_report(out, config, state, e)
        raise
    finally:
        _restore_sigterm(prev)
        lock.release()


def _run_stages(config, work_dir, out, refs_root, samplesheet_path, samples, paired, cons,
                runner, allow_qc_fail, state):
    for d in ["00_inputs", "01_qc_raw", "02_trimmed", "03_qc_trimmed", "04_align",
              "05_counts", "06_deseq", "qc"]:
        (out / d).mkdir(parents=True, exist_ok=True)
    strand = config.reference.strandedness
    threads = _threads(config)
    par, per_sample = _parallelism(config, threads, len(samples))
    params = {"strandedness": strand, "paired": paired, "layout": config.reads.layout,
              "threads": threads, "parallel_samples": par,
              "threads_per_sample": per_sample, "allow_qc_fail": allow_qc_fail}
    state["params"] = params
    provenance = state["provenance"]
    provenance["inputs"] = {"config": "00_inputs/config.yaml",
                            "samplesheet": "00_inputs/samplesheet.tsv"}

    # Provenance: exactly what this run was given.
    state["stage"] = "inputs"
    (out / "00_inputs" / "config.yaml").write_text(
        yaml.safe_dump(config.model_dump(), sort_keys=False))
    shutil.copyfile(samplesheet_path, out / "00_inputs" / "samplesheet.tsv")

    state["stage"] = "refs"
    ref = config.reference
    bundle = build_bundle(ref.species, refs_root=refs_root,
                          out_dir=_bundle_dir(config, work_dir), threads=threads,
                          fasta=ref.fasta, gff=ref.gff, feature_types=tuple(ref.feature_types),
                          id_attribute=ref.id_attribute, exclude_seqids=ref.exclude_seqids,
                          seqid_map=ref.seqid_map, structural_rna=ref.structural_rna)
    provenance["reference"] = {
        "bundle_dir": str(_bundle_dir(config, work_dir)),
        "reused": bundle.get("reused", False),
        "fasta_md5": md5(bundle["fasta"]) if Path(bundle["fasta"]).exists() else None,
        "gff_md5": bundle.get("stamp", {}).get("gff_md5"),
        "structural_rna_md5": bundle.get("stamp", {}).get("structural_rna_md5"),
        "saf_md5": bundle.get("saf_md5")}

    state["stage"] = "fastqc_raw"
    fastqs = [f for s in samples for f in (s.fastq_r1, s.fastq_r2) if f]
    _check(runner, C.fastqc_cmd(fastqs, str(out / "01_qc_raw"), threads))

    # Trim + align several samples at once, splitting the thread budget.
    state["stage"] = "trim_align"
    reads = _trim_align(samples, paired, out, bundle, per_sample, par, runner, state["cleaned"])
    provenance["reads"] = reads

    # Before counting, prove every sample's outputs are intact; re-process once.
    state["stage"] = "integrity"
    for s in samples:
        problems = _sample_integrity(s.sample_id, paired, out, runner, threads)
        if not problems:
            continue
        _remove_outputs(_final_paths(out, s.sample_id))
        reads[s.sample_id] = _process_sample(s, paired, out, bundle, threads, runner,
                                             cleaned=state["cleaned"])
        reads[s.sample_id]["reprocessed_after_integrity"] = problems
        again = _sample_integrity(s.sample_id, paired, out, runner, threads)
        if again:
            raise SampleError(s.sample_id, "integrity", "; ".join(again))
    align_logs = {s.sample_id: str(out / "04_align" / f"{s.sample_id}.bowtie2.log")
                  for s in samples}
    bams = [str(out / "04_align" / f"{s.sample_id}.bam") for s in samples]

    state["stage"] = "fastqc_trimmed"
    trimmed = sorted(str(p) for p in (out / "02_trimmed").glob("*.fq.gz"))
    if trimmed:
        _check(runner, C.fastqc_cmd(trimmed, str(out / "03_qc_trimmed"), threads))

    state["stage"] = "featurecounts"
    fc_tmp = out / "05_counts" / ".featurecounts.partial"
    fc = out / "05_counts" / "featurecounts.txt"
    _featurecounts(runner, fc_tmp, state["cleaned"], saf=bundle["saf"], out=str(fc),
                   bams=bams, threads=threads, strandedness=strand, paired=paired)
    _reshape_counts(fc, out / "05_counts" / "counts.tsv",
                    keep_ids=set(bundle.get("gene_ids") or []) or None,
                    other_tsv=out / "05_counts" / "ncrna_counts.tsv")
    state["stage"] = "strand_check"
    strand_fracs = _strand_check(runner, bundle["saf"], bams, threads, paired,
                                 out / "05_counts" / "strand_check", strand,
                                 str(fc) + ".summary", tmp=fc_tmp, cleaned=state["cleaned"])
    state["stage"] = "qc"
    ncrna = _ncrna_fracs(fc, bundle.get("structural", {}))
    qc = _collect_qc(align_logs, str(fc) + ".summary", strand_fracs, ncrna, strand)
    state["stage"] = "multiqc"
    _check(runner, C.multiqc_cmd(str(out), str(out / "qc" / "multiqc")))

    invariants = {"strandedness": strand, "n_features": bundle["n_features"],
                  "seqids": bundle["seqids"],
                  "structural_rnas_added": bundle.get("extra_ids", [])}
    outputs = {"counts": "05_counts/counts.tsv", "ncrna_counts": "05_counts/ncrna_counts.tsv",
               "strand_check": "05_counts/strand_check", "multiqc": "qc/multiqc"}

    # QC gate: a FAIL halts before DESeq2 so bad data cannot silently produce a DE table.
    failed = [sid for sid, v in qc.items() if v.get("verdict") == "FAIL"]
    if failed and not allow_qc_fail:
        rep = run_report.build_report(config.run_name, params, qc, invariants, [],
                                      outputs, "qc_fail", provenance=provenance,
                                      warnings=state["warnings"])
        run_report.write_report(rep, out / "00_run_report.json")
        return rep

    # DESeq2 writes into a staging folder; results replace the old ones only on success.
    state["stage"] = "deseq2"
    dstage = out / "06_deseq" / ".partial"
    shutil.rmtree(dstage, ignore_errors=True)
    dstage.mkdir(parents=True)
    ok = False
    try:
        use_batch = 1 if (config.design.batch_variable and any(s.batch for s in samples)) else 0
        coldata = dstage / "coldata.tsv"
        hdr = "sample\tcondition" + ("\tbatch" if use_batch else "") + "\n"
        coldata.write_text(hdr + "".join(
            f"{s.sample_id}\t{s.condition}" + (f"\t{s.batch}" if use_batch else "") + "\n"
            for s in samples))
        ctsv = dstage / "contrasts.tsv"
        ctsv.write_text("".join(f"{c.name}\t{c.numerator}\t{c.denominator}\n" for c in cons))
        _run_deseq2(runner, Path(__file__).parents[1] / "r" / "deseq2.R",
                    out / "05_counts" / "counts.tsv", coldata, ctsv, dstage, use_batch)
        _promote_dir(dstage, out / "06_deseq")
        ok = True
    finally:
        shutil.rmtree(dstage, ignore_errors=True)
        if not ok:
            state["cleaned"].append("06_deseq/.partial")

    outputs["deseq"] = "06_deseq/results"
    rep = run_report.build_report(
        config.run_name, params, qc, invariants, [c.name for c in cons], outputs, "ok",
        provenance=provenance, warnings=state["warnings"],
        de_summary=_de_summary(out / "06_deseq" / "results", cons,
                               config.thresholds.padj, config.thresholds.log2fc))
    run_report.write_report(rep, out / "00_run_report.json")
    return rep
