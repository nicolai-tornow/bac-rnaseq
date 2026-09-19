from __future__ import annotations
import hashlib
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import yaml
from .samplesheet import read_samplesheet, is_paired, apply_layout, check_samplesheet
from .contrasts import expand_contrasts
from .build_refs import build_bundle, md5
from . import commands as C
from . import qc_triage, run_report, siteconfig

STRANDS = ("reverse", "forward", "unstranded")
CORES_PER_SAMPLE = 8      # default split of the thread budget across parallel samples
SORT_MEM = "256M"         # samtools sort memory per thread (4 x 8 threads ~ 8 GB)


def _check(runner, cmd, **kw):
    res = runner(cmd, capture_output=True, text=True, **kw)
    if getattr(res, "returncode", 0) != 0:
        raise RuntimeError(f"stage failed: {cmd[0]}\n{getattr(res, 'stderr', '')}")
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
    res = runner(["bash", "-c", 'gzip -cdf "$1" | wc -l', "_", str(fastq)],
                 capture_output=True, text=True)
    txt = (getattr(res, "stdout", "") or "").strip()
    return int(txt) // 4 if txt else None


def _process_sample(s, paired, out, bundle, threads, runner):
    """fastp -> bowtie2 -> sorted BAM. Skips a sample whose BAM already passes
    samtools quickcheck (resume after a qc_fail or an interrupted run)."""
    t1 = out / "02_trimmed" / f"{s.sample_id}_R1.fq.gz"
    t2 = out / "02_trimmed" / f"{s.sample_id}_R2.fq.gz" if paired else None
    fjson = out / "02_trimmed" / f"{s.sample_id}.json"
    bam = out / "04_align" / f"{s.sample_id}.bam"
    log = out / "04_align" / f"{s.sample_id}.bowtie2.log"
    if bam.exists() and log.exists() and fjson.exists() and \
            getattr(runner(["samtools", "quickcheck", str(bam)], capture_output=True,
                           text=True), "returncode", 1) == 0:
        return {"resumed": True, **_fastp_stats(fjson)}

    # fastp can stop early on a truncated/corrupt input and still exit 0, so
    # compare what it read with the raw read count.
    raw = _count_reads(runner, s.fastq_r1)
    _check(runner, C.fastp_cmd(s.fastq_r1, str(t1), threads,
                               r2=s.fastq_r2, out2=str(t2) if t2 else None,
                               json=str(fjson), html=str(fjson.with_suffix(".html"))))
    stats = _fastp_stats(fjson)
    expected = raw * (2 if paired else 1) if raw is not None else None
    if expected is not None and stats.get("reads_in") is not None and stats["reads_in"] != expected:
        raise RuntimeError(f"{s.sample_id}: fastp read {stats['reads_in']} reads but the "
                           f"input has {expected}; the FASTQ may be truncated or corrupt")
    stats["raw_reads"] = expected

    sam = out / "04_align" / f"{s.sample_id}.sam"
    bt2 = C.bowtie2_cmd(bundle["index_prefix"], str(t1), threads, r2=str(t2) if t2 else None)
    res = runner(bt2 + ["-S", str(sam)], capture_output=True, text=True)
    Path(log).write_text(getattr(res, "stderr", "") or "")
    if getattr(res, "returncode", 0) != 0:
        raise RuntimeError(f"stage failed: bowtie2 ({s.sample_id}); see {log}")
    _check(runner, ["samtools", "sort", "-@", str(threads), "-m", SORT_MEM,
                    "-o", str(bam), str(sam)])
    _check(runner, ["samtools", "index", str(bam)])
    _check(runner, ["samtools", "quickcheck", str(bam)])
    sam.unlink(missing_ok=True)
    return {"resumed": False, **stats}


def _fastp_stats(fjson) -> dict:
    p = Path(fjson)
    if not p.exists():
        return {}
    s = json.loads(p.read_text())["summary"]
    return {"reads_in": s["before_filtering"]["total_reads"],
            "reads_passed": s["after_filtering"]["total_reads"]}


def _strand_check(runner, saf, bams, threads, paired, out_dir, declared, main_summary):
    """Assigned fraction per sample at all three -s settings."""
    per = {declared: qc_triage.parse_featurecounts_summary(main_summary)}
    out_dir.mkdir(parents=True, exist_ok=True)
    for st in STRANDS:
        if st == declared:
            continue
        fc = out_dir / f"fc_{st}.txt"
        _check(runner, C.featurecounts_cmd(saf, str(fc), bams, threads,
                                           strandedness=st, paired=paired))
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


def run_pipeline(config, work_dir, refs_root, samplesheet_path, runner=subprocess.run,
                 allow_qc_fail=False, check_files=True):
    work_dir = Path(work_dir)
    out = work_dir / "out" / config.run_name
    for d in ["00_inputs", "01_qc_raw", "02_trimmed", "03_qc_trimmed", "04_align",
              "05_counts", "06_deseq", "qc"]:
        (out / d).mkdir(parents=True, exist_ok=True)

    samples = apply_layout(read_samplesheet(samplesheet_path), config.reads.layout)
    paired = is_paired(samples)
    conds = [s.condition for s in samples]
    cons = expand_contrasts(config.contrasts.explicit, config.contrasts.all_vs_all, conds)
    problems = check_samplesheet(samples, cons, check_files=check_files)
    if problems:
        raise ValueError("sample sheet problems:\n  " + "\n  ".join(problems))
    strand = config.reference.strandedness
    threads = _threads(config)

    # Provenance: exactly what this run was given.
    (out / "00_inputs" / "config.yaml").write_text(
        yaml.safe_dump(config.model_dump(), sort_keys=False))
    shutil.copyfile(samplesheet_path, out / "00_inputs" / "samplesheet.tsv")

    ref = config.reference
    bundle = build_bundle(ref.species, refs_root=refs_root,
                          out_dir=_bundle_dir(config, work_dir), threads=threads,
                          fasta=ref.fasta, gff=ref.gff, feature_types=tuple(ref.feature_types),
                          id_attribute=ref.id_attribute, exclude_seqids=ref.exclude_seqids,
                          seqid_map=ref.seqid_map, structural_rna=ref.structural_rna)

    fastqs = [f for s in samples for f in (s.fastq_r1, s.fastq_r2) if f]
    _check(runner, C.fastqc_cmd(fastqs, str(out / "01_qc_raw"), threads))

    # Trim + align several samples at once, splitting the thread budget.
    par, per_sample = _parallelism(config, threads, len(samples))
    with ThreadPoolExecutor(max_workers=par) as pool:
        futures = {s.sample_id: pool.submit(_process_sample, s, paired, out, bundle,
                                            per_sample, runner) for s in samples}
        reads = {sid: f.result() for sid, f in futures.items()}   # re-raises a failure
    align_logs = {s.sample_id: str(out / "04_align" / f"{s.sample_id}.bowtie2.log")
                  for s in samples}
    bams = [str(out / "04_align" / f"{s.sample_id}.bam") for s in samples]
    trimmed = sorted(str(p) for p in (out / "02_trimmed").glob("*.fq.gz"))
    if trimmed:
        _check(runner, C.fastqc_cmd(trimmed, str(out / "03_qc_trimmed"), threads))

    fc = out / "05_counts" / "featurecounts.txt"
    _check(runner, C.featurecounts_cmd(bundle["saf"], str(fc), bams, threads,
                                       strandedness=strand, paired=paired))
    _reshape_counts(fc, out / "05_counts" / "counts.tsv",
                    keep_ids=set(bundle.get("gene_ids") or []) or None,
                    other_tsv=out / "05_counts" / "ncrna_counts.tsv")
    strand_fracs = _strand_check(runner, bundle["saf"], bams, threads, paired,
                                 out / "05_counts" / "strand_check", strand,
                                 str(fc) + ".summary")
    ncrna = _ncrna_fracs(fc, bundle.get("structural", {}))
    qc = _collect_qc(align_logs, str(fc) + ".summary", strand_fracs, ncrna, strand)
    _check(runner, C.multiqc_cmd(str(out), str(out / "qc" / "multiqc")))

    invariants = {"strandedness": strand, "n_features": bundle["n_features"],
                  "seqids": bundle["seqids"],
                  "structural_rnas_added": bundle.get("extra_ids", [])}
    params = {"strandedness": strand, "paired": paired, "layout": config.reads.layout,
              "threads": threads, "parallel_samples": par,
              "threads_per_sample": per_sample, "allow_qc_fail": allow_qc_fail}
    provenance = {"plugin": run_report.plugin_version(),
                  "inputs": {"config": "00_inputs/config.yaml",
                             "samplesheet": "00_inputs/samplesheet.tsv"},
                  "reference": {"bundle_dir": str(_bundle_dir(config, work_dir)),
                                "reused": bundle.get("reused", False),
                                "fasta_md5": md5(bundle["fasta"]) if Path(bundle["fasta"]).exists() else None,
                                "gff_md5": bundle.get("stamp", {}).get("gff_md5"),
                                "structural_rna_md5": bundle.get("stamp", {}).get("structural_rna_md5"),
                                "saf_md5": bundle.get("saf_md5")},
                  "reads": reads}
    outputs = {"counts": "05_counts/counts.tsv", "ncrna_counts": "05_counts/ncrna_counts.tsv",
               "strand_check": "05_counts/strand_check", "multiqc": "qc/multiqc"}

    # QC gate: a FAIL halts before DESeq2 so bad data cannot silently produce a DE table.
    failed = [sid for sid, v in qc.items() if v.get("verdict") == "FAIL"]
    if failed and not allow_qc_fail:
        rep = run_report.build_report(config.run_name, params, qc, invariants, [],
                                      outputs, "qc_fail", provenance=provenance)
        run_report.write_report(rep, out / "00_run_report.json")
        return rep

    use_batch = 1 if (config.design.batch_variable and any(s.batch for s in samples)) else 0
    coldata = out / "06_deseq" / "coldata.tsv"
    hdr = "sample\tcondition" + ("\tbatch" if use_batch else "") + "\n"
    coldata.write_text(hdr + "".join(
        f"{s.sample_id}\t{s.condition}" + (f"\t{s.batch}" if use_batch else "") + "\n"
        for s in samples))
    ctsv = out / "06_deseq" / "contrasts.tsv"
    ctsv.write_text("".join(f"{c.name}\t{c.numerator}\t{c.denominator}\n" for c in cons))
    _run_deseq2(runner, Path(__file__).parents[1] / "r" / "deseq2.R",
                out / "05_counts" / "counts.tsv", coldata, ctsv, out / "06_deseq", use_batch)

    outputs["deseq"] = "06_deseq/results"
    rep = run_report.build_report(
        config.run_name, params, qc, invariants, [c.name for c in cons], outputs, "ok",
        provenance=provenance,
        de_summary=_de_summary(out / "06_deseq" / "results", cons,
                               config.thresholds.padj, config.thresholds.log2fc))
    run_report.write_report(rep, out / "00_run_report.json")
    return rep
