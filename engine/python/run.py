from __future__ import annotations
import subprocess
from pathlib import Path
from .samplesheet import read_samplesheet, is_paired
from .contrasts import expand_contrasts
from .build_refs import build_bundle
from . import commands as C
from . import qc_triage, run_report


def _check(runner, cmd, **kw):
    res = runner(cmd, capture_output=True, text=True, **kw)
    if getattr(res, "returncode", 0) != 0:
        raise RuntimeError(f"stage failed: {cmd[0]}\n{getattr(res, 'stderr', '')}")
    return res


def _run_deseq2(runner, r_script, counts, coldata, contrasts_tsv, out_dir, use_batch):
    _check(runner, ["Rscript", str(r_script), str(counts), str(coldata),
                    str(contrasts_tsv), str(out_dir), str(use_batch)])


def _collect_qc(align_logs, fc_summary, strandedness):
    out = {}
    fc = qc_triage.parse_featurecounts_summary(fc_summary)
    frac = fc["assigned_frac"]
    for sid, log in align_logs.items():
        pct = qc_triage.parse_bowtie2_log(Path(log).read_text())
        out[sid] = qc_triage.triage_sample(pct, frac, 0.0, strandedness=strandedness)
    return out


def _reshape_counts(fc_txt, out_tsv):
    lines = [l for l in Path(fc_txt).read_text().splitlines() if not l.startswith("#")]
    header = lines[0].split("\t")
    keep = list(range(6, len(header)))  # sample columns follow Geneid,Chr,Start,End,Strand,Length
    samp = [Path(header[i]).name[:-4] if header[i].endswith(".bam") else header[i] for i in keep]
    with open(out_tsv, "w") as fh:
        fh.write("Gene\t" + "\t".join(samp) + "\n")
        for l in lines[1:]:
            f = l.split("\t")
            fh.write(f[0] + "\t" + "\t".join(f[i] for i in keep) + "\n")
    return out_tsv


def run_pipeline(config, work_dir, refs_root, samplesheet_path, runner=subprocess.run,
                 allow_qc_fail=False):
    work_dir = Path(work_dir)
    out = work_dir / "out" / config.run_name
    for d in ["01_qc_raw", "02_trimmed", "03_qc_trimmed", "04_align", "05_counts", "06_deseq"]:
        (out / d).mkdir(parents=True, exist_ok=True)

    samples = read_samplesheet(samplesheet_path)
    paired = is_paired(samples)
    strand = config.reference.strandedness
    threads = config.resources.threads or 4

    bundle = build_bundle(config.reference.species, refs_root=refs_root,
                          out_dir=out.parent / "refs" / config.reference.species,
                          threads=threads, fasta=config.reference.fasta,
                          gff=config.reference.gff)

    align_logs, bams = {}, []
    for s in samples:
        t1 = out / "02_trimmed" / f"{s.sample_id}_R1.fq.gz"
        t2 = out / "02_trimmed" / f"{s.sample_id}_R2.fq.gz" if paired else None
        _check(runner, C.fastp_cmd(s.fastq_r1, str(t1), threads,
                                   r2=s.fastq_r2, out2=str(t2) if t2 else None,
                                   json=str(out / "02_trimmed" / f"{s.sample_id}.json"),
                                   html=str(out / "02_trimmed" / f"{s.sample_id}.html")))
        bam = out / "04_align" / f"{s.sample_id}.bam"
        sam = out / "04_align" / f"{s.sample_id}.sam"
        log = out / "04_align" / f"{s.sample_id}.bowtie2.log"
        bt2 = C.bowtie2_cmd(bundle["index_prefix"], str(t1), threads, r2=str(t2) if t2 else None)
        res = runner(bt2 + ["-S", str(sam)], capture_output=True, text=True)
        Path(log).write_text(getattr(res, "stderr", "") or "")
        _check(runner, ["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)])
        _check(runner, ["samtools", "index", str(bam)])
        align_logs[s.sample_id] = str(log)
        bams.append(str(bam))

    fc = out / "05_counts" / "featurecounts.txt"
    _check(runner, C.featurecounts_cmd(bundle["saf"], str(fc), bams, threads,
                                       strandedness=strand, paired=paired))
    _reshape_counts(fc, out / "05_counts" / "counts.tsv")

    qc = _collect_qc(align_logs, str(fc) + ".summary", strand)
    invariants = {"strandedness": strand, "n_features": bundle["n_features"],
                  "seqids": bundle["seqids"]}
    params = {"strandedness": strand, "paired": paired, "threads": threads,
              "allow_qc_fail": allow_qc_fail}

    # QC gate: a FAIL halts before DESeq2 so bad data cannot silently produce a DE table.
    failed = [sid for sid, v in qc.items() if v.get("verdict") == "FAIL"]
    if failed and not allow_qc_fail:
        rep = run_report.build_report(config.run_name, params, qc, invariants, [],
                                      {"counts": "05_counts/counts.tsv"}, "qc_fail")
        run_report.write_report(rep, out / "00_run_report.json")
        return rep

    conds = [s.condition for s in samples]
    use_batch = 1 if (config.design.batch_variable and any(s.batch for s in samples)) else 0
    coldata = out / "06_deseq" / "coldata.tsv"
    hdr = "sample\tcondition" + ("\tbatch" if use_batch else "") + "\n"
    coldata.write_text(hdr + "".join(
        f"{s.sample_id}\t{s.condition}" + (f"\t{s.batch}" if use_batch else "") + "\n"
        for s in samples))
    cons = expand_contrasts(config.contrasts.explicit, config.contrasts.all_vs_all, conds)
    ctsv = out / "06_deseq" / "contrasts.tsv"
    ctsv.write_text("".join(f"{c.name}\t{c.numerator}\t{c.denominator}\n" for c in cons))
    _run_deseq2(runner, Path(__file__).parents[1] / "r" / "deseq2.R",
                out / "05_counts" / "counts.tsv", coldata, ctsv, out / "06_deseq", use_batch)

    rep = run_report.build_report(
        config.run_name, params, qc, invariants, [c.name for c in cons],
        {"counts": "05_counts/counts.tsv", "deseq": "06_deseq/results"}, "ok")
    run_report.write_report(rep, out / "00_run_report.json")
    return rep
