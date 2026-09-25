from __future__ import annotations

_STRAND = {"reverse": "2", "forward": "1", "unstranded": "0"}
# Settings that change the result (the resume fingerprint compares these).
FASTP_SETTINGS = ["--detect_adapter_for_pe", "--qualified_quality_phred", "20",
                  "--length_required", "36"]
BOWTIE2_SETTINGS = ["--sensitive", "--no-unal"]


def strand_flag(strandedness: str) -> str:
    return _STRAND[strandedness]


def fastp_cmd(r1, out1, threads, r2=None, out2=None, json=None, html=None):
    cmd = ["fastp", "--in1", r1, "--out1", out1, "--thread", str(threads), *FASTP_SETTINGS]
    if r2:
        cmd += ["--in2", r2, "--out2", out2]
    if json:
        cmd += ["--json", json]
    if html:
        cmd += ["--html", html]
    return cmd


def fastqc_cmd(files, out_dir, threads):
    return ["fastqc", "-t", str(threads), "-o", out_dir, *files]


def multiqc_cmd(in_dir, out_dir, config=None, replace_names=None):
    cmd = ["multiqc", "-f", "-o", out_dir]
    if config:
        cmd += ["-c", config]
    if replace_names:
        cmd += ["--replace-names", replace_names]
    return cmd + [in_dir]


def bowtie2_cmd(index_prefix, r1, threads, r2=None):
    cmd = ["bowtie2", "-x", index_prefix, *BOWTIE2_SETTINGS, "--threads", str(threads)]
    if r2:
        cmd += ["-1", r1, "-2", r2]
    else:
        cmd += ["-U", r1]
    return cmd


def featurecounts_cmd(saf, out, bams, threads, strandedness, paired, tmp_dir=None):
    cmd = ["featureCounts", "-a", saf, "-F", "SAF", "-o", out,
           "-T", str(threads), "-s", strand_flag(strandedness)]
    if tmp_dir:
        cmd += ["--tmpDir", tmp_dir]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += list(bams)
    return cmd
