from __future__ import annotations
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def _v(tool, args):
    exe = shutil.which(tool)
    if not exe:
        return None
    try:
        out = subprocess.run([exe, *args], capture_output=True, text=True)
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception:
        return None


def tool_versions() -> dict:
    return {"fastp": _v("fastp", ["--version"]),
            "bowtie2": _v("bowtie2", ["--version"]),
            "samtools": _v("samtools", ["--version"]),
            "featureCounts": _v("featureCounts", ["-v"])}


def build_report(run_name, params, samples_qc, invariants, contrasts, outputs, status):
    return {"schema_version": "1.0", "run_name": run_name,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "versions": tool_versions(), "params": params,
            "samples_qc": samples_qc, "invariants": invariants,
            "contrasts": contrasts, "outputs": outputs, "status": status}


def write_report(report, path):
    Path(path).write_text(json.dumps(report, indent=2))
    return path
