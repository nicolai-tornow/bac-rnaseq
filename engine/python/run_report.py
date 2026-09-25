from __future__ import annotations
import json
import os
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


PLUGIN_ROOT = Path(__file__).resolve().parents[2]


def plugin_version() -> dict:
    """Plugin version from its manifest, plus the git commit when installed from a clone."""
    info = {"version": None, "commit": None}
    try:
        info["version"] = json.loads((PLUGIN_ROOT / ".claude-plugin" / "plugin.json")
                                     .read_text()).get("version")
    except Exception:
        pass
    try:
        out = subprocess.run(["git", "-C", str(PLUGIN_ROOT), "rev-parse", "HEAD"],
                             capture_output=True, text=True)
        if out.returncode == 0:
            info["commit"] = out.stdout.strip()
    except Exception:
        pass
    return info


def build_report(run_name, params, samples_qc, invariants, contrasts, outputs, status,
                 provenance=None, de_summary=None, warnings=None, failure=None):
    """status: ok | qc_fail | failed (then `failure` says where and why)."""
    rep = {"schema_version": "1.2", "run_name": run_name,
           "timestamp": datetime.now(timezone.utc).isoformat(),
           "versions": tool_versions(), "provenance": provenance or {},
           "params": params, "samples_qc": samples_qc, "invariants": invariants,
           "contrasts": contrasts, "de_summary": de_summary or {},
           "outputs": outputs, "warnings": list(warnings or []), "status": status}
    if failure is not None:
        rep["failure"] = failure
    return rep


def write_report(report, path):
    """Atomic: a crash mid-write never leaves a truncated report."""
    p = Path(path)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(report, indent=2))
    os.replace(tmp, p)
    return path
