"""Remove a finished run's regenerable intermediates: never while a run works on the
folder, never a raw FASTQ of the sample sheet, and never anything outside the folder.
`plan` is the dry run; `execute` deletes what a fresh plan lists, under the run lock."""
from __future__ import annotations
import csv
import json
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from . import run_report, runlock

# tier -> (folder, patterns, recursive)
TIERS = {"trimmed": ("02_trimmed", ("*.fq.gz",), False),
         "sam": ("04_align", ("*.sam",), False),
         "fc_temp": ("05_counts", ("temp-core-*", "temp-sort-*"), True),
         "bam": ("04_align", ("*.bam", "*.bai"), False)}
DEFAULT_TIERS = ("trimmed", "sam", "fc_temp")
IDLE_MINUTES = 15.0
MANIFEST = "cleanup_manifest.tsv"


@dataclass
class Candidate:
    path: Path
    bytes: int
    tier: str


@dataclass
class Plan:
    run_dir: Path
    include_bams: bool
    candidates: list = field(default_factory=list)
    selected: list = field(default_factory=list)
    kept_files: int = 0
    kept_bytes: int = 0
    refusals: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    def by_tier(self) -> dict:
        out = {}
        for c in self.candidates:
            n, b = out.get(c.tier, (0, 0))
            out[c.tier] = (n + 1, b + c.bytes)
        return out

    @property
    def tiers(self) -> tuple:
        return DEFAULT_TIERS + (("bam",) if self.include_bams else ())


class CleanupRefused(RuntimeError):
    def __init__(self, reasons):
        super().__init__("; ".join(reasons))
        self.reasons = list(reasons)


def _candidates(run_dir) -> list[Candidate]:
    found = {}
    for tier, (sub, patterns, recursive) in TIERS.items():
        d = run_dir / sub
        if not d.is_dir():
            continue
        for pat in patterns:
            for p in (d.rglob(pat) if recursive else d.glob(pat)):
                if p in found or (p.is_dir() and not p.is_symlink()):
                    continue
                found[p] = Candidate(p, p.lstat().st_size, tier)
    return sorted(found.values(), key=lambda c: str(c.path))


def _sheet_paths(run_dir) -> list[str]:
    """FASTQ paths of 00_inputs/samplesheet.tsv. `run` copies a CSV sheet verbatim
    under that name, so the delimiter is read from the header, not the file name."""
    text = (run_dir / "00_inputs" / "samplesheet.tsv").read_text()
    lines = text.splitlines()
    delim = "\t" if lines and "\t" in lines[0] else ","
    rows = csv.DictReader(lines, delimiter=delim)
    if not rows.fieldnames or "fastq_r1" not in rows.fieldnames:
        raise ValueError("no fastq_r1 column")
    return [v.strip() for r in rows for k in ("fastq_r1", "fastq_r2")
            if (v := r.get(k) or "").strip()]


def _protected(run_dir, paths, cwd) -> tuple[set, set]:
    """Realpaths and (device, inode) of the sheet's FASTQs. A relative path is taken
    relative to every place `run` may have been started from."""
    bases = [run_dir.parent.parent, Path(cwd or os.getcwd()), run_dir / "00_inputs"]
    real, inodes = set(), set()
    for v in paths:
        for c in ([Path(v)] if os.path.isabs(v) else [b / v for b in bases]):
            real.add(os.path.realpath(c))
            try:
                st = os.stat(c)
                inodes.add((st.st_dev, st.st_ino))
            except OSError:
                pass
    return real, inodes


def _newest(run_dir) -> tuple[float | None, str | None]:
    newest, which = None, None
    for root, _, files in os.walk(run_dir):
        for f in files:
            if f.startswith(runlock.LOCK_NAME):
                continue            # the lock is judged by its heartbeat instead
            p = os.path.join(root, f)
            try:
                m = os.lstat(p).st_mtime
            except OSError:
                continue
            if newest is None or m > newest:
                newest, which = m, os.path.relpath(p, run_dir)
    return newest, which


def _all_files(run_dir):
    for root, _, files in os.walk(run_dir):
        for f in files:
            yield Path(root) / f


def plan(run_dir, include_bams=False, idle_minutes=IDLE_MINUTES, now=None, cwd=None,
         own_lock=None) -> Plan:
    """The dry run. Refusals block the whole clean-up."""
    run_dir = Path(run_dir).absolute()
    p = Plan(run_dir, include_bams)
    if not run_dir.is_dir():
        p.refusals.append(f"{run_dir} is not a folder")
        return p
    root = os.path.realpath(run_dir)

    try:
        status = json.loads((run_dir / "00_run_report.json").read_text()).get("status")
        if status != "ok":
            p.refusals.append(f"the run report says status '{status}', not 'ok': finish or "
                              "re-run the analysis first")
    except FileNotFoundError:
        p.refusals.append("there is no 00_run_report.json: the run never finished")
    except (OSError, ValueError):
        p.refusals.append("00_run_report.json cannot be read")

    st, info = runlock.lock_status(run_dir)
    own = own_lock or {}
    if st == "live" and not all(own.get(k) == (info or {}).get(k) for k in ("host", "pid", "started")):
        p.refusals.append(f"a run is using this folder ({runlock.describe(info)})")

    partial = sorted(str(q.relative_to(run_dir)) for q in run_dir.rglob("*.partial"))
    if partial:
        p.refusals.append(f"staging folders are present ({', '.join(partial[:3])}): a run is "
                          "working here or was killed; re-run the analysis first")

    now = time.time() if now is None else now
    newest, which = _newest(run_dir)
    if idle_minutes > 0 and newest is not None and now - newest < idle_minutes * 60:
        retry = datetime.fromtimestamp(newest + idle_minutes * 60).strftime("%H:%M")
        p.refusals.append(f"{which} was modified {(now - newest) / 60:.0f} min ago, so a run "
                          f"may still be writing here. Retry after {retry}, or pass "
                          "--idle-minutes 0 if you are sure nothing is running")

    p.candidates = _candidates(run_dir)
    p.selected = [c for c in p.candidates if c.tier in p.tiers]

    try:
        real, inodes = _protected(run_dir, _sheet_paths(run_dir), cwd)
    except FileNotFoundError:
        real, inodes = set(), set()
        p.refusals.append("00_inputs/samplesheet.tsv is missing, so the raw FASTQs cannot be "
                          "protected")
    except (OSError, ValueError) as e:
        real, inodes = set(), set()
        p.refusals.append(f"00_inputs/samplesheet.tsv cannot be read ({e}), so the raw "
                          "FASTQs cannot be protected")
    for c in p.selected:
        rel = c.path.relative_to(run_dir)
        if c.path.is_symlink():
            p.refusals.append(f"{rel} is a symlink")
            continue
        rp = os.path.realpath(c.path)
        if not rp.startswith(root + os.sep):
            p.refusals.append(f"{rel} resolves outside the run folder ({rp})")
            continue
        try:
            s = os.stat(c.path)
        except OSError:
            continue
        if rp in real or (s.st_dev, s.st_ino) in inodes:
            p.refusals.append(f"{rel} is a raw FASTQ listed in the sample sheet")

    cand = {c.path for c in p.candidates}
    for f in _all_files(run_dir):
        if f not in cand:
            p.kept_files += 1
            p.kept_bytes += f.lstat().st_size

    if include_bams and any(c.tier == "bam" for c in p.selected):
        p.notes.append("the next run re-trims and re-aligns every sample whose BAM is removed")
    if not list((run_dir / "04_align").glob("*.done.json")):
        p.notes.append("this run has no completion markers (made before bac-rnaseq 0.3.0): "
                       "the next run re-trims and re-aligns every sample regardless")
    return p


def _gb(b) -> str:
    return f"{b / 1e9:.2f}"


def format_plan(p: Plan) -> str:
    t = p.by_tier()
    lines = [f"clean-up of {p.run_dir}", "",
             f"{'tier':<9} {'action':<7} {'files':>6} {'GB':>8}"]
    for tier in TIERS:
        n, b = t.get(tier, (0, 0))
        action = "delete" if tier in p.tiers else "keep"
        hint = "   (--include-bams to delete)" if tier == "bam" and n and not p.include_bams else ""
        lines.append(f"{tier:<9} {action:<7} {n:>6} {_gb(b):>8}{hint}")
    lines.append(f"{'other':<9} {'keep':<7} {p.kept_files:>6} {_gb(p.kept_bytes):>8}")
    lines.append(f"{'total':<9} {'delete':<7} {len(p.selected):>6} "
                 f"{_gb(sum(c.bytes for c in p.selected)):>8}")
    lines += [f"note: {n}" for n in p.notes]
    if p.refusals:
        lines += ["", "REFUSED, nothing will be deleted:"] + [f"  - {r}" for r in p.refusals]
    return "\n".join(lines)


def _append_report(run_dir, entry):
    path = run_dir / "00_run_report.json"
    rep = json.loads(path.read_text())
    rep.setdefault("cleanup", []).append(entry)
    run_report.write_report(rep, path)


def execute(run_dir, include_bams=False, idle_minutes=IDLE_MINUTES, cwd=None) -> dict:
    """Delete what a fresh plan lists, holding the run lock so no run can start
    meanwhile. Every deleted file gets a manifest row; the run report gets an entry."""
    run_dir = Path(run_dir).absolute()
    lock = runlock.RunLock(run_dir)
    try:
        lock.acquire()
    except runlock.RunLocked as e:
        raise CleanupRefused([str(e)]) from e
    try:
        p = plan(run_dir, include_bams, idle_minutes, cwd=cwd, own_lock=lock.info)
        if p.refusals:
            raise CleanupRefused(p.refusals)
        plugin = run_report.plugin_version()
        stamp = datetime.now(timezone.utc).isoformat()
        manifest = run_dir / MANIFEST
        new = not manifest.exists()
        deleted = freed = 0
        error = None
        with open(manifest, "a") as fh:
            if new:
                fh.write("path\tbytes\ttier\ttimestamp\tplugin_commit\n")
            for c in p.selected:
                rel = c.path.relative_to(run_dir)
                try:
                    c.path.unlink()
                except FileNotFoundError:
                    continue
                except OSError as e:
                    error = f"{rel}: {e}"
                    break
                fh.write(f"{rel}\t{c.bytes}\t{c.tier}\t{stamp}\t{plugin.get('commit') or ''}\n")
                fh.flush()
                deleted += 1
                freed += c.bytes
        entry = {"timestamp": stamp, "tiers": list(p.tiers), "files": deleted, "bytes": freed,
                 "manifest": MANIFEST, "plugin": plugin}
        if error:
            entry["error"] = error
        _append_report(run_dir, entry)
        if error:
            raise OSError(f"clean-up stopped at {error}; the {deleted} files deleted before "
                          f"are listed in {MANIFEST}")
        return entry
    finally:
        lock.release()
