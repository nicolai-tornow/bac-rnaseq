"""Expected md5 of the raw FASTQs, from sample-sheet columns (md5_r1, md5_r2) or an
md5sum file as sequencing facilities deliver it (`<md5>  <path>`)."""
from __future__ import annotations
import os
import re
from pathlib import Path
from .build_refs import md5

MD5_RE = re.compile(r"^[0-9a-fA-F]{32}$")


def read_md5_manifest(path) -> dict:
    """Keys: each path as written, resolved against the file's folder, and its file
    name. A key listed twice with different sums maps to None (ambiguous)."""
    base = Path(path).resolve().parent
    man: dict = {}

    def put(k, v):
        if k in man and man[k] != v:
            man[k] = None
        else:
            man.setdefault(k, v)

    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2 or not MD5_RE.match(parts[0]):
            raise ValueError(f"{path}: not an md5sum line: {line!r}")
        digest, name = parts[0].lower(), parts[1].strip().lstrip("*")
        put(name, digest)
        put(os.path.normpath(name if os.path.isabs(name) else base / name), digest)
        put(os.path.basename(name), digest)
    return man


def lookup(man, fastq) -> str | None:
    for k in (fastq, os.path.abspath(fastq), os.path.realpath(fastq), os.path.basename(fastq)):
        if k in man:
            if man[k] is None:
                raise ValueError(f"{os.path.basename(fastq)} is listed more than once with "
                                 "different md5 sums; list full paths in the md5 file")
            return man[k]
    return None


def expected_md5s(samples, manifest=None) -> dict[str, str]:
    """FASTQ path -> expected md5. Sample-sheet columns win over the md5 file."""
    out = {}
    for s in samples:
        for fq, col in ((s.fastq_r1, s.md5_r1), (s.fastq_r2, s.md5_r2)):
            if not fq:
                continue
            want = col or (lookup(manifest, fq) if manifest else None)
            if want:
                out[fq] = want.lower()
    return out


def check_md5s(expected: dict[str, str]) -> list[str]:
    problems = []
    for fq, want in expected.items():
        try:
            got = md5(fq)
        except OSError as e:
            problems.append(f"{fq}: cannot be read ({e})")
            continue
        if got != want:
            problems.append(f"{fq}: md5 {got} does not match the expected {want}")
    return problems
