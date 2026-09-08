from __future__ import annotations
import argparse
import os
import shutil
import sys
from .config import load_config
from .build_refs import build_bundle
from .coredetect import suggest_threads
from . import siteconfig


def _validate(args):
    try:
        load_config(args.config)
    except Exception as e:
        print(f"config invalid: {e}", file=sys.stderr)
        return 2
    print("config ok")
    return 0


def _build_refs(args):
    res = build_bundle(args.species, refs_root=args.refs_root, out_dir=args.out,
                       threads=args.threads or suggest_threads(),
                       fasta=args.fasta, gff=args.gff)
    print(f"built {args.species}: {res['n_features']} features on {res['seqids']}")
    return 0


def _doctor(args):
    total = os.cpu_count() or 1
    sug = suggest_threads(total)
    tools = {t: shutil.which(t) for t in
             ["fastp", "bowtie2", "samtools", "featureCounts", "multiqc", "Rscript"]}
    print(f"cores detected: {total}; suggested budget: {sug}")
    for t, p in tools.items():
        print(f"  {t}: {'OK' if p else 'MISSING'}")
    site = siteconfig.read_site()
    site["threads"] = sug
    siteconfig.write_site(site)
    return 0 if all(tools.values()) else 1


def main(argv=None):
    ap = argparse.ArgumentParser(prog="bac-rnaseq")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("config")
    v.set_defaults(fn=_validate)
    b = sub.add_parser("build-refs")
    b.add_argument("--species", required=True)
    b.add_argument("--refs-root", default="refs")
    b.add_argument("--out", required=True)
    b.add_argument("--fasta")
    b.add_argument("--gff")
    b.add_argument("--threads", type=int)
    b.set_defaults(fn=_build_refs)
    d = sub.add_parser("doctor")
    d.set_defaults(fn=_doctor)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
