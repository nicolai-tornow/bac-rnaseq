from __future__ import annotations
import argparse
import os
import shutil
import sys
from .config import load_config
from .build_refs import build_bundle
from .coredetect import suggest_threads
from . import siteconfig
from . import run as run_module


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


def _run(args):
    cfg = load_config(args.config)
    rep = run_module.run_pipeline(cfg, args.work_dir, args.refs_root, args.samplesheet)
    print(f"run '{rep.get('run_name')}': {rep.get('status')}")
    return 0 if rep.get("status") == "ok" else 1


def _visualize(args):
    from ..viz.volcano import load_results, volcano
    df = load_results(args.results)
    sel = ({"mode": "genes", "genes": args.genes.split(",")} if args.genes
           else {"mode": "top", "n": args.top, "by": "l2fc"})
    pdf, png = volcano(df, args.out, selection=sel)
    print(f"wrote {pdf}, {png}")
    return 0


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
    rn = sub.add_parser("run")
    rn.add_argument("config")
    rn.add_argument("--work-dir", required=True)
    rn.add_argument("--samplesheet", required=True)
    rn.add_argument("--refs-root", default="refs")
    rn.set_defaults(fn=_run)
    vz = sub.add_parser("visualize")
    vz.add_argument("results")
    vz.add_argument("--out", required=True)
    vz.add_argument("--top", type=int, default=10)
    vz.add_argument("--genes")
    vz.set_defaults(fn=_visualize)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
