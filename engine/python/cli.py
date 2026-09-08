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
    rep = run_module.run_pipeline(cfg, args.work_dir, args.refs_root, args.samplesheet,
                                  allow_qc_fail=args.allow_qc_fail)
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


def _enrich(args):
    from ..viz.volcano import load_results
    from ..enrichment.genesets import load_genesets, align_hits
    from ..enrichment.ora import run_ora
    from ..enrichment.bubble import bubble
    df = load_results(args.results)
    sets = load_genesets(args.categories, args.sheet, args.id_col, args.cat_col)
    geneset_ids = set().union(*sets.values())
    mapping = align_hits(list(df.index), geneset_ids, max_unmatched=args.max_unmatched)
    universe = {mapping[g] for g in df.index if g in mapping}
    hits = {mapping[g] for g in df.index
            if g in mapping and df.loc[g, "padj"] < args.padj and abs(df.loc[g, "log2FoldChange"]) >= args.log2fc}
    res = run_ora(hits, sets, universe, alternative=args.alternative, min_size=args.min_size)
    res.to_csv(args.out + ".tsv", sep="\t", index=False)
    pdf, png = bubble(res, args.out)
    print(f"wrote {args.out}.tsv, {pdf}, {png}")
    return 0


def _batch(args):
    import pandas as pd
    from ..python.batch import assert_same_reference, read_common_vst, run_combat
    assert_same_reference(args.vst)
    mat = read_common_vst(args.vst)
    meta = pd.read_csv(args.meta, sep="\t", index_col=0)
    res = run_combat(mat, meta, args.out_dir)
    print(f"corrected: {res['corrected']}; PCA: {res['pca_before']}, {res['pca_after']}")
    return 0


def _export(args):
    from ..python.gene_names import gene_names_from_gff
    from ..python.export_xlsx import export_workbook
    gn = gene_names_from_gff(args.gff) if args.gff else {}
    export_workbook(args.results_dir, args.out, gene_names=gn,
                    normalized_tsv=args.normalized, vst_tsv=args.vst)
    print(f"wrote {args.out}")
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
    rn.add_argument("--allow-qc-fail", dest="allow_qc_fail", action="store_true",
                    help="proceed to DESeq2 even if a sample FAILs QC (recorded in the run report)")
    rn.set_defaults(fn=_run)
    vz = sub.add_parser("visualize")
    vz.add_argument("results")
    vz.add_argument("--out", required=True)
    vz.add_argument("--top", type=int, default=10)
    vz.add_argument("--genes")
    vz.set_defaults(fn=_visualize)
    en = sub.add_parser("enrich")
    en.add_argument("results")
    en.add_argument("--categories", required=True)
    en.add_argument("--out", required=True)
    en.add_argument("--sheet", default="2 Gene Annotations")
    en.add_argument("--id-col", dest="id_col", default="Gene ID")
    en.add_argument("--cat-col", dest="cat_col", default="Module")
    en.add_argument("--padj", type=float, default=0.05)
    en.add_argument("--log2fc", type=float, default=1.0)
    en.add_argument("--alternative", default="two-sided")
    en.add_argument("--min-size", dest="min_size", type=int, default=2)
    en.add_argument("--max-unmatched", dest="max_unmatched", type=float, default=0.05)
    en.set_defaults(fn=_enrich)
    bt = sub.add_parser("batch")
    bt.add_argument("--vst", nargs="+", required=True)
    bt.add_argument("--meta", required=True)
    bt.add_argument("--out-dir", dest="out_dir", required=True)
    bt.set_defaults(fn=_batch)
    ex = sub.add_parser("export")
    ex.add_argument("--results-dir", dest="results_dir", required=True)
    ex.add_argument("--out", required=True)
    ex.add_argument("--gff")
    ex.add_argument("--normalized")
    ex.add_argument("--vst")
    ex.set_defaults(fn=_export)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
