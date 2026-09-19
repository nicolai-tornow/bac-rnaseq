from __future__ import annotations
import subprocess
from pathlib import Path
import pandas as pd
from .ids import normalize_id


def read_common_vst(paths):
    """Genes present in every input, matched case-insensitively (MAB0001 = mab0001).
    The output keeps the first input's spelling of each ID."""
    dfs = [pd.read_csv(p, sep="\t", index_col=0) for p in paths]
    maps = [{normalize_id(str(g)): g for g in d.index} for d in dfs]
    common = set(maps[0])
    for m in maps[1:]:
        common &= set(m)
    if not common:
        raise ValueError("no genes common to all VST inputs — different references?")
    common = sorted(common)
    parts = [d.loc[[m[k] for k in common]].set_axis([maps[0][k] for k in common])
             for d, m in zip(dfs, maps)]
    return pd.concat(parts, axis=1)


def assert_same_reference(paths, min_overlap=0.95):
    sets = [{normalize_id(str(g)) for g in pd.read_csv(p, sep="\t", index_col=0).index}
            for p in paths]
    base = sets[0]
    for s in sets[1:]:
        overlap = len(base & s) / max(len(base | s), 1)
        if overlap < min_overlap:
            raise ValueError(f"VST inputs share only {overlap:.0%} of genes (<{min_overlap:.0%}) "
                             "— likely aligned to different references; do not merge")


def run_combat(matrix, meta, out_dir, r_script=None, runner=subprocess.run):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mpath, metapath = out_dir / "combined_vst.tsv", out_dir / "meta.tsv"
    matrix.to_csv(mpath, sep="\t")
    meta.to_csv(metapath, sep="\t")
    r_script = r_script or (Path(__file__).parents[1] / "r" / "batch_combat.R")
    res = runner(["Rscript", str(r_script), str(mpath), str(metapath), str(out_dir)],
                 capture_output=True, text=True)
    if getattr(res, "returncode", 0) != 0:
        raise RuntimeError(f"ComBat failed:\n{getattr(res, 'stderr', '')}")
    return {"corrected": str(out_dir / "corrected_vst.tsv"),
            "pca_before": str(out_dir / "pca_before.pdf"),
            "pca_after": str(out_dir / "pca_after.pdf")}
