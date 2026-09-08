import glob
import os
import pytest
from pathlib import Path
from engine.viz.volcano import load_results
from engine.enrichment.genesets import load_genesets, align_hits
from engine.enrichment.ora import run_ora
from engine.enrichment.bubble import bubble

SHEET = Path("/path/to/projects/figures/supplemental/tables/SuppTable_gene_categorization.xlsx")
CANDS = sorted(glob.glob("/path/to/projects/rnaseq_mabs_media/06_deseq/results/*.tsv"))
pytestmark = pytest.mark.skipif(not (SHEET.exists() and CANDS), reason="need sheet + results")


def test_enrichment_end_to_end(tmp_path):
    df = load_results(CANDS[0])
    sets = load_genesets(str(SHEET), "2 Gene Annotations", "Gene ID", "Module")
    geneset_ids = set().union(*sets.values())
    mapping = align_hits(list(df.index), geneset_ids, max_unmatched=0.10)
    universe = {mapping[g] for g in df.index if g in mapping}
    hits = {mapping[g] for g in df.index
            if g in mapping and df.loc[g, "padj"] < 0.05 and abs(df.loc[g, "log2FoldChange"]) >= 1}
    res = run_ora(hits, sets, universe, alternative="two-sided", min_size=2)
    assert not res.empty
    assert {"set", "p", "fdr", "odds_ratio", "effect"}.issubset(res.columns)
    pdf, png = bubble(res, str(tmp_path / "bub"))
    assert os.path.exists(pdf) and os.path.exists(png)
