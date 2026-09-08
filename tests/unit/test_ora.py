from engine.enrichment.ora import run_ora


def test_ora_detects_enrichment():
    universe = {f"g{i}" for i in range(100)}
    hits = {"g0", "g1", "g2", "g3"}
    gene_sets = {"enriched": {"g0", "g1", "g2", "g3", "g4"},
                 "random": {f"g{i}" for i in range(50, 60)}}
    res = run_ora(hits, gene_sets, universe, min_size=2)
    row = res.set_index("set").loc["enriched"]
    assert row["n_hits"] == 4 and row["odds_ratio"] > 1 and row["p"] < 0.05
    assert "fdr" in res.columns
