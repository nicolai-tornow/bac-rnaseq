import pandas as pd
from engine.viz.tpm import compute_tpm


def test_tpm_sums_to_1e6():
    counts = pd.DataFrame({"s1": [100, 200, 0], "s2": [50, 50, 10]},
                          index=["g1", "g2", "MABr5051"])
    lengths = {"g1": 1000, "g2": 2000, "MABr5051": 1500}
    tpm = compute_tpm(counts, lengths, ncrna_ids=["MABr5051"])
    assert "MABr5051" not in tpm.index
    assert abs(tpm["s1"].sum() - 1e6) < 1
