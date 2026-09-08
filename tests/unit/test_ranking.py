import pandas as pd
from engine.viz.ranking import top_bottom


def test_top_bottom():
    tpm = pd.DataFrame({"s1": [100, 1, 50, 5], "s2": [120, 2, 40, 3]},
                       index=["hi", "lo", "mid", "lo2"])
    top, bottom = top_bottom(tpm, 1)
    assert top == ["hi"]
    assert bottom == ["lo"]
