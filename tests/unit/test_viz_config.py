import os
import matplotlib
from engine.viz import config as V


def test_style_and_save(tmp_path):
    V.apply_style()
    assert matplotlib.rcParams["pdf.fonttype"] == 42
    import matplotlib.pyplot as plt
    fig = plt.figure()
    plt.plot([0, 1], [0, 1])
    pdf, png = V.save(fig, str(tmp_path / "x"))
    assert pdf.endswith(".pdf") and png.endswith(".png")
    assert os.path.exists(pdf) and os.path.exists(png)
