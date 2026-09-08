from engine.python.config import Contrast
from engine.python.contrasts import expand_contrasts


def test_explicit_only():
    ex = [Contrast(name="S_vs_7", numerator="SCFM2", denominator="7H9")]
    out = expand_contrasts(ex, False, ["7H9", "SCFM2", "ALI"])
    assert [(c.numerator, c.denominator) for c in out] == [("SCFM2", "7H9")]


def test_all_vs_all_dedups_explicit():
    ex = [Contrast(name="S_vs_7", numerator="SCFM2", denominator="7H9")]
    out = expand_contrasts(ex, True, ["7H9", "SCFM2"])
    pairs = [(c.numerator, c.denominator) for c in out]
    assert ("SCFM2", "7H9") in pairs
    assert pairs.count(("SCFM2", "7H9")) == 1
    assert ("7H9", "SCFM2") in pairs
