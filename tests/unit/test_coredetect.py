from engine.python.coredetect import suggest_threads


def test_leaves_headroom():
    assert suggest_threads(16) == 14
    assert suggest_threads(4) == 3
    assert suggest_threads(1) == 1
    assert suggest_threads(2) == 1
