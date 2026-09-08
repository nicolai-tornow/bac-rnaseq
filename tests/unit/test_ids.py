from engine.python.ids import normalize_id, match_ids


def test_normalize_case_and_underscore():
    assert normalize_id("mab3648") == "MAB3648"
    assert normalize_id("MAB_0812") == "MAB0812"
    assert normalize_id("MAB0813c") == "MAB0813C"


def test_match_ids_reports_unmatched():
    mapping, frac = match_ids(["mab0001", "mab9999"], ["MAB0001", "MAB0002c"])
    assert mapping["mab0001"] == "MAB0001"
    assert 0.49 < frac < 0.51
