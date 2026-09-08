from engine.python.config import load_config


def test_pipeline_config_defaults():
    cfg = load_config({"run_name": "r", "reference": {"species": "mabs"},
                       "contrasts": {"explicit": [{"name": "S_vs_7", "numerator": "SCFM2", "denominator": "7H9"}]}})
    assert cfg.design.variable == "condition"
    assert cfg.design.batch_variable is None
    assert cfg.thresholds.padj == 0.05 and cfg.thresholds.log2fc == 1.0
    assert cfg.contrasts.explicit[0].numerator == "SCFM2"
    assert cfg.contrasts.all_vs_all is False
