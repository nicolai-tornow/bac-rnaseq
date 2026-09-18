import json
from engine.python.run_report import build_report, write_report, plugin_version


def test_report_shape(tmp_path):
    rep = build_report("run1", {"strandedness": "reverse"},
                       {"s1": {"verdict": "PASS", "reasons": []}},
                       {"strandedness_ok": True, "n_features": 4970},
                       ["SCFM2_vs_7H9"], {"counts": "05_counts/counts.tsv"}, "ok")
    assert rep["schema_version"] == "1.1"
    assert rep["params"]["strandedness"] == "reverse"
    assert rep["invariants"]["n_features"] == 4970
    p = tmp_path / "00_run_report.json"
    write_report(rep, p)
    assert json.loads(p.read_text())["run_name"] == "run1"


def test_plugin_version_reads_manifest():
    v = plugin_version()
    assert v["version"]            # from .claude-plugin/plugin.json
