import subprocess
import sys
import yaml
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def test_validate_ok(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "r", "reference": {"species": "mabs"}}))
    r = subprocess.run([sys.executable, "-m", "engine.python.cli", "validate", str(cfg)],
                       cwd=REPO, capture_output=True)
    assert r.returncode == 0


def test_validate_bad(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "r", "reference": {"species": "custom"}}))
    r = subprocess.run([sys.executable, "-m", "engine.python.cli", "validate", str(cfg)],
                       cwd=REPO, capture_output=True)
    assert r.returncode == 2
