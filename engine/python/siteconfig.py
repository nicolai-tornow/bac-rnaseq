import os
from pathlib import Path
import yaml


def site_path() -> Path:
    base = Path(os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config")))
    return base / "bac-rnaseq" / "site.yaml"


def read_site() -> dict:
    p = site_path()
    return yaml.safe_load(p.read_text()) if p.exists() else {}


def write_site(data: dict) -> Path:
    p = site_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(yaml.safe_dump(data, sort_keys=True))
    return p
