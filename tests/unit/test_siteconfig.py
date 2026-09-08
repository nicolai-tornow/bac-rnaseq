from engine.python import siteconfig


def test_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    siteconfig.write_site({"env_prefix": "/shared/env", "threads": 8})
    got = siteconfig.read_site()
    assert got["env_prefix"] == "/shared/env" and got["threads"] == 8
