"""Exercise test cache isolation in fresh pytest processes, as CI does."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
CACHE_VARIABLES = ("OSTEOSARC_CACHE", "OPENVAX_DATA_CACHE")


@pytest.mark.parametrize("ambient", [(), CACHE_VARIABLES[:1], CACHE_VARIABLES[1:], CACHE_VARIABLES])
def test_ambient_cache_environment_does_not_redirect_tests(tmp_path, ambient):
    shutil.copyfile(ROOT / "tests/conftest.py", tmp_path / "conftest.py")
    (tmp_path / "test_cache_consumer.py").write_text('''
from hashlib import sha256
import json
import os
from pathlib import Path

import pytest
from vaxrank import download_test_data as downloader


def test_declared_fallback_reads_its_own_cached_object(tmp_path, monkeypatch):
    source = tmp_path / "source.bam"
    source.write_bytes(b"local fixture")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(dict(
        schema_version=1, dataset="fixture", data_version="v1", cases=[],
        assets=[dict(filename=source.name, url=source.as_uri(),
                     sha256=sha256(source.read_bytes()).hexdigest(),
                     size_bytes=source.stat().st_size)])))
    root = tmp_path / "test-cache"
    expected = downloader.download_test_data(manifest_path=manifest, cache_root=root)
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(root))
    assert downloader.download_test_data(manifest_path=manifest, offline=True) == expected


def test_next_test_uses_its_platform_default(tmp_path, monkeypatch):
    root = tmp_path / "default"
    monkeypatch.setattr("datacache.common.appdirs.user_cache_dir", lambda namespace: str(root))
    assert downloader.cache_root_for() == root
    assert not root.exists()


@pytest.mark.shared_read_cache
def test_opt_in_retains_the_callers_cache(tmp_path, monkeypatch):
    ambient = json.loads(os.environ["VAXRANK_TEST_AMBIENT_CACHE"])
    fallback = tmp_path / "default"
    monkeypatch.setattr("datacache.common.appdirs.user_cache_dir", lambda namespace: str(fallback))
    expected = ambient.get("OSTEOSARC_CACHE") or ambient.get("OPENVAX_DATA_CACHE") or fallback
    assert downloader.cache_root_for() == Path(expected)


def test_explicit_test_environment_still_controls_precedence(tmp_path, monkeypatch):
    preferred, fallback, explicit = [tmp_path / name for name in ("preferred", "fallback", "explicit")]
    monkeypatch.setenv("OSTEOSARC_CACHE", str(preferred))
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(fallback))
    assert downloader.cache_root_for() == preferred
    assert downloader.cache_root_for(explicit) == explicit
    monkeypatch.delenv("OSTEOSARC_CACHE")
    assert downloader.cache_root_for() == fallback
''')
    env = os.environ.copy()
    roots = {}
    for name in CACHE_VARIABLES:
        env.pop(name, None)
        if name in ambient:
            path = tmp_path / name.lower()
            path.mkdir()
            roots[name] = str(path)
    env.update(roots)
    temp_root = tmp_path / "pytest-temp"
    temp_root.mkdir()
    env.update(
        VAXRANK_TEST_AMBIENT_CACHE=json.dumps(roots),
        PYTEST_DEBUG_TEMPROOT=str(temp_root),
        PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
        PYTHONPATH=str(ROOT),
    )
    env.pop("PYTEST_ADDOPTS", None)
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "test_cache_consumer.py"],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "4 passed" in result.stdout
    assert all(not list(Path(root).iterdir()) for root in roots.values())
