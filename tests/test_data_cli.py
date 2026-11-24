"""Tests for ``vaxrank data dirs|status|path|prune`` and read-only cache resolution."""

import hashlib
import json

import datacache
import pytest

from vaxrank.cli.discoverability import resolve_command_line


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for key in (
        "VAXRANK_REF_PEPTIDES_DIR",
        "OPENVAX_DATA_CACHE",
        "DATACACHE_OFFLINE",
        "OPENVAX_OFFLINE",
    ):
        monkeypatch.delenv(key, raising=False)


def test_data_dirs_path_status_readonly_and_json(tmp_path, monkeypatch, capsys):
    shared = tmp_path / "shared"
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")

    handled, _ = resolve_command_line(["data", "path"])
    assert handled is True
    assert capsys.readouterr().out.strip() == str(shared / "vaxrank")
    assert not shared.exists()

    handled, _ = resolve_command_line(["data", "dirs", "--json"])
    assert handled is True
    payload = json.loads(capsys.readouterr().out)
    assert payload["cache_dir"] == str(shared / "vaxrank")
    assert payload["origin"] == "shared_env"
    assert payload["detail"] == "OPENVAX_DATA_CACHE"
    assert payload["exists"] is False
    assert payload["offline"] is True
    assert not shared.exists()

    handled, _ = resolve_command_line(["data", "status", "--json"])
    assert handled is True
    status = json.loads(capsys.readouterr().out)
    assert status["kmer_index_count"] == 0
    assert not shared.exists()


def test_data_status_and_prune_dry_run_then_apply(tmp_path, monkeypatch, capsys):
    cache_dir = tmp_path / "vaxrank_cache"
    monkeypatch.setenv("VAXRANK_REF_PEPTIDES_DIR", str(cache_dir))

    # Seed one k-mer index file and one 2-generation datacache bundle store
    cache_dir.mkdir(parents=True)
    idx = cache_dir / "content_demo_kmer_set_8_11.pkl.gz"
    idx.write_bytes(b"kmer-index-bytes" * 100)

    store = cache_dir / "bundles" / "demo"
    src1 = tmp_path / "v1.bin"
    b1 = b"gen1" * 40
    src1.write_bytes(b1)
    datacache.install_bundle(
        store,
        {"a.bin": {"url": src1.as_uri(), "sha256": hashlib.sha256(b1).hexdigest(), "size": len(b1)}},
    )
    src2 = tmp_path / "v2.bin"
    b2 = b"gen2" * 50
    src2.write_bytes(b2)
    datacache.install_bundle(
        store,
        {"a.bin": {"url": src2.as_uri(), "sha256": hashlib.sha256(b2).hexdigest(), "size": len(b2)}},
        force=True,
    )

    resolve_command_line(["data", "status", "--json"])
    status = json.loads(capsys.readouterr().out)
    assert status["kmer_index_count"] == 1
    assert status["total_bytes"] == idx.stat().st_size

    # Dry-run default removes nothing
    resolve_command_line(["data", "prune", "--kmer-indexes", "--json"])
    dry = json.loads(capsys.readouterr().out)
    assert dry["dry_run"] is True
    assert dry["removed_bundle_count"] == 1
    assert len(dry["removed_kmer_indexes"]) == 1
    assert idx.exists()

    # --yes deletes stale bundle generation and k-mer index
    resolve_command_line(["data", "prune", "--kmer-indexes", "--yes", "--json"])
    applied = json.loads(capsys.readouterr().out)
    assert applied["dry_run"] is False
    assert applied["removed_bundle_count"] == 1
    assert not idx.exists()
