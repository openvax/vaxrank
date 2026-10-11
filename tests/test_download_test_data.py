"""Offline downloader integrity, cache sharing and original-read contracts."""

from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor
import errno
import json
import shutil
import stat
from pathlib import Path
from types import SimpleNamespace

from datacache import Cache, FileValidationError
from isovar import ReadCollector
import pytest

from vaxrank import download_test_data as downloader
from vaxrank.sid_test_data import sid_test_data, sid_reads


ROOT = Path(__file__).resolve().parents[1]
DATA = sid_test_data() / "osteosarc/shared-v1"
MANIFEST = downloader.load_manifest()


@pytest.fixture
def tiny_manifest(tmp_path):
    upstream = tmp_path / "upstream"
    upstream.mkdir()
    assets = []
    for name, content in (("first.bam", b"source one"), ("second.bam", b"source two")):
        source = upstream / name
        source.write_bytes(content)
        assets.append(dict(filename=name, url=source.as_uri(), sha256=sha256(content).hexdigest(), size_bytes=len(content)))
    manifest = dict(schema_version=1, dataset="fixture", data_version="v1", assets=assets, cases=[])
    path = tmp_path / "source-manifest.json"
    path.write_text(json.dumps(manifest))
    return path, manifest


def forbid_fetch(*args, **kwargs):
    raise AssertionError("Unexpected network/cache mutation")


def test_shared_cache_and_offline_generation(tmp_path, tiny_manifest, monkeypatch):
    path, manifest = tiny_manifest
    root = tmp_path / "shared"
    paths = downloader.download_test_data(manifest_path=path, cache_root=root)
    for asset in manifest["assets"]:
        expected = root / "objects/sha256" / downloader.cache_filename(asset)
        assert paths[asset["filename"]] == expected
    monkeypatch.setattr(Cache, "fetch", forbid_fetch)
    # A second application/process has only the same root and pinned manifest.
    monkeypatch.setenv(downloader.CACHE_ENVIRONMENT, str(root))
    output = downloader.download_test_data(tmp_path / "export", manifest_path=path, offline=True)
    assert downloader.verify_dataset(output, manifest) == output
    # Existing valid exports need neither cache nor network, even on read-only use.
    monkeypatch.setattr(downloader, "Cache", forbid_fetch)
    assert downloader.download_test_data(output, manifest_path=path) == output


def test_missing_offline_cache_has_no_side_effects(tmp_path, tiny_manifest, monkeypatch):
    path, _ = tiny_manifest
    monkeypatch.setattr(Cache, "fetch", forbid_fetch)
    with pytest.raises(FileNotFoundError):
        downloader.download_test_data(tmp_path / "export", manifest_path=path, cache_root=tmp_path / "absent", offline=True)
    assert not (tmp_path / "absent").exists()
    assert not (tmp_path / "export").exists()


def test_cache_corruption_requires_explicit_repair(tmp_path, tiny_manifest):
    path, manifest = tiny_manifest
    root = tmp_path / "cache"
    paths = downloader.download_test_data(manifest_path=path, cache_root=root)
    paths["first.bam"].write_bytes(b"corrupt")
    with pytest.raises(FileValidationError):
        downloader.download_test_data(manifest_path=path, cache_root=root)
    with pytest.raises(FileValidationError):
        downloader.download_test_data(manifest_path=path, cache_root=root, offline=True)
    downloader.download_test_data(manifest_path=path, cache_root=root, repair_cache=True)
    assert sha256(paths["first.bam"].read_bytes()).hexdigest() == manifest["assets"][0]["sha256"]


def test_concurrent_exports_converge_without_mixed_assets(tmp_path, tiny_manifest):
    path, manifest = tiny_manifest
    output, root = tmp_path / "export", tmp_path / "cache"
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(downloader.download_test_data, output,
                    manifest_path=path, cache_root=root) for _ in range(2)]
        assert [f.result() for f in futures] == [output, output]
    assert downloader.verify_dataset(output, manifest) == output


def test_concurrent_empty_destination_is_preserved(tmp_path, tiny_manifest, monkeypatch):
    path, _ = tiny_manifest
    output = tmp_path / "export"
    original_verify = downloader.verify_dataset
    created = []

    def create_foreign_destination_after_staging(directory, manifest):
        result = original_verify(directory, manifest)
        if directory != output:
            output.mkdir(mode=0o700)
            created.append(output.stat())
        return result

    monkeypatch.setattr(downloader, "verify_dataset", create_foreign_destination_after_staging)
    with pytest.raises(ValueError):
        downloader.download_test_data(output, manifest_path=path, cache_root=tmp_path / "cache")
    assert len(created) == 1
    assert output.stat().st_ino == created[0].st_ino
    assert stat.S_IMODE(output.stat().st_mode) == 0o700
    assert list(output.iterdir()) == []
    assert not list(tmp_path.glob(".osteosarc-*"))


def test_concurrent_valid_destination_is_reused_without_replacement(tmp_path, tiny_manifest, monkeypatch):
    path, manifest = tiny_manifest
    output = tmp_path / "export"
    original_publish = downloader._publish_no_replace
    created = []

    def publish_after_other_export(staged, destination):
        shutil.copytree(staged, destination)
        destination.chmod(0o700)
        created.append(destination.stat())
        original_publish(staged, destination)

    monkeypatch.setattr(downloader, "_publish_no_replace", publish_after_other_export)
    assert downloader.download_test_data(output, manifest_path=path, cache_root=tmp_path / "cache") == output
    assert downloader.verify_dataset(output, manifest) == output
    assert output.stat().st_ino == created[0].st_ino
    assert stat.S_IMODE(output.stat().st_mode) == 0o700
    assert not list(tmp_path.glob(".osteosarc-*"))


@pytest.mark.parametrize("missing_symbol", [False, True])
def test_unavailable_exclusive_rename_fails_without_publishing(tmp_path, tiny_manifest, monkeypatch, missing_symbol):
    path, _ = tiny_manifest
    output = tmp_path / "export"

    def unsupported(*args):
        downloader.ctypes.set_errno(errno.ENOTSUP)
        return -1

    libc = SimpleNamespace() if missing_symbol else SimpleNamespace(renameat2=unsupported, renamex_np=unsupported)
    monkeypatch.setattr(downloader.ctypes, "CDLL", lambda *args, **kwargs: libc)
    with pytest.raises(OSError) as error:
        downloader.download_test_data(output, manifest_path=path, cache_root=tmp_path / "cache")
    assert error.value.errno == errno.ENOTSUP
    assert not output.exists()
    assert not list(tmp_path.glob(".osteosarc-*"))


def test_failed_last_download_does_not_publish_partial_subset(tmp_path, tiny_manifest):
    path, manifest = tiny_manifest
    (tmp_path / "upstream/second.bam").write_bytes(b"wrong bytes")
    output, root = tmp_path / "export", tmp_path / "cache"
    with pytest.raises(FileValidationError):
        downloader.download_test_data(output, manifest_path=path, cache_root=root)
    assert not output.exists()
    assert (root / "objects/sha256" / downloader.cache_filename(manifest["assets"][0])).exists()
    assert not (root / "objects/sha256" / downloader.cache_filename(manifest["assets"][1])).exists()


def test_modified_or_foreign_output_is_never_replaced(tmp_path, tiny_manifest, monkeypatch):
    path, manifest = tiny_manifest
    output = downloader.download_test_data(tmp_path / "export", manifest_path=path, cache_root=tmp_path / "cache")
    (output / "first.bam").write_bytes(b"user change")
    monkeypatch.setattr(Cache, "fetch", forbid_fetch)
    with pytest.raises(FileValidationError):
        downloader.download_test_data(output, manifest_path=path)
    assert (output / "first.bam").read_bytes() == b"user change"
    (output / "notes.txt").write_text("user notes")
    with pytest.raises(ValueError):
        downloader.verify_dataset(output, manifest)
    assert (output / "notes.txt").read_text() == "user notes"


@pytest.mark.parametrize("filename", ["../escape.bam", "/escape.bam", "a/b.bam", "a\\b.bam", "manifest.json", "C:escape", ".."])
def test_unsafe_asset_paths_rejected_before_io(tmp_path, tiny_manifest, filename):
    path, manifest = tiny_manifest
    manifest["assets"][0]["filename"] = filename
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        downloader.download_test_data(tmp_path / "export", manifest_path=path, cache_root=tmp_path / "cache")
    assert not (tmp_path / "cache").exists()


def test_offline_and_verify_cli_do_not_fetch(tmp_path, tiny_manifest, monkeypatch):
    path, _ = tiny_manifest
    output = downloader.download_test_data(tmp_path / "export", manifest_path=path, cache_root=tmp_path / "cache")
    monkeypatch.setattr(Cache, "fetch", forbid_fetch)
    downloader.main(["--manifest", str(path), "--output", str(output), "--verify-only"])
    with pytest.raises(SystemExit) as error:
        downloader.main(["--manifest", str(path), "--offline", "--cache-root", str(tmp_path / "missing")])
    assert error.value.code == 1


def test_bundled_subset_covers_all_original_vaccine_loci_offline(monkeypatch):
    monkeypatch.setattr(Cache, "fetch", forbid_fetch)
    downloader.verify_dataset(DATA)
    identities = set()
    for case in MANIFEST["cases"]:
        allele = case["variant"].get("original_identity", case["variant"])
        identities.add(tuple(allele[k] for k in ("chrom", "pos", "ref", "alt")))
    assert len(MANIFEST["cases"]) == 49
    assert len(identities) == MANIFEST["original_vaccine_loci"] == 44
    assert len({i for i in identities if i[0] == "chr14" and i[1] in (101980529, 102030200)}) == 2


@pytest.mark.parametrize("case", MANIFEST["cases"], ids=lambda c: c["case_id"])
def test_original_bams_and_indexes_are_usable(case):
    with sid_reads("osteosarc/shared-v1/" + case["bam"]).open() as bam:
        assert bam.check_index()
        reads = list(bam)
        assert len(reads) == case["selected_record_count"]
        allele = case["variant"]
        contigs = [c for c in bam.references if c.removeprefix("chr").replace("MT", "M") == allele["chrom"].removeprefix("chr").replace("MT", "M")]
        assert len(contigs) == 1
        # Indexed native-coordinate retrieval must find original evidence, not
        # merely an apparently well-formed empty BAM or an index for another file.
        assert list(bam.fetch(contigs[0], allele["pos"] - 1, allele["pos"] + len(allele["ref"])))


def test_isovar_consumes_original_ntf3_compound_rna_without_network(tmp_path):
    from varcode import Variant
    from .osteosarc_fixture_support import indexed_genome

    # Existing pinned reference fixture includes NTF3; no genome download/cache
    # from the user's machine participates in this integration test.
    reference = sid_test_data() / "osteosarc/selection_validation/isovar/reference"
    genome = indexed_genome(reference, reference_name="GRCh38-shared-subset-ntf3-test",
        annotation_name="shared-subset-ensembl", annotation_version=87,
        cache_directory=tmp_path / "reference")
    case, = [c for c in MANIFEST["cases"] if c["variant"]["gene"] == "NTF3"]
    allele = case["variant"]
    variant = Variant(allele["chrom"].removeprefix("chr"), allele["pos"], allele["ref"], allele["alt"], ensembl=genome)
    with sid_reads("osteosarc/shared-v1/" + case["bam"]).open() as bam:
        evidence = ReadCollector().read_evidence_for_variant(variant, bam)
    assert evidence.alt_reads
    # At this locus the adjacent reference G is T in the RNA: focal A>G alone
    # cannot describe the observed compound AG>GT sequence.
    assert any(read.allele == "G" and read.suffix.startswith("T") for read in evidence.alt_reads)


def test_export_failure_preserves_destination_absence(tmp_path, tiny_manifest, monkeypatch):
    path, _ = tiny_manifest
    def fail_copy(*args, **kwargs):
        raise OSError("simulated disk full")
    monkeypatch.setattr(shutil, "copyfile", fail_copy)
    with pytest.raises(OSError):
        downloader.download_test_data(tmp_path / "export", manifest_path=path, cache_root=tmp_path / "cache")
    assert not (tmp_path / "export").exists()
    assert not list(tmp_path.glob(".osteosarc-*"))


def test_default_export_uses_only_the_sid_test_data(tmp_path, monkeypatch):
    import socket
    monkeypatch.setattr(downloader, "Cache", forbid_fetch)
    monkeypatch.setattr(socket.socket, "connect", forbid_fetch)
    missing_cache = tmp_path / "absent-cache"
    output = downloader.download_test_data(tmp_path / "bundle", cache_root=missing_cache)
    assert downloader.verify_dataset(output) == output
    assert not missing_cache.exists()


def test_offline_default_export_needs_the_cached_shared_reads(tmp_path, monkeypatch):
    from osteosarc import OfflineError
    monkeypatch.setenv("OSTEOSARC_CACHE", str(tmp_path / "empty-cache"))
    with pytest.raises(OfflineError):
        downloader.download_test_data(tmp_path / "export", offline=True)
    with pytest.raises(SystemExit) as error:
        downloader.main(["--output", str(tmp_path / "export"), "--verify-only"])
    assert error.value.code == 1
    assert not (tmp_path / "export").exists()


def test_inspect_dataset_reports_status_without_raising(tmp_path, tiny_manifest):
    """A caller that wants to know *what* is wrong gets a per-file status.

    verify_dataset answers "is it usable"; inspect_dataset answers "what is
    wrong with it", which is what the CLI prints and what a caller needs to
    decide between refetching and giving up.
    """
    path, manifest = tiny_manifest
    output = downloader.download_test_data(
        tmp_path / "export", manifest_path=path, cache_root=tmp_path / "cache")

    healthy = downloader.inspect_dataset(output, manifest)
    assert healthy["status"] == "available"
    assert healthy["verified"] is True
    assert {name: entry["status"] for name, entry in healthy["files"].items()} == {
        "first.bam": "available", "second.bam": "available"}
    assert all(entry["verified"] for entry in healthy["files"].values())

    (output / "first.bam").write_bytes(b"tampered")
    corrupt = downloader.inspect_dataset(output, manifest)
    assert corrupt["status"] == "corrupt"
    assert corrupt["verified"] is False
    assert corrupt["files"]["first.bam"] == dict(status="corrupt", verified=False)
    # The intact sibling is still reported as fine, rather than the whole
    # dataset collapsing to one exception.
    assert corrupt["files"]["second.bam"] == dict(status="available", verified=True)
    assert "first.bam" in corrupt["detail"]

    (output / "stray.txt").write_text("x")
    extra = downloader.inspect_dataset(output, manifest)
    assert extra["status"] == "unexpected_contents"
    assert "stray.txt" in extra["detail"]

    assert downloader.inspect_dataset(tmp_path / "absent", manifest)["status"] == "not_a_directory"


def test_inspect_dataset_rejects_symlinked_assets(tmp_path, tiny_manifest):
    """datacache follows symlinks, so the refusal has to stay ours.

    inspect_files would call a symlink pointing at correct bytes "available";
    an exported dataset must hold real files.
    """
    path, manifest = tiny_manifest
    output = downloader.download_test_data(
        tmp_path / "export", manifest_path=path, cache_root=tmp_path / "cache")
    real = (output / "first.bam").read_bytes()
    elsewhere = tmp_path / "elsewhere.bam"
    elsewhere.write_bytes(real)
    (output / "first.bam").unlink()
    (output / "first.bam").symlink_to(elsewhere)

    report = downloader.inspect_dataset(output, manifest)
    assert report["status"] == "symlink"
    assert "first.bam" in report["detail"]
    with pytest.raises(ValueError, match="must not be symlinks"):
        downloader.verify_dataset(output, manifest)


def test_cli_prints_one_json_object_for_success_and_failure(tmp_path, tiny_manifest, capsys):
    """The CLI's output is a contract, not debug prose.

    Nothing pinned it before, so the JSON could change shape unnoticed while
    exit codes stayed the same.
    """
    path, manifest = tiny_manifest
    output = tmp_path / "export"

    downloader.main(["--manifest", str(path), "--output", str(output),
                     "--cache-root", str(tmp_path / "cache")])
    exported = json.loads(capsys.readouterr().out)
    assert exported["action"] == "export"
    assert exported["status"] == "available"
    assert exported["dataset"] == manifest["dataset"]
    assert exported["data_version"] == manifest["data_version"]
    assert exported["assets"] == len(manifest["assets"])
    assert exported["output"] == str(output)
    assert exported["error"] is None
    assert {name: entry["status"] for name, entry in exported["files"].items()} == {
        "first.bam": "available", "second.bam": "available"}

    downloader.main(["--manifest", str(path), "--output", str(output), "--verify-only"])
    verified = json.loads(capsys.readouterr().out)
    assert verified["action"] == "verify"
    assert verified["status"] == "available"
    assert verified["error"] is None

    (output / "second.bam").write_bytes(b"tampered")
    with pytest.raises(SystemExit) as failure:
        downloader.main(["--manifest", str(path), "--output", str(output), "--verify-only"])
    assert failure.value.code == 1
    reported = json.loads(capsys.readouterr().out)
    assert reported["status"] == "corrupt"
    assert reported["files"]["second.bam"]["status"] == "corrupt"
    assert "second.bam" in reported["error"]


def test_cli_reports_cached_paths_without_an_output(tmp_path, tiny_manifest, capsys):
    """Without --output the result is cache locations, and says so in `status`.

    The old CLI inferred this by type-sniffing its own return value.
    """
    path, _ = tiny_manifest
    downloader.main(["--manifest", str(path), "--cache-root", str(tmp_path / "cache")])
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "cached"
    assert report["output"] is None
    assert sorted(report["cached_paths"]) == ["first.bam", "second.bam"]
    assert report["files"] == {}


def test_osteosarc_cache_environment_is_honoured(tmp_path, tiny_manifest, monkeypatch):
    """One invocation must not read bundled reads and custom assets from two roots.

    osteosarc resolves OSTEOSARC_CACHE first, so the downloader does too;
    previously it saw only OPENVAX_DATA_CACHE.
    """
    path, _ = tiny_manifest
    preferred = tmp_path / "osteosarc-root"
    fallback = tmp_path / "openvax-root"
    monkeypatch.setenv(downloader.OSTEOSARC_ENVIRONMENT, str(preferred))
    monkeypatch.setenv(downloader.CACHE_ENVIRONMENT, str(fallback))
    assert downloader.cache_root_for() == preferred

    downloader.download_test_data(manifest_path=path)
    assert (preferred / "objects" / "sha256").is_dir()
    assert not fallback.exists()

    # An explicit root still wins over both variables.
    explicit = tmp_path / "explicit"
    assert downloader.cache_root_for(explicit) == explicit
    monkeypatch.delenv(downloader.OSTEOSARC_ENVIRONMENT)
    assert downloader.cache_root_for() == fallback


def test_test_data_cli_is_installed_as_a_console_script():
    """`vaxrank-test-data` is the documented entry point, so keep it declared."""
    setup = (ROOT / "setup.py").read_text()
    assert "vaxrank-test-data = vaxrank.download_test_data:main" in setup


@pytest.mark.parametrize("envkey", ["OSTEOSARC_CACHE", "OPENVAX_DATA_CACHE"])
@pytest.mark.parametrize("form", ["absolute", "padded", "tilde", "blank"])
def test_cache_environment_normalization_matches_osteosarc(
        tmp_path, monkeypatch, platform_cache_home, envkey, form):
    from uuid import uuid4
    from osteosarc import Cache as OsteosarcCache

    for name in ("OSTEOSARC_CACHE", "OPENVAX_DATA_CACHE"):
        monkeypatch.delenv(name, raising=False)
    expected = tmp_path / "cache"
    if form == "absolute":
        value = str(expected)
    elif form == "padded":
        value = " \t" + str(expected) + " \n"
    elif form == "tilde":
        basename = "__vaxrank_cache_test_" + uuid4().hex
        value, expected = "~/" + basename, Path.home() / basename
    else:
        value, expected = " \t\n", platform_cache_home
    monkeypatch.setenv(envkey, value)
    assert downloader.cache_root_for() == OsteosarcCache().root == expected
    assert not expected.exists()
    assert not platform_cache_home.exists()


@pytest.mark.parametrize("blank", ["", " \t\n"])
def test_blank_preferred_cache_environment_uses_fallback(tmp_path, monkeypatch, blank):
    from osteosarc import Cache as OsteosarcCache

    expected = tmp_path / "fallback"
    monkeypatch.setenv("OSTEOSARC_CACHE", blank)
    monkeypatch.setenv("OPENVAX_DATA_CACHE", " " + str(expected) + " ")
    assert downloader.cache_root_for() == OsteosarcCache().root == expected
    assert not expected.exists()


def test_explicit_cache_path_is_preserved_despite_normalized_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("OSTEOSARC_CACHE", " ~/environment ")
    explicit = tmp_path / " literal path with spaces "
    assert downloader.cache_root_for(explicit) == explicit
    assert not explicit.exists()


def test_both_blank_cache_environment_values_use_platform_root(monkeypatch, platform_cache_home):
    from osteosarc import Cache as OsteosarcCache

    expected = platform_cache_home
    monkeypatch.setenv("OSTEOSARC_CACHE", " \t")
    monkeypatch.setenv("OPENVAX_DATA_CACHE", " \n")
    assert downloader.cache_root_for() == OsteosarcCache().root == expected
    assert not expected.exists()
