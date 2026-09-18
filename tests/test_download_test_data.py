"""Offline downloader integrity, cache sharing and original-read contracts."""

from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil

from datacache import Cache, FileValidationError
from isovar import ReadCollector
import pysam
import pytest

from vaxrank import download_test_data as downloader


DATA = Path(__file__).parent / "data/osteosarc/shared-v1"
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


def test_checked_in_subset_covers_all_original_vaccine_loci_offline(monkeypatch):
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
    with pysam.AlignmentFile(DATA / case["bam"]) as bam:
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
    reference = Path(__file__).parent / "data/osteosarc/selection_validation/isovar/reference"
    genome = indexed_genome(reference, reference_name="GRCh38-shared-subset-ntf3-test",
        annotation_name="shared-subset-ensembl", annotation_version=87,
        cache_directory=tmp_path / "reference")
    case, = [c for c in MANIFEST["cases"] if c["variant"]["gene"] == "NTF3"]
    allele = case["variant"]
    variant = Variant(allele["chrom"].removeprefix("chr"), allele["pos"], allele["ref"], allele["alt"], ensembl=genome)
    with pysam.AlignmentFile(DATA / case["bam"]) as bam:
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
