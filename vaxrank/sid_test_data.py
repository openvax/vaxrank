"""The deliberately small, offline Sid regression bundle.

Read acquisition and source/variant identities come from osteosarc. This module
only materializes the selected package resources and verifies their integrity;
opening a bundle never downloads reads, references or dataset metadata.
"""

from functools import lru_cache
from importlib.resources import files
import json
from pathlib import Path, PurePosixPath
import tempfile
import zipfile

from osteosarc import ReadSubset, Variant, Variants, digest


BUNDLE_NAME = "sid-test-data.zip"


def extract_bundle(archive, destination):
    """Verify an explicit allowlist before materializing into an empty directory."""
    destination = Path(destination)
    if destination.is_symlink() or not destination.is_dir() or any(destination.iterdir()):
        raise ValueError("Bundle destination must be an empty, real directory")
    with zipfile.ZipFile(archive) as source:
        names = source.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate bundle member")
        manifest = json.loads(source.read("bundle.json"))
        if manifest.get("schema_version") != 1:
            raise ValueError("Unsupported Sid bundle schema")
        expected = manifest["files"]
        if set(names) != set(expected) | {"bundle.json"}:
            raise ValueError("Missing or unexpected bundle member")
        for name in names:
            path = PurePosixPath(name)
            if (path.is_absolute() or ".." in path.parts or not path.parts
                    or "\\" in name or ":" in name or str(path) != name):
                raise ValueError("Unsafe bundle path: " + name)
            info = source.getinfo(name)
            if info.is_dir() or (info.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("Bundle contains a non-file member")
            if name != "bundle.json" and info.file_size != expected[name]["size"]:
                raise ValueError("Bundle size mismatch: " + name)
        for name in names:
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(source.read(name))
            if name != "bundle.json" and digest(path) != expected[name]["sha256"]:
                raise ValueError("Bundle checksum mismatch: " + name)
    return destination


@lru_cache(maxsize=1)
def _materialized():
    directory = tempfile.TemporaryDirectory(prefix="vaxrank-sid-tests-")
    try:
        with files("vaxrank").joinpath("data", BUNDLE_NAME).open("rb") as archive:
            extract_bundle(archive, directory.name)
    except BaseException:
        directory.cleanup()
        raise
    # Holding the TemporaryDirectory keeps resources alive for this process.
    return directory


def sid_test_data():
    """Return the verified package subset, with no network or persistent cache."""
    return Path(_materialized().name)


def sid_variants(ids):
    """Native osteosarc entries with explicit published-allele snapshot identity."""
    provenance = json.loads((sid_test_data() / "provenance.json").read_text())
    variants = []
    for variant_id in ids:
        data = dict(provenance["variants"][variant_id])
        data["alleles"] = tuple(tuple(a) for a in data["alleles"])
        data["vaccines"] = tuple(data["vaccines"])
        data["pipelines"] = tuple(data["pipelines"])
        variants.append(Variant(**data))
    return Variants(variants, source={"snapshot_id": provenance["snapshot_id"], "corrections": False})


def sid_reads(relative_path):
    """Open a bundled BAM as an osteosarc ReadSubset with extraction lineage."""
    root = sid_test_data()
    provenance = json.loads((root / "provenance.json").read_text())
    cohort, = [c for c in provenance["cohorts"] if c["path"] == relative_path and c["format"] == "bam"]
    return ReadSubset(root / relative_path, root / (relative_path + ".bai"),
                      dict(cohort=cohort, records=len(cohort["records"]), scope="explicit_test_records",
                           source=provenance["sources"][cohort["source"]],
                           snapshot_id=provenance["snapshot_id"]))
