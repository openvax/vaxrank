"""The deliberately small, offline Sid regression bundle.

Read acquisition and source/variant identities come from osteosarc. This module
only materializes the selected package resources and verifies their integrity;
opening a bundle never downloads reads, references or dataset metadata.
"""

from functools import lru_cache
from importlib.resources import files
import json
from pathlib import Path
import tempfile

from osteosarc import ReadSubset, Variant, Variants
from osteosarc.cohort_bundle import extract_bundle as extract_bundle


BUNDLE_NAME = "sid-test-data.zip"


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
