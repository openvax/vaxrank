"""The deliberately small Sid regression data, built from shared test reads.

``vaxrank/data/sid-recipe`` lists every read record the Sid tests use, with the
fixture headers, fusion templates and non-read reference, expectation and
prediction files. The reads come from openvax-v1, the OpenVax libraries'
shared Sid test data published by osteosarc (iskandr/osteosarc#56). The first
build downloads and verifies openvax-v1 (28 MB) into the osteosarc cache
(``OSTEOSARC_CACHE``, else the shared OpenVax cache); later builds are offline.
"""

from functools import lru_cache
import gzip
from importlib.resources import files
import json
from pathlib import Path
import shutil
import tempfile

import osteosarc
from osteosarc import ReadSubset, Variant, Variants, digest
from osteosarc.cohort_bundle import select_records, update_manifests, write_cohort, write_json
from osteosarc.shared import published


READS = "openvax-v1"
MEMBER_PREFIX = "vaxrank/"


def recipe_directory():
    """The reviewed Sid recipe shipped with Vaxrank."""
    return Path(str(files("vaxrank").joinpath("data", "sid-recipe")))


def build_sid_test_data(destination):
    """Write the Sid test files into an empty directory and return it.

    Each cohort's records are selected by exact record digest from its
    openvax-v1 member and written with the recipe's reviewed header and order,
    so the files match the recipe's allowlist byte for byte.
    """
    recipe, root = recipe_directory(), Path(destination)
    if any(root.iterdir()):
        raise ValueError("Sid test data destination must be empty")
    plan = json.loads(gzip.decompress((recipe / "selection.json.gz").read_bytes()))
    bundle = Path(osteosarc.fetch_bundle(READS))
    members = json.loads((bundle / "manifest.json").read_text())["members"]
    sources = json.loads((bundle / "recipe.json").read_text())["sources"]
    shutil.copytree(recipe / "support", root, dirs_exist_ok=True)
    provenance_sources = {}
    for cohort in plan["cohorts"]:
        member = MEMBER_PREFIX + cohort["path"]
        source_id = members[member]["source"]
        identity = sources[source_id]["identity"]
        if identity["url"] != cohort["source"]:
            raise ValueError("%s comes from %s, not %s" % (member, identity["url"], cohort["source"]))
        provenance_sources[cohort["source"]] = dict(member_source=source_id, identity=identity)
        path = Path(osteosarc.bundle_file(bundle, member))
        records = select_records(ReadSubset(path, Path(str(path) + ".bai"), {}), cohort)
        write_cohort(root, cohort, records, recipe)
    update_manifests(root)
    catalog = json.loads((recipe / "catalog.json").read_text())
    write_json(root / "provenance.json", dict(
        schema_version=2, snapshot_id=plan["snapshot_id"], corrections=False,
        osteosarc_version=osteosarc.__version__,
        reads=dict(bundle=READS, manifest_sha256=published(READS)["manifest_sha256"]),
        recipe_files={p.relative_to(recipe).as_posix(): digest(p)
                      for p in sorted(recipe.rglob("*")) if p.is_file()},
        variants=catalog["variants"], sources=provenance_sources,
        cohorts=[{k: v for k, v in c.items() if k not in ("header", "record_metadata")}
                 for c in plan["cohorts"]]))
    return root


@lru_cache(maxsize=1)
def _materialized():
    directory = tempfile.TemporaryDirectory(prefix="vaxrank-sid-tests-")
    try:
        build_sid_test_data(directory.name)
    except BaseException:
        directory.cleanup()
        raise
    # Holding the TemporaryDirectory keeps the files alive for this process.
    return directory


def sid_test_data():
    """Return the Sid test files, built once per process."""
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
    """Open a Sid test BAM as an osteosarc ReadSubset with its lineage."""
    root = sid_test_data()
    provenance = json.loads((root / "provenance.json").read_text())
    cohort, = [c for c in provenance["cohorts"] if c["path"] == relative_path and c["format"] == "bam"]
    return ReadSubset(root / relative_path, root / (relative_path + ".bai"),
                      dict(cohort=cohort, records=len(cohort["records"]), scope="explicit_test_records",
                           source=provenance["sources"][cohort["source"]], reads=provenance["reads"],
                           snapshot_id=provenance["snapshot_id"]))
