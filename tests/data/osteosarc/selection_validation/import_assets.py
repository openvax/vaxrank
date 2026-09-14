"""Explicit maintainer-only import; never invoked by tests or package import.

Usage: python tests/data/osteosarc/selection_validation/import_assets.py ISOVAR_REPO
Reads Git objects, not the possibly edited Isovar working tree.
"""

import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


COMMIT = "9297bb7bf19def6c33dd8295030de2b040d7dee3"
PREFIX = "tests/data/osteosarc/expansion/"
VARIANTS = {
    "DYNC1H1-chr14-101980529", "DYNC1H1-chr14-102030200",
    "EXOC4-chr7-133274996", "H1_2-chr6-26055824", "MAP2-chr2-209694768",
}


def main(repo):
    destination = Path(__file__).resolve().parent / "isovar"
    destination.mkdir(exist_ok=True)

    def read(path):
        return subprocess.check_output(["git", "-C", repo, "show", COMMIT + ":" + PREFIX + path])

    manifest_bytes = read("corpus/manifest.json")
    manifest = json.loads(manifest_bytes)
    cases = [c for c in manifest["cases"] if c["variant"]["variant_id"] in VARIANTS]
    # Identity, not cardinality: a duplicated case for one target plus a
    # missing case for another still totals len(VARIANTS), and would import
    # a fixture set that silently drops a variant and double-counts another.
    matched = [c["variant"]["variant_id"] for c in cases]
    if sorted(matched) != sorted(VARIANTS):
        raise ValueError(
            "Upstream cases do not match the target variants exactly; "
            "missing %s, duplicated %s" % (
                sorted(VARIANTS.difference(matched)) or "none",
                sorted({v for v in matched if matched.count(v) > 1}) or "none"))
    reference_path = "corpus/references/GRCh38/"
    reference_bytes = read(reference_path + "manifest.json")
    reference = json.loads(reference_bytes)
    files = {}

    def copy(path, relative, expected=None):
        data = read(path)
        digest = sha256(data).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError("Upstream checksum mismatch: " + path)
        output = destination / relative
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(data)
        files[relative] = digest

    for case in cases:
        # All original selected alignments, not the primary-only sensitivity run.
        for filename in (case["bam"], case["bam"] + ".bai"):
            copy("corpus/" + filename, filename, case["files"][filename])
    copy(reference_path + "manifest.json", "reference/manifest.json")
    for filename, digest in reference["files"].items():
        copy(reference_path + filename, "reference/" + filename, digest)
    index = json.loads(gzip.decompress(read("audit/counts-index.json.gz")))
    source_ids = {case["source_id"] for case in cases}
    sources = index["sources"]
    if isinstance(sources, dict):
        sources = list(sources.values())
    sources = [s for s in sources if s["source_id"] in source_ids]
    # Retain source metadata including header evidence verbatim. Do not infer
    # biological replication or cell identity from vendor/reprocessing labels.
    result = dict(
        upstream_repo="https://github.com/openvax/isovar", upstream_commit=COMMIT,
        upstream_release="1.8.1", upstream_manifest_sha256=sha256(manifest_bytes).hexdigest(),
        selection_baseline=manifest["selection_baseline"],
        scope="Selected-read fixtures are not unbiased VAF or cohort-prevalence estimates.",
        files=files, cases=cases, sources=sources)
    (destination / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("Imported", len(cases), "cases and", len(files), "checksummed files")


if __name__ == "__main__":
    main(sys.argv[1])
