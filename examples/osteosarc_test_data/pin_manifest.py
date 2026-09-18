"""Maintainer-only: pin the upstream Isovar subset from immutable Git objects.

Ordinary downloads need only the packaged manifest, not an Isovar checkout.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


UPSTREAM_COMMIT = "0cad5b275c852263a1c77722aa463fe5568b2c76"
UPSTREAM_DIRECTORY = "tests/data/osteosarc/expansion/corpus"


def pin_manifest(repository, output):
    def blob(path):
        return subprocess.check_output(["git", "-C", str(repository), "show", UPSTREAM_COMMIT + ":" + path])

    raw = blob(UPSTREAM_DIRECTORY + "/manifest.json")
    upstream = json.loads(raw)
    assets, cases = [], []
    for case in upstream["cases"]:
        for filename in (case["bam"], case["bam"] + ".bai"):
            data = blob(UPSTREAM_DIRECTORY + "/" + filename)
            sha256 = hashlib.sha256(data).hexdigest()
            if sha256 != case["files"][filename]:
                raise ValueError("Upstream manifest and Git asset disagree: " + filename)
            assets.append(dict(filename=filename, sha256=sha256, size_bytes=len(data),
                url="https://raw.githubusercontent.com/openvax/isovar/%s/%s/%s" % (UPSTREAM_COMMIT, UPSTREAM_DIRECTORY, filename)))
        record = {key: case[key] for key in ("case_id", "bam", "variant", "source_id", "source_url", "selection",
                  "source_region_sha256", "source_receipt_sha256", "selection_row_sha256")}
        record["selected_record_count"] = len(case["selected_records"])
        cases.append(record)
    manifest = dict(schema_version=1, dataset="osteosarc", data_version="vaccine-rna-v1",
        species="Homo sapiens", taxon_id=9606, source_data_license="CC0-1.0",
        license_source="https://registry.opendata.aws/sid-osteosarc/",
        scope=upstream["scope"], original_vaccine_loci=44,
        interpretation="Deliberately selected regression reads, not full-source coverage, independent platform samples, somatic truth or validated vaccine targets.",
        cache_namespace="openvax", cache_object_layout="objects/sha256/{sha256}{original_suffixes}",
        upstream=dict(repository="https://github.com/openvax/isovar", commit=UPSTREAM_COMMIT,
            manifest_url="https://raw.githubusercontent.com/openvax/isovar/%s/%s/manifest.json" % (UPSTREAM_COMMIT, UPSTREAM_DIRECTORY),
            manifest_sha256=hashlib.sha256(raw).hexdigest(), generator_path="tests/data/osteosarc/expansion/fixtures.py"),
        assets=assets, cases=cases)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print("Pinned %d cases / %d assets / %d bytes" % (len(cases), len(assets), sum(a["size_bytes"] for a in assets)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isovar-repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("vaxrank/data/osteosarc-test-data-v1.json"))
    args = parser.parse_args()
    pin_manifest(args.isovar_repo, args.output)
