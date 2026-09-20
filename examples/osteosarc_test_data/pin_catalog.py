#!/usr/bin/env python3
"""Pin the recipe's selected native assets/variants from an osteosarc snapshot.

Only metadata and BAM indexes are downloaded. Read selection remains the
separately reviewed selection.json.gz allowlist; this never resamples it.
"""

import argparse
from dataclasses import asdict
import gzip
import json
from pathlib import Path

from osteosarc import Cache, Dataset, digest


RECIPE = Path(__file__).with_name("recipe")


def pin_catalog(recipe, cache):
    plan = json.loads(gzip.decompress((recipe / "selection.json.gz").read_bytes()))
    dataset = Dataset.open(plan["snapshot_name"], cache=cache, offline=cache.offline, corrections=False)
    if dataset.id != plan["snapshot_id"]:
        raise ValueError("Snapshot differs from the reviewed selection recipe")
    variants = dataset.variants("all")
    catalog = dict(snapshot=dataset.manifest, corrections=False, variants={}, assets={})
    for variant_id in sorted({v for c in plan["cohorts"] for v in c["variants"]}):
        catalog["variants"][variant_id] = {k: v for k, v in asdict(variants[variant_id]).items()
                                           if k != "annotations"}
    for url in sorted({c["source"] for c in plan["cohorts"]}):
        asset = dataset.asset(url)
        index = dataset.asset(asset.index_urls[0])
        path = dataset.download(index)
        receipt = cache.fetch(index.url, sha256=digest(path), size=path.stat().st_size)
        catalog["assets"][url] = dict(asset=asdict(asset), index=asdict(index), index_receipt=receipt.to_dict())
    (recipe / "catalog.json").write_text(json.dumps(catalog, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, default=RECIPE)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    pin_catalog(args.recipe, Cache(args.cache))
