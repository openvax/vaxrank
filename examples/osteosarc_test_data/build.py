#!/usr/bin/env python3
"""Regenerate the reviewed Sid cohorts through Osteosarc's shared implementation."""
import argparse
from pathlib import Path
from osteosarc import Cache
from osteosarc.cohort_bundle import (
    generate_cohort_bundle as build,
    record_digest as record_digest,
    retrieval_regions as retrieval_regions,
    select_records as select_records,
    write_cohort as write_cohort,
    update_manifests as update_manifests,
)

RECIPE = Path(__file__).with_name("recipe")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, default=RECIPE)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("vaxrank/data/sid-test-data.zip"))
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--panel-recipe", type=Path, help="Shared Osteosarc v1 panel recipe")
    parser.add_argument("--panel-source", action="append", default=[], metavar="ID=LOCAL_BAM")
    args = parser.parse_args()
    if args.panel_recipe:
        from osteosarc.bundles import generate_panel
        generate_panel(args.panel_recipe, args.output, sources=args.panel_source,
                       cache=Cache(args.cache, offline=args.offline))
        return

    build(args.recipe, Cache(args.cache, offline=args.offline), args.output)


if __name__ == "__main__":
    main()
