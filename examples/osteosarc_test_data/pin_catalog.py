#!/usr/bin/env python3
"""Pin the recipe's selected native assets/variants from an osteosarc snapshot.

Only metadata and BAM indexes are downloaded. Read selection remains the
separately reviewed selection.json.gz allowlist; this never resamples it.
"""

import argparse
from pathlib import Path

from osteosarc import Cache
from osteosarc.cohort_bundle import pin_catalog as pin_catalog

from vaxrank.sid_test_data import recipe_directory


RECIPE = recipe_directory()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, default=RECIPE)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    pin_catalog(args.recipe, Cache(args.cache))
