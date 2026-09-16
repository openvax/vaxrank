#!/usr/bin/env python3
"""Rebuild the compact osteosarcoma complex-variant Vaxrank results."""

import argparse
from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from vaxrank.complex_variant_visualization import (  # noqa: E402
    generate_complex_variant_results,
)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root", default="examples/osteosarc_complex_results/runs")
    parser.add_argument("--timestamp")
    parser.add_argument("--combined-output")
    args = parser.parse_args(argv)
    source = Path(__file__).parent / "source" / "results.json"
    run = generate_complex_variant_results(
        source,
        args.output_root,
        timestamp=args.timestamp,
        combined_output=args.combined_output,
    )
    print(run)


if __name__ == "__main__":
    main()
