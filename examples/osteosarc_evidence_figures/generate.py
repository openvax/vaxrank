#!/usr/bin/env python3
"""Rebuild the osteosarc.com cross-platform evidence figures."""

import argparse
from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from vaxrank.evidence_visualization import generate_evidence_figures  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root", default="examples/osteosarc_evidence_figures/runs")
    parser.add_argument("--timestamp")
    args = parser.parse_args(argv)
    source = Path(__file__).parent / "source" / "evidence.json"
    run = generate_evidence_figures(
        source, args.output_root, timestamp=args.timestamp)
    print(run)


if __name__ == "__main__":
    main()
