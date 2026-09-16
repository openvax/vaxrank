"""Command-line interface for cross-platform sequence-evidence figures."""

import argparse

from ..evidence_visualization import generate_evidence_figures


def make_parser():
    parser = argparse.ArgumentParser(
        description="Render publication-ready DNA/RNA sequence-evidence figures")
    parser.add_argument("input_json", help="Evidence record JSON")
    parser.add_argument("--output-root", default="evidence-figures")
    parser.add_argument(
        "--format", action="append", choices=("svg", "pdf", "png"),
        dest="formats", help="Repeat to select formats (default: all)")
    parser.add_argument("--timestamp", help="UTC run name: YYYY-MM-DDTHHMMSSZ")
    parser.add_argument("--png-scale", type=float, default=3)
    return parser


def main(argv=None):
    args = make_parser().parse_args(argv)
    run = generate_evidence_figures(
        args.input_json,
        args.output_root,
        formats=args.formats or ("svg", "pdf", "png"),
        timestamp=args.timestamp,
        png_scale=args.png_scale,
    )
    print(run)
