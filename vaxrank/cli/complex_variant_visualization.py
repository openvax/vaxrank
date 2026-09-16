"""Command-line interface for assembly-dependent Vaxrank result figures."""

import argparse

from ..complex_variant_visualization import generate_complex_variant_results


def make_parser():
    parser = argparse.ArgumentParser(
        description="Render complex-variant Vaxrank decisions and rankings")
    parser.add_argument("input_json", help="Complex-variant result JSON")
    parser.add_argument("--output-root", default="complex-variant-results")
    parser.add_argument("--timestamp", help="UTC run name: YYYY-MM-DDTHHMMSSZ")
    parser.add_argument("--png-scale", type=float, default=3)
    parser.add_argument("--combined-output", help="Optional copy of the combined PDF")
    return parser


def main(argv=None):
    args = make_parser().parse_args(argv)
    run = generate_complex_variant_results(
        args.input_json,
        args.output_root,
        timestamp=args.timestamp,
        png_scale=args.png_scale,
        combined_output=args.combined_output,
    )
    print(run)
