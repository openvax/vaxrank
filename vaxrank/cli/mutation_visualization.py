"""Command-line entry point for mutation-context figures."""

import argparse

from ..mutation_visualization import generate_mutation_figures


def make_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Render reference, annotation-only, and RNA-assembled transcript and "
            "protein context from a Vaxrank --output-isovar-csv file."
        ))
    parser.add_argument("input_isovar_csv", help="CSV produced by --output-isovar-csv")
    parser.add_argument(
        "--output-root", default="mutation-figures",
        help="Parent directory for the timestamped run (default: %(default)s)")
    parser.add_argument(
        "--variant", action="append", default=[],
        help="Case-insensitive row substring to include; may be repeated")
    parser.add_argument(
        "--format", action="append", choices=("svg", "pdf", "png"), dest="formats",
        help="Output format; may be repeated (default: svg, pdf, and png)")
    parser.add_argument(
        "--png-scale", type=float, default=3,
        help="PNG scale relative to 1200x700 (default: 3, or 3600x2100)")
    parser.add_argument(
        "--timestamp",
        help="Reproducible UTC run name in YYYY-MM-DDTHHMMSSZ form")
    return parser


def main(argv=None):
    args = make_parser().parse_args(argv)
    run_directory = generate_mutation_figures(
        args.input_isovar_csv,
        args.output_root,
        variants=args.variant,
        formats=args.formats or ("svg", "pdf", "png"),
        timestamp=args.timestamp,
        png_scale=args.png_scale,
    )
    print(run_directory)


if __name__ == "__main__":
    main()
