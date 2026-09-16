#!/usr/bin/env python3
"""Rebuild the osteosarc.com mutation-context figure examples."""

import argparse
import logging
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from isovar import isovar_results_to_dataframe, run_isovar  # noqa: E402
from isovar.cli import (  # noqa: E402
    protein_sequence_creator_from_args,
    read_collector_from_args,
)
import pysam  # noqa: E402

from tests.osteosarc_helpers import load_osteosarc  # noqa: E402
from vaxrank.cli import make_vaxrank_arg_parser  # noqa: E402
from vaxrank.cli.isovar_config_args import resolve_isovar_args  # noqa: E402
from vaxrank.mutation_visualization import generate_mutation_figures  # noqa: E402
from vaxrank.vaccine_config import VaccineConfig  # noqa: E402


CASES = [
    {
        "sample": "bulk_star_t0",
        "gene": "DYNC1H1",
        "label": "DYNC1H1",
        "note": "T0 short-read RNA: assembly confirms p.Val314Ile",
        "url": "https://osteosarc.com/variant/DYNC1H1-chr14-101980529/",
    },
    {
        "sample": "bulk_star_t0",
        "gene": "EXOC4",
        "label": "EXOC4",
        "note": "T0 short-read RNA: assembly confirms p.Ser34Ile",
        "url": "https://osteosarc.com/variant/EXOC4-chr7-133274996/",
    },
    {
        "sample": "ont_t1",
        "gene": "H1-2",
        "label": "H1-2",
        "note": "T1 long-read RNA: assembly refines a repetitive deletion junction",
        "url": "https://osteosarc.com/variant/H1_2-chr6-26055824/",
    },
    {
        "sample": "ont_t1",
        "gene": "GTF3C5",
        "label": "GTF3C5",
        "note": "T1 long-read RNA: assembly resolves a repetitive glutamate deletion context",
        "url": "https://osteosarc.com/variant/GTF3C5-chr9-133057893/",
    },
    {
        "sample": "bulk_star_t0",
        "gene": "MAP2",
        "label": "MAP2",
        "note": "T0 short-read RNA: annotation predicts a frameshift, but no alternate fragments assemble",
        "url": "https://osteosarc.com/variant/MAP2-chr2-209694768/",
    },
    {
        "sample": "bulk_star_t0",
        "gene": "PIP5K1A",
        "label": "PIP5K1A",
        "note": "T0 short-read RNA: annotation predicts a frameshift, but no alternate fragments assemble",
        "url": "https://osteosarc.com/variant/PIP5K1A-chr1-151242178/",
    },
]


def reconstruct(osteosarc, sample, gene):
    variants, bams, _ = osteosarc
    args = make_vaxrank_arg_parser().parse_args([
        "--vcf", "unused",
        "--bam", bams[sample],
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    resolve_isovar_args(args, VaccineConfig(), {})
    with pysam.AlignmentFile(bams[sample]) as alignment_file:
        result, = run_isovar(
            variants=[variants[gene]],
            alignment_file=alignment_file,
            read_collector=read_collector_from_args(args),
            protein_sequence_creator=protein_sequence_creator_from_args(args),
        )
    return result


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        default="examples/osteosarc_mutation_figures",
        help="Example directory (default: %(default)s)",
    )
    parser.add_argument(
        "--timestamp",
        help="UTC run name in YYYY-MM-DDTHHMMSSZ form; defaults to current time",
    )
    args = parser.parse_args(argv)
    logging.disable(logging.INFO)
    output_root = Path(args.output_root)
    source_directory = output_root / "source"
    source_directory.mkdir(parents=True, exist_ok=True)
    source_csv = source_directory / "osteosarc-isovar.csv"
    with TemporaryDirectory(prefix="vaxrank-osteosarc-figures-") as directory:
        osteosarc = load_osteosarc(Path(directory))
        dataframe = isovar_results_to_dataframe([
            reconstruct(osteosarc, case["sample"], case["gene"])
            for case in CASES
        ])
    dataframe["figure_label"] = [case["label"] for case in CASES]
    dataframe["figure_note"] = [case["note"] for case in CASES]
    dataframe["figure_source_url"] = [case["url"] for case in CASES]
    dataframe["figure_rna_sample"] = [case["sample"] for case in CASES]
    dataframe.to_csv(source_csv, index=False)
    run = generate_mutation_figures(
        source_csv,
        output_root / "runs",
        timestamp=args.timestamp,
    )
    print(run)


if __name__ == "__main__":
    main()
