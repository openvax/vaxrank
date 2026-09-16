#!/usr/bin/env python3
"""Rebuild the osteosarc.com mutation-context figure examples."""

import argparse
import logging
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from isovar import isovar_results_to_dataframe, run_isovar  # noqa: E402
from isovar.cli import (  # noqa: E402
    protein_sequence_creator_from_args,
    read_collector_from_args,
)
import pysam  # noqa: E402
import pandas as pd  # noqa: E402
from varcode import Variant  # noqa: E402

from tests.osteosarc_helpers import load_osteosarc  # noqa: E402
from tests.osteosarc_selection_helpers import (  # noqa: E402
    load_selection_inputs,
    reconstruct_selection,
)
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

SELECTION_CASES = [
    {
        "variant_id": "DYNC1H1-chr14-102030200",
        "label": "DYNC1H1",
        "note": "T0 short-read RNA: assembly confirms the second DYNC1H1 locus, p.Gln3267His",
        "url": "https://osteosarc.com/variant/DYNC1H1-chr14-102030200/",
        "sample": "bulk_bostongene_t0",
    },
]

NAV2_CASE = {
    "variant_id": "NAV2-chr11-20080142",
    "transcript_id": "ENST00000396087",
    "label": "NAV2",
    "note": "DNA catalog variant: transcript-specific p.Ala1809Val; RNA not assessed in the pinned fixture",
    "url": "https://osteosarc.com/variants/",
    "sample": "not_assessed",
}


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


def _transcript_alt(effect):
    alt = effect.variant.alt
    if effect.transcript.strand == "-":
        from isovar.dna import reverse_complement_dna

        alt = reverse_complement_dna(alt)
    return alt


def cdna_context_columns(result):
    """Build transcript-oriented reference, annotation, and assembly tracks."""
    from isovar.variant_helpers import (
        interbase_range_affected_by_variant_on_transcript,
    )

    effect = result.predicted_effect
    variant_start, variant_end = interbase_range_affected_by_variant_on_transcript(
        effect.variant, effect.transcript)
    reference = effect.transcript.sequence
    alt = _transcript_alt(effect)
    annotation = reference[:variant_start] + alt + reference[variant_end:]
    columns = {
        "reference_cdna_sequence": reference,
        "reference_cdna_variant_start": variant_start,
        "reference_cdna_variant_end": variant_end,
        "annotation_cdna_sequence": annotation,
        "annotation_cdna_variant_start": variant_start,
        "annotation_cdna_variant_end": variant_start + len(alt),
        "rna_assembled_cdna_sequence": None,
        "rna_assembled_cdna_variant_start": None,
        "rna_assembled_cdna_variant_end": None,
    }
    if result.top_protein_sequence:
        translations = sorted(
            result.top_protein_sequence.translations,
            key=lambda translation: (
                effect.transcript_id not in {
                    transcript.id
                    for transcript in translation.reference_context.transcripts
                },
                translation.cdna_sequence,
            ),
        )
        translation = translations[0]
        columns.update({
            "rna_assembled_cdna_sequence": translation.cdna_sequence,
            "rna_assembled_cdna_variant_start": (
                translation.variant_cdna_interval_start),
            "rna_assembled_cdna_variant_end": translation.variant_cdna_interval_end,
        })
    return columns


def annotation_only_row(effect):
    """Represent a DNA catalog variant without implying RNA-negative evidence."""
    return {
        "variant": effect.variant.short_description,
        "predicted_effect": effect.short_description,
        "predicted_effect_gene_name": effect.gene_name,
        "predicted_effect_gene_id": effect.gene_id,
        "predicted_effect_transcript_id": effect.transcript_id,
        "predicted_effect_transcript_name": effect.transcript_name,
        "predicted_effect_original_protein_sequence": (
            effect.original_protein_sequence),
        "predicted_effect_mutant_protein_sequence": effect.mutant_protein_sequence,
        "predicted_effect_aa_mutation_start_offset": (
            effect.aa_mutation_start_offset),
        "protein_sequence": None,
        "protein_sequence_gene_names": None,
        "protein_sequence_gene_ids": None,
        "protein_sequence_transcript_names": None,
        "protein_sequence_transcript_ids": None,
        "num_alt_fragments": 0,
        "num_ref_fragments": 0,
        "num_other_fragments": 0,
        "num_fragments_supporting_top_protein_sequence": 0,
    }


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
        results = [
            reconstruct(osteosarc, case["sample"], case["gene"])
            for case in CASES
        ]
        selection_inputs = load_selection_inputs(Path(directory) / "selection")
        for case in SELECTION_CASES:
            result, _ = reconstruct_selection(
                selection_inputs, case["variant_id"], length=25)
            results.append(result)
        dataframe = isovar_results_to_dataframe(results)
        cdna_rows = [cdna_context_columns(result) for result in results]

        genome, _, _ = selection_inputs
        nav2 = Variant("11", 20080142, "C", "T", ensembl=genome)
        nav2_effect, = [
            effect for effect in nav2.effects()
            if effect.transcript_id == NAV2_CASE["transcript_id"]
            and effect.modifies_protein_sequence
        ]
        nav2_row = annotation_only_row(nav2_effect)
        dataframe = pd.concat(
            [dataframe, pd.DataFrame([nav2_row])], ignore_index=True)
        cdna_rows.append(cdna_context_columns(SimpleNamespace(
            predicted_effect=nav2_effect,
            top_protein_sequence=None,
        )))

    all_cases = [*CASES, *SELECTION_CASES, NAV2_CASE]
    for column in cdna_rows[0]:
        dataframe[column] = [row[column] for row in cdna_rows]
    dataframe["figure_label"] = [case["label"] for case in all_cases]
    dataframe["figure_note"] = [case["note"] for case in all_cases]
    dataframe["figure_source_url"] = [case["url"] for case in all_cases]
    dataframe["figure_rna_sample"] = [case["sample"] for case in all_cases]
    dataframe["figure_rna_status"] = [
        "not_assessed" if case["sample"] == "not_assessed" else "assessed"
        for case in all_cases
    ]
    dataframe.to_csv(source_csv, index=False)
    run = generate_mutation_figures(
        source_csv,
        output_root / "runs",
        timestamp=args.timestamp,
    )
    print(run)


if __name__ == "__main__":
    main()
