"""Osteosarc-acquired historical cohorts for actual final-selection comparisons."""

import json

from isovar import run_isovar
from isovar.cli import protein_sequence_creator_from_args, read_collector_from_args
import pysam
from vaxrank.sid_test_data import sid_test_data, sid_variants

from vaxrank.cli import make_vaxrank_arg_parser
from vaxrank.cli.isovar_config_args import resolve_isovar_args
from vaxrank.vaccine_config import VaccineConfig

from .osteosarc_fixture_support import indexed_genome, verify_manifest_files


DATA = sid_test_data() / "osteosarc" / "selection_validation"


def load_selection_inputs(directory):
    root = DATA / "isovar"
    manifest = json.loads((root / "manifest.json").read_text())
    verify_manifest_files(root, manifest["files"], label="selection_validation/isovar")
    reference = root / "reference"
    metadata = json.loads((reference / "manifest.json").read_text())
    # Distinct from the older six-transcript and upstream stress fixture IDs:
    # pyensembl/varcode cache reference-derived information by identity.
    genome = indexed_genome(
        reference,
        reference_name="GRCh38-vaxrank-sid-selection-isovar181",
        annotation_name="sid-cohort-ensembl", annotation_version=87,
        cache_directory=directory / "reference")
    return genome, {c["variant"]["variant_id"]: c for c in manifest["cases"]}, metadata


def reconstruct_selection(inputs, variant_id, length=25, flags=(), variant=None):
    """Reconstruct a pinned case, optionally comparing another allele on its BAM."""
    genome, cases, metadata = inputs
    case = cases[variant_id]
    if variant is None:
        variant, = sid_variants([variant_id]).to_varcode(genome=genome, assembly="GRCh38")
    config = VaccineConfig(preferred_peptide_length=length,
                           min_peptide_length=length, max_peptide_length=length)
    bam_path = DATA / "isovar" / case["bam"]
    args = make_vaxrank_arg_parser().parse_args([
        "--vcf", "unused", "--bam", str(bam_path), "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-B*27:05", *flags])
    resolve_isovar_args(args, config, {})
    with pysam.AlignmentFile(str(bam_path)) as bam:
        result, = run_isovar(
            variants=[variant], alignment_file=bam,
            transcript_id_whitelist=metadata["variant_transcripts"][variant_id],
            read_collector=read_collector_from_args(args),
            protein_sequence_creator=protein_sequence_creator_from_args(args))
    return result, config
