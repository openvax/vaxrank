"""Pinned 1.8.1 original-read inputs for actual final-selection comparisons."""

from hashlib import sha256
import json
from pathlib import Path

from isovar import run_isovar
from isovar.cli import protein_sequence_creator_from_args, read_collector_from_args
from pyensembl import Genome
import pysam
from varcode import Variant

from vaxrank.cli import make_vaxrank_arg_parser
from vaxrank.cli.isovar_config_args import resolve_isovar_args
from vaxrank.vaccine_config import VaccineConfig


DATA = Path(__file__).parent / "data" / "osteosarc" / "selection_validation"


def load_selection_inputs(directory):
    root = DATA / "isovar"
    manifest = json.loads((root / "manifest.json").read_text())
    for filename, digest in manifest["files"].items():
        assert sha256((root / filename).read_bytes()).hexdigest() == digest, filename
    reference = root / "reference"
    metadata = json.loads((reference / "manifest.json").read_text())
    # Distinct from the older six-transcript and upstream stress fixture IDs:
    # pyensembl/varcode cache reference-derived information by identity.
    genome = Genome(
        reference_name="GRCh38-vaxrank-sid-selection-isovar181",
        annotation_name="sid-cohort-ensembl", annotation_version=87,
        gtf_path_or_url=str(reference / "reference.gtf.gz"),
        transcript_fasta_paths_or_urls=[str(reference / "reference.cdna.fa.gz")],
        protein_fasta_paths_or_urls=[str(reference / "reference.pep.fa.gz")],
        copy_local_files_to_cache=True, cache_directory_path=str(directory / "reference"))
    genome.index()
    return genome, {c["variant"]["variant_id"]: c for c in manifest["cases"]}, metadata


def reconstruct_selection(inputs, variant_id, length=25, flags=()):
    genome, cases, metadata = inputs
    case = cases[variant_id]
    record = case["variant"]
    variant = Variant(record["chrom"].removeprefix("chr"), record["pos"],
                      record["ref"], record["alt"], ensembl=genome)
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
