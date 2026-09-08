"""Original RNA -> Isovar policy -> Vaxrank fragments and targetable windows.

These are reconstruction contracts, not historical vaccine-ranking goldens.
"""

from unittest.mock import Mock

from isovar import run_isovar
from isovar.cli import protein_sequence_creator_from_args, read_collector_from_args
import pysam
import pytest

from vaxrank.cli import make_vaxrank_arg_parser
from vaxrank.cli.isovar_config_args import resolve_isovar_args
from vaxrank.core_logic import vaccine_peptides_for_variant
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_config import VaccineConfig

from .osteosarc_helpers import load_osteosarc


@pytest.fixture(scope="module")
def osteosarc(tmp_path_factory):
    return load_osteosarc(tmp_path_factory.mktemp("osteosarc_context"))


def reconstruct(osteosarc, sample, gene, flags=()):
    variants, bams, _ = osteosarc
    args = make_vaxrank_arg_parser().parse_args([
        "--vcf", "unused", "--bam", bams[sample], "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01", *flags])
    resolve_isovar_args(args, VaccineConfig(), {})
    creator = protein_sequence_creator_from_args(args)
    with pysam.AlignmentFile(bams[sample]) as bam:
        result, = run_isovar(
            variants=[variants[gene]], alignment_file=bam,
            read_collector=read_collector_from_args(args), protein_sequence_creator=creator)
    return result


@pytest.mark.parametrize("sample", ["bulk_star_t0", "ont_t1"])
@pytest.mark.parametrize("gene", ["DYNC1H1", "EXOC4", "H1-2", "GTF3C5"])
def test_original_rna_yields_expected_mutation_windows(osteosarc, sample, gene):
    result = reconstruct(osteosarc, sample, gene)
    protein = result.top_protein_sequence
    assert protein is not None
    assert protein.amino_acids in osteosarc[2][gene]
    fragment = MutantProteinFragment.from_isovar_result(result)
    assert fragment.amino_acids == protein.amino_acids
    assert fragment.n_alt_reads == result.num_alt_reads
    assert fragment.n_alt_fragments == result.num_alt_fragments
    assert fragment.n_alt_fragments_supporting_protein_sequence == protein.num_supporting_fragments
    assert fragment.n_alt_reads_supporting_protein_sequence == protein.num_supporting_reads
    assert fragment.supporting_reference_transcripts

    # Count full mutation-overlapping windows independently of Isovar's ranker.
    # Strict inequalities require both flanks of a zero-width deletion junction.
    start, end = protein.mutation_start_idx, protein.mutation_end_idx
    expected_offsets = {i for i in range(len(protein.amino_acids) - 25 + 1)
                        if i < end and i + 25 > start}
    actual_offsets = {i for i, candidate in fragment.generate_subsequences(25)
                      if len(candidate) == 25 and candidate.interval_overlaps_mutation(0, 25)}
    assert actual_offsets == expected_offsets
    assert actual_offsets  # a genuine full-length RNA-backed window
    for i in actual_offsets:
        assert fragment.amino_acids[i:i + 25] in osteosarc[2][gene]


@pytest.mark.parametrize("sample", ["bulk_star_t0", "ont_t1"])
def test_historical_context_cannot_claim_a_full_25mer(osteosarc, sample):
    result = reconstruct(osteosarc, sample, "DYNC1H1", [
        "--protein-sequence-length", "20", "--protein-sequence-preference", "support"])
    fragment = MutantProteinFragment.from_isovar_result(result)
    assert fragment.amino_acids in osteosarc[2]["DYNC1H1"]
    assert len(fragment) <= 20
    assert not [candidate for _, candidate in fragment.generate_subsequences(25)
                if len(candidate) == 25]


@pytest.mark.parametrize("sample,gene", [
    ("bulk_star_t0", "MAP2"), ("ont_t1", "MAP2"), ("bulk_star_t0", "PIP5K1A")])
def test_absent_alt_rna_does_not_become_a_selected_peptide(osteosarc, sample, gene):
    result = reconstruct(osteosarc, sample, gene)
    assert result.num_alt_reads == 0
    assert result.top_protein_sequence is None
    predictor = Mock(side_effect=AssertionError("No peptide should reach prediction"))
    assert vaccine_peptides_for_variant(result, predictor) == []
    assert not predictor.mock_calls
