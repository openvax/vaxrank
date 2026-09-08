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
from vaxrank.core_logic import vaccine_peptides_for_variant, vaccine_peptides_from_epitopes
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_config import VaccineConfig

from .osteosarc_helpers import load_osteosarc
from .test_core_logic_config import _make_epitope


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


@pytest.mark.parametrize("sample,gene,window_count", [
    ("bulk_star_t0", "DYNC1H1", 25), ("ont_t1", "DYNC1H1", 0),
    ("bulk_star_t0", "EXOC4", 25), ("ont_t1", "EXOC4", 1),
    ("bulk_star_t0", "H1-2", 0), ("ont_t1", "H1-2", 6),
    ("bulk_star_t0", "GTF3C5", 24), ("ont_t1", "GTF3C5", 5),
])
def test_original_rna_yields_expected_mutation_windows(osteosarc, sample, gene, window_count):
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
    assert len(actual_offsets) == window_count
    for i in actual_offsets:
        assert fragment.amino_acids[i:i + 25] in osteosarc[2][gene]


def test_small_ont_fixture_cannot_silently_relax_relative_support_budget(osteosarc):
    result = reconstruct(osteosarc, "ont_t1", "DYNC1H1", [
        "--max-protein-sequences-per-variant", "0"])
    top = result.top_protein_sequence
    assert len(top.amino_acids) == 20
    assert top.num_supporting_fragments == 11
    assert result.num_alt_fragments == 16  # total alternate names != compatible names
    full = [p for p in result.sorted_protein_sequences if len(p.amino_acids) >= 25]
    assert full
    assert all(p.num_supporting_fragments / 11 < .85 for p in full)
    # An explicitly looser budget admits 37 aa / 9 names; context-first
    # admits 49 aa / 7 names. Neither is the default or a confidence claim.
    looser = reconstruct(osteosarc, "ont_t1", "DYNC1H1", [
        "--min-protein-sequence-support-fraction", "0.8"])
    context = reconstruct(osteosarc, "ont_t1", "DYNC1H1", [
        "--protein-sequence-preference", "context"])
    assert (len(looser.top_protein_sequence.amino_acids),
            looser.top_protein_sequence.num_supporting_fragments) == (37, 9)
    assert (len(context.top_protein_sequence.amino_acids),
            context.top_protein_sequence.num_supporting_fragments) == (49, 7)


def test_absolute_coverage_floor_is_independent_of_relative_budget(osteosarc):
    default = reconstruct(osteosarc, "bulk_star_t0", "H1-2")
    assert len(default.top_protein_sequence.amino_acids) == 24
    assert default.top_protein_sequence.num_supporting_fragments == 2
    # Three raw alternate read objects include paired mates; after merging
    # only two objects cover the mutant RNA. A floor of three rejects it,
    # even if we explicitly remove the relative budget with context-first.
    assert default.num_alt_reads == 3
    strict = reconstruct(osteosarc, "bulk_star_t0", "H1-2", [
        "--min-variant-sequence-coverage", "3", "--protein-sequence-preference", "context"])
    assert strict.num_alt_reads == 3
    assert strict.top_protein_sequence is None


@pytest.mark.parametrize("sample,gene,length", [
    ("ont_t1", "DYNC1H1", 20), ("bulk_star_t0", "H1-2", 24)])
def test_short_rna_context_is_not_selected_below_configured_minimum(osteosarc, sample, gene, length):
    result = reconstruct(osteosarc, sample, gene)
    fragment = MutantProteinFragment.from_isovar_result(result)
    assert len(fragment) == length
    # Deliberately synthetic scores isolate the length boundary; these are
    # not historical patient predictions or goldens for vaccine agreement.
    start = fragment.mutant_amino_acid_start_offset - 4
    epitope = _make_epitope(
        fragment.amino_acids[start:start + 9], ic50=100., wt_ic50=1000.,
        source_sequence=fragment.amino_acids, offset=start)
    assert vaccine_peptides_from_epitopes(
        result.variant, fragment, [epitope], vaccine_config=VaccineConfig()) == []
    permitted, = vaccine_peptides_from_epitopes(
        result.variant, fragment, [epitope],
        vaccine_config=VaccineConfig(min_peptide_length=length))
    assert permitted.mutant_protein_fragment.amino_acids == fragment.amino_acids
    assert len(permitted.mutant_protein_fragment) == length


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
