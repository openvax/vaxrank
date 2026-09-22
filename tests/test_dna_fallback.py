# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the DNA-only fallback (--allow-dna-only-fallback)."""

import pytest
import varcode

from vaxrank.mutant_protein_fragment import (
    MutantProteinFragment,
    find_mutation_region,
)


@pytest.fixture
def genome(human_genome_grch37):
    """Alias of the session-scoped GRCh37 handle for this module's tests."""
    return human_genome_grch37


# ---------------------------------------------------------------------------
# find_mutation_region unit tests
# ---------------------------------------------------------------------------

def test_find_mutation_region_substitution():
    ref = "ABCDEFGH"
    mut = "ABCXEFGH"
    start, end = find_mutation_region(ref, mut)
    assert start == 3
    assert end == 4


def test_find_mutation_region_insertion():
    ref = "ABCDEFGH"
    mut = "ABCXXDEFGH"
    start, end = find_mutation_region(ref, mut)
    assert start == 3
    assert end == 5  # two inserted AAs


def test_find_mutation_region_deletion():
    ref = "ABCDEFGH"
    mut = "ABCFGH"
    start, end = find_mutation_region(ref, mut)
    assert start == 3
    assert end == 3  # zero-length in mutant


def test_find_mutation_region_frameshift():
    ref = "ABCDEFGH"
    mut = "ABCXYZ"
    start, end = find_mutation_region(ref, mut)
    assert start == 3
    assert end == 6  # everything from pos 3 onward differs


def test_find_mutation_region_identical():
    ref = "ABCDEFGH"
    start, end = find_mutation_region(ref, ref)
    # start == len(ref), end == len(ref) — no mutation found
    assert start == len(ref)


def test_find_mutation_region_extension():
    ref = "ABCDEFGH"
    mut = "ABCDEFGHIJK"
    start, end = find_mutation_region(ref, mut)
    assert start == 8
    assert end == 11


# ---------------------------------------------------------------------------
# from_variant_dna tests
# ---------------------------------------------------------------------------

def test_from_variant_dna_substitution(genome):
    """BRAF V600E: a clean substitution should produce a valid fragment."""
    v = varcode.Variant(contig='7', start=140453136, ref='A', alt='T', ensembl=genome)
    frag = MutantProteinFragment.from_variant_dna(v, protein_sequence_length=35)
    assert frag is not None
    assert len(frag.amino_acids) == 35
    # The mutation region should be within the fragment
    assert frag.mutant_amino_acid_start_offset < frag.mutant_amino_acid_end_offset
    assert frag.mutant_amino_acid_end_offset <= len(frag.amino_acids)
    # DNA-only: zero reads
    assert frag.n_alt_reads == 0
    assert frag.n_alt_reads_supporting_protein_sequence == 0
    assert frag.gene_name == "BRAF"
    assert len(frag.supporting_reference_transcripts) == 1


def test_from_variant_dna_frameshift(genome):
    """TP53 frameshift should produce a valid fragment with extended mutation region."""
    v = varcode.Variant(contig='17', start=7577121, ref='G', alt='', ensembl=genome)
    frag = MutantProteinFragment.from_variant_dna(v, protein_sequence_length=35)
    assert frag is not None
    assert len(frag.amino_acids) <= 35
    assert frag.n_alt_reads == 0
    assert frag.gene_name == "TP53"


def test_from_variant_dna_noncoding_returns_none(genome):
    """An intronic variant should return None (no protein-modifying effect)."""
    # Use a known intronic variant
    v = varcode.Variant(contig='7', start=140500000, ref='A', alt='T', ensembl=genome)
    effects = v.effects()
    # Make sure it's actually non-coding
    coding = [e for e in effects if hasattr(e, 'transcript') and e.modifies_protein_sequence]
    if not coding:
        frag = MutantProteinFragment.from_variant_dna(v, protein_sequence_length=35)
        assert frag is None
    else:
        pytest.skip("Variant is unexpectedly coding")


def test_from_variant_dna_short_protein(genome):
    """When the protein is shorter than protein_sequence_length, the window
    should be clamped to the actual protein length."""
    v = varcode.Variant(contig='7', start=140453136, ref='A', alt='T', ensembl=genome)
    # Request a very long fragment
    frag = MutantProteinFragment.from_variant_dna(v, protein_sequence_length=5000)
    assert frag is not None
    # Fragment should be the full protein (766 AAs for BRAF)
    assert len(frag.amino_acids) < 5000
    assert len(frag.amino_acids) > 0


# ---------------------------------------------------------------------------
# DNA fallback integration with core_logic
# ---------------------------------------------------------------------------

def test_dna_fallback_disabled_returns_empty(genome):
    """With allow_dna_only_fallback=False, a variant without RNA should
    produce no vaccine peptides."""
    from types import SimpleNamespace
    from vaxrank.core_logic import vaccine_peptides_for_variant

    v = varcode.Variant(contig='7', start=140453136, ref='A', alt='T', ensembl=genome)
    isovar_result = SimpleNamespace(
        variant=v,
        passes_all_filters=False,
    )
    result = vaccine_peptides_for_variant(
        isovar_result=isovar_result,
        mhc_predictor=None,  # won't be reached
        allow_dna_only_fallback=False,
    )
    assert result == []


def test_dna_fallback_enabled_attempts_construction(genome):
    """With allow_dna_only_fallback=True, a variant without RNA should
    attempt DNA-based protein fragment construction."""
    from types import SimpleNamespace
    from unittest.mock import patch, MagicMock
    from vaxrank.core_logic import vaccine_peptides_for_variant

    v = varcode.Variant(contig='7', start=140453136, ref='A', alt='T', ensembl=genome)
    isovar_result = SimpleNamespace(
        variant=v,
        passes_all_filters=False,
    )

    # Mock predict_epitopes to avoid needing a real MHC predictor
    with patch('vaxrank.core_logic.predict_epitopes') as mock_predict:
        mock_predict.return_value = {}
        vaccine_peptides_for_variant(
            isovar_result=isovar_result,
            mhc_predictor=MagicMock(),
            allow_dna_only_fallback=True,
        )
    # With no epitope predictions, result is empty — but predict_epitopes
    # was called, proving the fallback path was taken
    mock_predict.assert_called_once()
    call_kwargs = mock_predict.call_args[1]
    frag = call_kwargs['protein_fragment']
    assert frag is not None
    assert frag.n_alt_reads == 0
    assert frag.gene_name == "BRAF"


@pytest.fixture
def fusion_candidates(genome, monkeypatch):
    """Real Varcode objects with synthetic proteins, not patient observations."""
    from varcode import MutantTranscript, StructuralVariant
    from varcode.effects import GeneFusion
    from varcode.effect_candidates import EffectCandidate

    transcript = genome.transcript_by_id("ENST00000003084")  # CFTR
    partner = genome.transcript_by_id("ENST00000269305")  # TP53
    variant = StructuralVariant(
        "7", 117171168, "BND", mate_contig="17", mate_start=7577121,
        genome=genome)
    reference = transcript.protein_sequence
    changed = reference[:40] + "WQWQWQWQWQ" + reference[50:]

    def build(primary_protein, alternative_protein=changed, alternative_partner=partner):
        def outcome(protein, partner):
            model = None if protein is None else MutantTranscript(
                reference_transcript=transcript, mutant_protein_sequence=protein,
                annotator_name="synthetic-selection-fixture")
            return GeneFusion(variant, transcript, partner, mutant_transcript=model)

        primary = outcome(primary_protein, partner)
        alternative = outcome(alternative_protein, alternative_partner)
        primary._extra_candidates = (EffectCandidate(
            alternative, source="selection-fixture",
            evidence={"junction": "synthetic", "rna_support": None}),)
        monkeypatch.setattr(StructuralVariant, "effects", lambda self: [primary])
        return variant, primary, alternative

    return build, reference, changed, transcript, partner


@pytest.mark.parametrize("primary_kind", ["unchanged", "unresolved", "empty", "truncated", "ambiguous"])
def test_dna_fusion_selects_usable_alternative(fusion_candidates, primary_kind):
    from vaxrank.varcode_effects import select_varcode_effect_outcome

    build, reference, changed, transcript, partner = fusion_candidates
    primary_protein = {"unchanged": reference, "unresolved": None,
                       "empty": "", "truncated": reference[:40],
                       "ambiguous": reference[:40] + "X" + reference[41:]}[primary_kind]
    variant, primary, alternative = build(primary_protein)
    # The generic priority selector intentionally retains Varcode's semantics.
    assert select_varcode_effect_outcome(primary) is primary
    fragment = MutantProteinFragment.from_variant_dna(variant, 35)
    assert fragment is not None
    assert fragment.amino_acids == changed[28:63]
    assert fragment.mutant_amino_acid_start_offset == 12
    assert fragment.mutant_amino_acid_end_offset == 22
    assert fragment.predicted_effect() is alternative
    assert fragment.predicted_effect("multi_outcome") is primary
    assert fragment.predicted_effect("most_likely") is primary
    assert fragment.global_start_pos() == 28
    assert fragment.supporting_reference_transcripts == [transcript, partner]
    selection = fragment.dna_effect_selection
    assert selection["transcript_id"] == transcript.id
    assert selection["partner_transcript_id"] == partner.id
    assert selection["five_prime_transcript_id"] == transcript.id
    assert selection["three_prime_transcript_id"] == partner.id
    assert selection["candidate_path"] == [dict(
        source="selection-fixture", evidence={"junction": "synthetic", "rna_support": None})]
    assert fragment.sequence_source == "varcode_translation"
    assert fragment.rna_evidence_method == fragment.rna_evidence_subject == ""
    assert fragment.n_alt_reads == fragment.n_alt_reads_supporting_protein_sequence == 0
    assert fragment.n_alt_fragments is None


@pytest.mark.parametrize("alternative_kind", ["unchanged", "unresolved", "empty", "truncated"])
def test_dna_fusion_without_usable_protein_returns_none(fusion_candidates, alternative_kind):
    build, reference, _changed, _transcript, _partner = fusion_candidates
    protein = {"unchanged": reference, "unresolved": None,
               "empty": "", "truncated": reference[:40]}[alternative_kind]
    variant, _, _ = build(None, protein)
    assert MutantProteinFragment.from_variant_dna(variant, 35) is None


def test_dna_fusion_selection_survives_json_and_subsequences(fusion_candidates):
    build, reference, _changed, _transcript, _partner = fusion_candidates
    variant, primary, alternative = build(reference)
    fragment = MutantProteinFragment.from_variant_dna(variant, 35)
    restored = MutantProteinFragment.from_json(fragment.to_json())
    assert restored.dna_effect_selection == fragment.dna_effect_selection
    assert restored.predicted_effect() is alternative
    assert restored.supporting_reference_transcripts == fragment.supporting_reference_transcripts
    for offset, part in restored.generate_subsequences(25):
        assert part.dna_effect_selection == restored.dna_effect_selection
        assert part.predicted_effect() is alternative
        assert part.global_start_pos() == restored.global_start_pos() + offset
    # A changed annotation must not attach a different effect to saved residues.
    primary._extra_candidates = ()
    assert restored.predicted_effect() is None


def test_dna_fusion_ties_use_partner_identity_not_producer_order(fusion_candidates, genome):
    build, _reference, changed, _transcript, partner = fusion_candidates
    other = genome.transcript_by_id("ENST00000357654")  # BRCA1
    variant, primary, alternative = build(changed, changed, other)
    assert partner.id < other.id
    assert MutantProteinFragment.from_variant_dna(variant, 35).predicted_effect() is primary
    from varcode.effect_candidates import EffectCandidate
    primary._primary_effects = (alternative,)
    primary._extra_candidates = (EffectCandidate(primary),)
    assert MutantProteinFragment.from_variant_dna(variant, 35).predicted_effect() is primary


def test_dna_fusion_prefers_longest_usable_protein(fusion_candidates):
    build, _reference, changed, _transcript, _partner = fusion_candidates
    variant, _primary, alternative = build(changed, changed + "WQWQ")
    assert MutantProteinFragment.from_variant_dna(variant, 35).predicted_effect() is alternative


def test_dna_unresolved_high_priority_does_not_hide_usable_lower_priority(fusion_candidates):
    from varcode.effects import Substitution
    from varcode.effect_candidates import EffectCandidate
    build, reference, _changed, transcript, _partner = fusion_candidates
    variant, primary, _alternative = build(None, None)
    substitution = Substitution(variant, transcript, 40, reference[40], "W")
    primary._extra_candidates = (EffectCandidate(substitution),)
    fragment = MutantProteinFragment.from_variant_dna(variant, 35)
    assert fragment is not None
    assert fragment.predicted_effect() is substitution


def test_partial_fusion_observation_is_not_compared_as_a_whole_protein(fusion_candidates):
    from dataclasses import replace
    build, _reference, changed, _transcript, _partner = fusion_candidates
    variant, _primary, alternative = build(None, changed[30:65])
    alternative.mutant_transcript = replace(
        alternative.mutant_transcript, evidence={"protein_completeness": "partial"})
    assert MutantProteinFragment.from_variant_dna(variant, 35) is None


def test_dna_wide_change_preserves_actual_window_coordinates(fusion_candidates):
    build, reference, _changed, _transcript, _partner = fusion_candidates
    changed = reference[:40] + "WQ" * 80 + reference[200:]
    variant, _, alternative = build(None, changed)
    fragment = MutantProteinFragment.from_variant_dna(variant, 35)
    start = fragment.global_start_pos()
    assert start == 103
    assert fragment.amino_acids == changed[start:start + 35]
    for offset, part in fragment.generate_subsequences(25):
        assert part.global_start_pos() == start + offset
        assert part.amino_acids == changed[start + offset:start + offset + 25]
        assert part.predicted_effect() is alternative
