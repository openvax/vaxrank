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
from pyensembl import EnsemblRelease
import varcode

from vaxrank.mutant_protein_fragment import MutantProteinFragment, _find_mutation_region


genome = EnsemblRelease(75)


# ---------------------------------------------------------------------------
# _find_mutation_region unit tests
# ---------------------------------------------------------------------------

def test_find_mutation_region_substitution():
    ref = "ABCDEFGH"
    mut = "ABCXEFGH"
    start, end = _find_mutation_region(ref, mut)
    assert start == 3
    assert end == 4


def test_find_mutation_region_insertion():
    ref = "ABCDEFGH"
    mut = "ABCXXDEFGH"
    start, end = _find_mutation_region(ref, mut)
    assert start == 3
    assert end == 5  # two inserted AAs


def test_find_mutation_region_deletion():
    ref = "ABCDEFGH"
    mut = "ABCFGH"
    start, end = _find_mutation_region(ref, mut)
    assert start == 3
    assert end == 3  # zero-length in mutant


def test_find_mutation_region_frameshift():
    ref = "ABCDEFGH"
    mut = "ABCXYZ"
    start, end = _find_mutation_region(ref, mut)
    assert start == 3
    assert end == 6  # everything from pos 3 onward differs


def test_find_mutation_region_identical():
    ref = "ABCDEFGH"
    start, end = _find_mutation_region(ref, ref)
    # start == len(ref), end == len(ref) — no mutation found
    assert start == len(ref)


def test_find_mutation_region_extension():
    ref = "ABCDEFGH"
    mut = "ABCDEFGHIJK"
    start, end = _find_mutation_region(ref, mut)
    assert start == 8
    assert end == 11


# ---------------------------------------------------------------------------
# from_variant_dna tests
# ---------------------------------------------------------------------------

def test_from_variant_dna_substitution():
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


def test_from_variant_dna_frameshift():
    """TP53 frameshift should produce a valid fragment with extended mutation region."""
    v = varcode.Variant(contig='17', start=7577121, ref='G', alt='', ensembl=genome)
    frag = MutantProteinFragment.from_variant_dna(v, protein_sequence_length=35)
    assert frag is not None
    assert len(frag.amino_acids) <= 35
    assert frag.n_alt_reads == 0
    assert frag.gene_name == "TP53"


def test_from_variant_dna_noncoding_returns_none():
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


def test_from_variant_dna_short_protein():
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

def test_dna_fallback_disabled_returns_empty():
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


def test_dna_fallback_enabled_attempts_construction():
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
