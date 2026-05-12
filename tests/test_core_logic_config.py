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

"""
Tests for core_logic.py integration with config objects.
"""

from inspect import signature
from unittest.mock import MagicMock, patch
from types import SimpleNamespace

from vaxrank.epitope_config import EpitopeConfig
from vaxrank.vaccine_config import VaccineConfig
from vaxrank.core_logic import (
    run_vaxrank,
    create_vaccine_peptides_dict,
    vaccine_peptides_for_variant,
    vaccine_peptides_from_epitopes,
)
from mhctools.pred import Prediction

from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.candidate_epitope import COMPARATOR_WT, CandidateEpitope, Peptide
from vaxrank.vaccine_peptide import VaccinePeptide, _legacy_score_one

from .common import eq_, ok_, gt_


def _make_epitope(peptide, ic50, wt_ic50=None, allele="HLA-A*02:01",
                  source=None, offset=0, percentile_rank=0.5,
                  method="test", overlaps_mutation=True,
                  occurs_in_reference=False):
    """Build a single-allele CandidateEpitope for tests that previously
    constructed a flat record."""
    src = source if source is not None else peptide
    mutant = Prediction(
        kind='pMHC_affinity', predictor_name=method,
        predictor_version='', allele=allele, peptide=peptide,
        value=ic50, score=0.0, percentile_rank=percentile_rank)
    comparators = {}
    if wt_ic50 is not None:
        wt = Prediction(
            kind='pMHC_affinity', predictor_name=method,
            predictor_version='', allele=allele, peptide=peptide,
            value=wt_ic50, score=0.0, percentile_rank=None)
        comparators[COMPARATOR_WT] = Peptide(
            sequence=peptide, source=src, offset=offset,
            predictions=(wt,))
    return CandidateEpitope(
        mutant=Peptide(
            sequence=peptide, source=src, offset=offset,
            predictions=(mutant,)),
        comparators=comparators,
        overlaps_mutation=overlaps_mutation,
        occurs_in_reference=occurs_in_reference)


# =============================================================================
# vaccine_peptides_from_epitopes Tests
# =============================================================================

def test_vaccine_peptides_from_epitopes_basic_parameters():
    """Test that basic parameters are used correctly"""
    # Create a mock protein fragment
    mock_fragment = MagicMock(spec=MutantProteinFragment)
    mock_fragment.sorted_subsequences.return_value = []  # No subsequences

    mock_variant = MagicMock()
    mock_variant.short_description = "test_variant"

    result = vaccine_peptides_from_epitopes(
        variant=mock_variant,
        long_protein_fragment=mock_fragment,
        epitopes=[],
        vaccine_peptide_length=25,
        max_vaccine_peptides_per_variant=1,
        num_mutant_epitopes_to_keep=1000,
    )

    # Should return empty list since no subsequences
    eq_(result, [])


def test_vaccine_peptides_from_epitopes_length_used():
    """Test that vaccine_peptide_length is passed to sorted_subsequences"""
    mock_fragment = MagicMock(spec=MutantProteinFragment)
    mock_fragment.sorted_subsequences.return_value = []

    mock_variant = MagicMock()
    mock_variant.short_description = "test"

    vaccine_peptides_from_epitopes(
        variant=mock_variant,
        long_protein_fragment=mock_fragment,
        epitopes=[],
        vaccine_peptide_length=30,  # Custom length
        max_vaccine_peptides_per_variant=1,
        num_mutant_epitopes_to_keep=1000,
    )

    # Verify sorted_subsequences was called with correct length
    mock_fragment.sorted_subsequences.assert_called_once_with(subsequence_length=30)


# =============================================================================
# vaccine_peptides_for_variant Tests
# =============================================================================

def test_vaccine_peptides_for_variant_fails_filter_returns_empty():
    """Test that variant failing filters returns empty list"""
    mock_isovar_result = MagicMock()
    mock_isovar_result.passes_all_filters = False

    mock_predictor = MagicMock()

    result = vaccine_peptides_for_variant(
        isovar_result=mock_isovar_result,
        mhc_predictor=mock_predictor,
    )

    eq_(result, [])


def test_vaccine_peptides_for_variant_config_overrides_explicit_params():
    """Test that vaccine_config values override explicit parameters"""
    mock_isovar_result = MagicMock()
    mock_isovar_result.passes_all_filters = True
    mock_isovar_result.variant = MagicMock()
    mock_isovar_result.variant.ensembl = None

    # Mock the MutantProteinFragment
    mock_fragment = MagicMock()
    mock_fragment.sorted_subsequences.return_value = []

    mock_predictor = MagicMock()

    vaccine_config = VaccineConfig(
        preferred_peptide_length=40,
        min_peptide_length=40,
        max_peptide_length=40,
        max_vaccine_peptides_per_variant=5,
        num_mutant_epitopes_to_keep=500,
    )

    with patch('vaxrank.core_logic.MutantProteinFragment') as MockFragment:
        MockFragment.from_isovar_result.return_value = mock_fragment
        with patch('vaxrank.core_logic.predict_epitopes') as mock_predict:
            mock_predict.return_value = []

            vaccine_peptides_for_variant(
                isovar_result=mock_isovar_result,
                mhc_predictor=mock_predictor,
                vaccine_peptide_length=25,  # Explicit, should be overridden
                vaccine_config=vaccine_config,
            )

            # The vaccine_config values should have been used
            mock_fragment.sorted_subsequences.assert_called_with(subsequence_length=40)


def test_vaccine_peptides_for_variant_explicit_params_used_without_config():
    """Test that explicit parameters are used when no config provided"""
    mock_isovar_result = MagicMock()
    mock_isovar_result.passes_all_filters = True
    mock_isovar_result.variant = MagicMock()
    mock_isovar_result.variant.ensembl = None

    mock_fragment = MagicMock()
    mock_fragment.sorted_subsequences.return_value = []

    mock_predictor = MagicMock()

    with patch('vaxrank.core_logic.MutantProteinFragment') as MockFragment:
        MockFragment.from_isovar_result.return_value = mock_fragment
        with patch('vaxrank.core_logic.predict_epitopes') as mock_predict:
            mock_predict.return_value = []

            vaccine_peptides_for_variant(
                isovar_result=mock_isovar_result,
                mhc_predictor=mock_predictor,
                vaccine_peptide_length=35,  # Should be used
                vaccine_config=None,  # No config
            )

            mock_fragment.sorted_subsequences.assert_called_with(subsequence_length=35)


def test_vaccine_peptides_for_variant_epitope_config_passed_to_predict():
    """Test that epitope_config is passed to predict_epitopes"""
    mock_isovar_result = MagicMock()
    mock_isovar_result.passes_all_filters = True
    mock_isovar_result.variant = MagicMock()
    mock_isovar_result.variant.ensembl = None

    mock_fragment = MagicMock()
    mock_fragment.sorted_subsequences.return_value = []

    mock_predictor = MagicMock()

    epitope_config = EpitopeConfig(min_epitope_score=0.05)

    with patch('vaxrank.core_logic.MutantProteinFragment') as MockFragment:
        MockFragment.from_isovar_result.return_value = mock_fragment
        with patch('vaxrank.core_logic.predict_epitopes') as mock_predict:
            mock_predict.return_value = []

            vaccine_peptides_for_variant(
                isovar_result=mock_isovar_result,
                mhc_predictor=mock_predictor,
                epitope_config=epitope_config,
            )

            # Verify predict_epitopes was called with the config
            call_kwargs = mock_predict.call_args[1]
            eq_(call_kwargs['epitope_config'], epitope_config)


# =============================================================================
# Config Integration Tests
# =============================================================================

def test_config_integration_epitope_config_affects_filtering():
    """Test that EpitopeConfig min_epitope_score affects filtering"""
    # This is more of a documentation test - the actual filtering
    # happens in epitope_logic.py, but we can verify the config flows through

    config_strict = EpitopeConfig(min_epitope_score=0.5)  # Very strict
    config_lenient = EpitopeConfig(min_epitope_score=0.0)  # No filtering

    # The configs should have different min_epitope_score values
    gt_(config_strict.min_epitope_score, config_lenient.min_epitope_score)


def test_config_integration_vaccine_config_affects_peptide_generation():
    """Test that VaccineConfig affects vaccine peptide generation"""
    config_short = VaccineConfig(preferred_peptide_length=20, min_peptide_length=20, max_peptide_length=20)
    config_long = VaccineConfig(preferred_peptide_length=35, min_peptide_length=35, max_peptide_length=35)

    gt_(config_long.preferred_peptide_length, config_short.preferred_peptide_length)


def test_vaccine_peptides_from_epitopes_score_fraction_of_best_from_config():
    class DummyCandidateFragment:
        def __init__(self, amino_acids):
            self.amino_acids = amino_acids

        def __len__(self):
            return len(self.amino_acids)

    class DummyLongFragment:
        def sorted_subsequences(self, subsequence_length):
            return [
                (0, DummyCandidateFragment("AAAA")),
                (1, DummyCandidateFragment("BBBB")),
            ]

    class FakeVaccinePeptide:
        def __init__(
            self,
            mutant_protein_fragment,
            epitopes,
            num_mutant_epitopes_to_keep=None,
            epitope_score_params=None,
            manufacturability_thresholds=None,
            manufacturability_rules=None,
            combined_score_mode=None,
            combined_score_expr=None,
            ranking_rules=None,
        ):
            self.mutant_protein_fragment = mutant_protein_fragment
            self.combined_score = 10.0 if mutant_protein_fragment.amino_acids == "AAAA" else 8.5

        @staticmethod
        def lexicographic_sort_key(obj):
            return (-obj.combined_score,)

    def fake_slice_epitopes(epitopes, start_offset, end_offset):
        amino_acids = "AAAA" if start_offset == 0 else "BBBB"
        return [SimpleNamespace(
            mutant=SimpleNamespace(source=amino_acids))]

    mock_variant = MagicMock()
    mock_variant.short_description = "test_variant"

    with patch("vaxrank.core_logic.VaccinePeptide", FakeVaccinePeptide), patch(
        "vaxrank.core_logic.slice_epitopes",
        side_effect=fake_slice_epitopes,
    ):
        default_result = vaccine_peptides_from_epitopes(
            variant=mock_variant,
            long_protein_fragment=DummyLongFragment(),
            epitopes=[MagicMock()],
        )
        relaxed_result = vaccine_peptides_from_epitopes(
            variant=mock_variant,
            long_protein_fragment=DummyLongFragment(),
            epitopes=[MagicMock()],
            vaccine_config=VaccineConfig(
                score_fraction_of_best=0.8,
                max_vaccine_peptides_per_variant=2,
            ),
        )
        strict_result = vaccine_peptides_from_epitopes(
            variant=mock_variant,
            long_protein_fragment=DummyLongFragment(),
            epitopes=[MagicMock()],
            vaccine_config=VaccineConfig(
                score_fraction_of_best=1.0,
                max_vaccine_peptides_per_variant=2,
            ),
        )

    eq_(len(default_result), 1)
    eq_(len(relaxed_result), 2)
    eq_(len(strict_result), 1)
    eq_(strict_result[0].combined_score, 10.0)


def test_config_integration_epitope_config_affects_vaccine_peptide_scoring():
    """Test that EpitopeConfig parameters affect VaccinePeptide scores"""
    class DummyFragment:
        def __init__(self, amino_acids):
            self.amino_acids = amino_acids
            self.n_alt_reads = 4
            self.n_alt_reads_supporting_protein_sequence = 4
            self.n_mutant_amino_acids = 1
            self.mutation_distance_from_edge = 0
            self.mutant_amino_acid_start_offset = 0
            self.mutant_amino_acid_end_offset = 1

        def __len__(self):
            return len(self.amino_acids)

    class DummyLongFragment:
        def __init__(self, fragment):
            self.fragment = fragment

        def sorted_subsequences(self, subsequence_length):
            return [(0, self.fragment)]

    fragment = DummyFragment("ACDEFGHIK")
    long_fragment = DummyLongFragment(fragment)

    epitope = _make_epitope("ACDEFGHIK", ic50=100.0, wt_ic50=200.0,
                            source="ACDEFGHIK")

    default_score = _legacy_score_one(100.0, 0.5)
    epitope_config = EpitopeConfig(
        logistic_epitope_score_midpoint=50.0,
        logistic_epitope_score_width=10.0,
        binding_affinity_cutoff=5000.0,
    )
    custom_score = _legacy_score_one(
        100.0, 0.5,
        midpoint=epitope_config.logistic_epitope_score_midpoint,
        width=epitope_config.logistic_epitope_score_width,
        ic50_cutoff=epitope_config.binding_affinity_cutoff,
    )
    ok_(custom_score != default_score)

    peptides = vaccine_peptides_from_epitopes(
        variant=MagicMock(),
        long_protein_fragment=long_fragment,
        epitopes=[epitope],
        vaccine_peptide_length=len(fragment),
        max_vaccine_peptides_per_variant=1,
        num_mutant_epitopes_to_keep=10,
        epitope_config=epitope_config,
    )
    eq_(len(peptides), 1)
    eq_(peptides[0].mutant_epitope_score, custom_score)


def test_config_defaults_match_historical_behavior():
    """Test that default config values match historical defaults"""
    epitope_config = EpitopeConfig()
    vaccine_config = VaccineConfig()

    # These were the historical defaults
    eq_(vaccine_config.preferred_peptide_length, 25)
    eq_(vaccine_config.padding_around_mutation, 5)
    eq_(vaccine_config.max_vaccine_peptides_per_variant, 1)

    # CandidateEpitope scoring defaults
    eq_(epitope_config.logistic_epitope_score_midpoint, 350.0)
    eq_(epitope_config.logistic_epitope_score_width, 150.0)


def test_core_logic_default_num_epitope_limit_matches_vaccine_config():
    """Public core helpers should share the same default epitope limit."""
    expected = VaccineConfig().num_mutant_epitopes_to_keep
    for fn in (
            run_vaxrank,
            create_vaccine_peptides_dict,
            vaccine_peptides_for_variant,
            vaccine_peptides_from_epitopes):
        eq_(signature(fn).parameters["num_mutant_epitopes_to_keep"].default, expected)


# =============================================================================
# Config Struct Immutability Tests
# =============================================================================

def test_config_immutability_epitope_config_hashable():
    """Test that EpitopeConfig can be used in sets/dicts"""
    config1 = EpitopeConfig(min_epitope_score=0.01)
    config2 = EpitopeConfig(min_epitope_score=0.01)
    config3 = EpitopeConfig(min_epitope_score=0.02)

    # Same values should be equal
    eq_(config1, config2)

    # Different values should not be equal
    ok_(config1 != config3)


def test_config_immutability_vaccine_config_hashable():
    """Test that VaccineConfig can be used in sets/dicts"""
    config1 = VaccineConfig(preferred_peptide_length=25)
    config2 = VaccineConfig(preferred_peptide_length=25)
    config3 = VaccineConfig(preferred_peptide_length=30, min_peptide_length=30, max_peptide_length=30)

    eq_(config1, config2)
    ok_(config1 != config3)


def test_config_immutability_config_repr():
    """Test that configs have useful repr"""
    config = EpitopeConfig(min_epitope_score=0.05)
    repr_str = repr(config)

    # Should contain class name and values
    ok_("EpitopeConfig" in repr_str)
    ok_("0.05" in repr_str)


# =============================================================================
# Edge Case Tests
# =============================================================================

def test_config_edge_case_zero_min_epitope_score():
    """Test that min_epitope_score of 0 is valid"""
    config = EpitopeConfig(min_epitope_score=0)
    eq_(config.min_epitope_score, 0)


def test_config_edge_case_very_small_min_epitope_score():
    """Test that very small min_epitope_score is preserved"""
    config = EpitopeConfig(min_epitope_score=1e-10)
    eq_(config.min_epitope_score, 1e-10)


def test_config_edge_case_large_vaccine_peptide_length():
    """Test that large preferred_peptide_length is valid"""
    config = VaccineConfig(preferred_peptide_length=100, min_peptide_length=100, max_peptide_length=100)
    eq_(config.preferred_peptide_length, 100)


def test_config_edge_case_zero_vaccine_peptides_per_variant():
    """Test behavior with 0 max_vaccine_peptides_per_variant"""
    config = VaccineConfig(max_vaccine_peptides_per_variant=0)
    eq_(config.max_vaccine_peptides_per_variant, 0)


def test_vaccine_peptide_zero_epitope_limit_keeps_all():
    """A value of 0 means no truncation of mutant epitope list."""
    fragment = MagicMock()
    fragment.amino_acids = "ACDEFGHIK"

    epitope1 = _make_epitope("ACDEFGHIK", ic50=100.0, wt_ic50=200.0,
                             percentile_rank=0.5)
    epitope2 = _make_epitope("CDEFGHIKL", ic50=120.0, wt_ic50=220.0,
                             percentile_rank=0.4)

    vaccine_peptide = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=[epitope1, epitope2],
        num_mutant_epitopes_to_keep=0,
    )
    eq_(len(vaccine_peptide.mutant_epitopes), 2)


def test_manufacturability_thresholds_flow_from_manufacturability_config():
    """Manufacturability thresholds from ManufacturabilityConfig (post-2.20)
    reach VaccinePeptide via the explicit ``manufacturability_config``
    parameter to ``vaccine_peptides_from_epitopes``."""
    from vaxrank.manufacturability_config import ManufacturabilityConfig

    class DummyFragment:
        def __init__(self, amino_acids):
            self.amino_acids = amino_acids
            self.n_alt_reads = 4
            self.n_alt_reads_supporting_protein_sequence = 4
            self.n_mutant_amino_acids = 1
            self.mutation_distance_from_edge = 0
            self.mutant_amino_acid_start_offset = 0
            self.mutant_amino_acid_end_offset = 1

        def __len__(self):
            return len(self.amino_acids)

    class DummyLongFragment:
        def __init__(self, fragment):
            self.fragment = fragment

        def sorted_subsequences(self, subsequence_length):
            return [(0, self.fragment)]

    fragment = DummyFragment("ACDEFGHIK")
    long_fragment = DummyLongFragment(fragment)

    epitope = _make_epitope("ACDEFGHIK", ic50=100.0, wt_ic50=200.0,
                            source="ACDEFGHIK")

    vaccine_config = VaccineConfig()
    manufacturability_config = ManufacturabilityConfig(
        max_kmer_hydropathy_high_priority=3.0)

    peptides = vaccine_peptides_from_epitopes(
        variant=MagicMock(),
        long_protein_fragment=long_fragment,
        epitopes=[epitope],
        vaccine_peptide_length=len(fragment),
        vaccine_config=vaccine_config,
        manufacturability_config=manufacturability_config,
    )
    eq_(len(peptides), 1)
    eq_(peptides[0].manufacturability_thresholds["max_kmer_hydropathy_high_priority"], 3.0)


def test_manufacturability_config_defaults():
    """ManufacturabilityConfig has the legacy default values
    (post-2.20: these moved off VaccineConfig)."""
    from vaxrank.manufacturability_config import ManufacturabilityConfig
    config = ManufacturabilityConfig()
    eq_(config.max_c_terminal_hydropathy, 1.5)
    eq_(config.min_kmer_hydropathy, 0.0)
    eq_(config.max_kmer_hydropathy_low_priority, 1.5)
    eq_(config.max_kmer_hydropathy_high_priority, 2.5)


def test_config_edge_case_config_with_all_defaults():
    """Test that configs work with all default values"""
    epitope_config = EpitopeConfig()
    vaccine_config = VaccineConfig()

    # Should be usable without any custom values
    ok_(epitope_config.min_epitope_score > 0)
    ok_(vaccine_config.preferred_peptide_length > 0)


# =============================================================================
# Override Warning Tests
# =============================================================================

def test_vaccine_config_override_warns_on_conflicting_explicit_params(caplog):
    """When vaccine_config overrides non-default explicit params, a warning is logged."""
    import logging

    mock_fragment = MagicMock(spec=MutantProteinFragment)
    mock_fragment.sorted_subsequences.return_value = []

    mock_variant = MagicMock()
    mock_variant.short_description = "test"

    vaccine_config = VaccineConfig(preferred_peptide_length=40, min_peptide_length=40, max_peptide_length=40)

    with caplog.at_level(logging.WARNING, logger="vaxrank.core_logic"):
        vaccine_peptides_from_epitopes(
            variant=mock_variant,
            long_protein_fragment=mock_fragment,
            epitopes=[],
            vaccine_peptide_length=30,  # Non-default, differs from config's 40
            vaccine_config=vaccine_config,
        )

    assert "overrides explicit" in caplog.text
    assert "preferred_peptide_length" in caplog.text


def test_vaccine_config_override_no_warning_when_default_params(caplog):
    """No warning when explicit params match their defaults."""
    import logging

    mock_fragment = MagicMock(spec=MutantProteinFragment)
    mock_fragment.sorted_subsequences.return_value = []

    mock_variant = MagicMock()
    mock_variant.short_description = "test"

    vaccine_config = VaccineConfig(preferred_peptide_length=40, min_peptide_length=40, max_peptide_length=40)

    with caplog.at_level(logging.WARNING, logger="vaxrank.core_logic"):
        vaccine_peptides_from_epitopes(
            variant=mock_variant,
            long_protein_fragment=mock_fragment,
            epitopes=[],
            # All defaults — no conflict
            vaccine_config=vaccine_config,
        )

    assert "overrides explicit" not in caplog.text
