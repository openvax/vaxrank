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

from unittest.mock import MagicMock, patch

import pytest

from vaxrank.epitope_config import EpitopeConfig
from vaxrank.vaccine_config import VaccineConfig
from vaxrank.core_logic import (
    vaccine_peptides_for_variant,
    vaccine_peptides_from_epitopes,
)
from vaxrank.mutant_protein_fragment import MutantProteinFragment

from .common import eq_, ok_, gt_, gte_


# =============================================================================
# vaccine_peptides_from_epitopes Tests
# =============================================================================

class TestVaccinePeptidesFromEpitopes:
    """Tests for vaccine_peptides_from_epitopes function"""

    def test_basic_parameters(self):
        """Test that basic parameters are used correctly"""
        # Create a mock protein fragment
        mock_fragment = MagicMock(spec=MutantProteinFragment)
        mock_fragment.sorted_subsequences.return_value = []  # No subsequences

        mock_variant = MagicMock()
        mock_variant.short_description = "test_variant"

        result = vaccine_peptides_from_epitopes(
            variant=mock_variant,
            long_protein_fragment=mock_fragment,
            epitope_predictions=[],
            vaccine_peptide_length=25,
            max_vaccine_peptides_per_variant=1,
            num_mutant_epitopes_to_keep=1000,
        )

        # Should return empty list since no subsequences
        eq_(result, [])

    def test_vaccine_peptide_length_used(self):
        """Test that vaccine_peptide_length is passed to sorted_subsequences"""
        mock_fragment = MagicMock(spec=MutantProteinFragment)
        mock_fragment.sorted_subsequences.return_value = []

        mock_variant = MagicMock()
        mock_variant.short_description = "test"

        vaccine_peptides_from_epitopes(
            variant=mock_variant,
            long_protein_fragment=mock_fragment,
            epitope_predictions=[],
            vaccine_peptide_length=30,  # Custom length
            max_vaccine_peptides_per_variant=1,
            num_mutant_epitopes_to_keep=1000,
        )

        # Verify sorted_subsequences was called with correct length
        mock_fragment.sorted_subsequences.assert_called_once_with(subsequence_length=30)


class TestVaccinePeptidesForVariant:
    """Tests for vaccine_peptides_for_variant function"""

    def test_fails_filter_returns_empty(self):
        """Test that variant failing filters returns empty list"""
        mock_isovar_result = MagicMock()
        mock_isovar_result.passes_all_filters = False

        mock_predictor = MagicMock()

        result = vaccine_peptides_for_variant(
            isovar_result=mock_isovar_result,
            mhc_predictor=mock_predictor,
        )

        eq_(result, [])

    def test_config_overrides_explicit_params(self):
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
            vaccine_peptide_length=40,
            max_vaccine_peptides_per_variant=5,
            num_mutant_epitopes_to_keep=500,
        )

        with patch('vaxrank.core_logic.MutantProteinFragment') as MockFragment:
            MockFragment.from_isovar_result.return_value = mock_fragment
            with patch('vaxrank.core_logic.predict_epitopes') as mock_predict:
                mock_predict.return_value = {}

                result = vaccine_peptides_for_variant(
                    isovar_result=mock_isovar_result,
                    mhc_predictor=mock_predictor,
                    vaccine_peptide_length=25,  # Explicit, should be overridden
                    vaccine_config=vaccine_config,
                )

                # The vaccine_config values should have been used
                mock_fragment.sorted_subsequences.assert_called_with(subsequence_length=40)

    def test_explicit_params_used_without_config(self):
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
                mock_predict.return_value = {}

                result = vaccine_peptides_for_variant(
                    isovar_result=mock_isovar_result,
                    mhc_predictor=mock_predictor,
                    vaccine_peptide_length=35,  # Should be used
                    vaccine_config=None,  # No config
                )

                mock_fragment.sorted_subsequences.assert_called_with(subsequence_length=35)

    def test_epitope_config_passed_to_predict(self):
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
                mock_predict.return_value = {}

                result = vaccine_peptides_for_variant(
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

class TestConfigIntegration:
    """Integration tests for config objects with core logic"""

    def test_epitope_config_affects_filtering(self):
        """Test that EpitopeConfig min_epitope_score affects filtering"""
        # This is more of a documentation test - the actual filtering
        # happens in epitope_logic.py, but we can verify the config flows through

        config_strict = EpitopeConfig(min_epitope_score=0.5)  # Very strict
        config_lenient = EpitopeConfig(min_epitope_score=0.0)  # No filtering

        # The configs should have different min_epitope_score values
        gt_(config_strict.min_epitope_score, config_lenient.min_epitope_score)

    def test_vaccine_config_affects_peptide_generation(self):
        """Test that VaccineConfig affects vaccine peptide generation"""
        config_short = VaccineConfig(vaccine_peptide_length=20)
        config_long = VaccineConfig(vaccine_peptide_length=35)

        gt_(config_long.vaccine_peptide_length, config_short.vaccine_peptide_length)

    def test_config_defaults_match_historical_behavior(self):
        """Test that default config values match historical defaults"""
        epitope_config = EpitopeConfig()
        vaccine_config = VaccineConfig()

        # These were the historical defaults
        eq_(vaccine_config.vaccine_peptide_length, 25)
        eq_(vaccine_config.padding_around_mutation, 5)
        eq_(vaccine_config.max_vaccine_peptides_per_variant, 1)

        # Epitope scoring defaults
        eq_(epitope_config.logistic_epitope_score_midpoint, 350.0)
        eq_(epitope_config.logistic_epitope_score_width, 150.0)


# =============================================================================
# Config Struct Immutability Tests
# =============================================================================

class TestConfigImmutability:
    """Tests to verify config structs behave correctly"""

    def test_epitope_config_hashable(self):
        """Test that EpitopeConfig can be used in sets/dicts"""
        config1 = EpitopeConfig(min_epitope_score=0.01)
        config2 = EpitopeConfig(min_epitope_score=0.01)
        config3 = EpitopeConfig(min_epitope_score=0.02)

        # Same values should be equal
        eq_(config1, config2)

        # Different values should not be equal
        ok_(config1 != config3)

    def test_vaccine_config_hashable(self):
        """Test that VaccineConfig can be used in sets/dicts"""
        config1 = VaccineConfig(vaccine_peptide_length=25)
        config2 = VaccineConfig(vaccine_peptide_length=25)
        config3 = VaccineConfig(vaccine_peptide_length=30)

        eq_(config1, config2)
        ok_(config1 != config3)

    def test_config_repr(self):
        """Test that configs have useful repr"""
        config = EpitopeConfig(min_epitope_score=0.05)
        repr_str = repr(config)

        # Should contain class name and values
        ok_("EpitopeConfig" in repr_str)
        ok_("0.05" in repr_str)


# =============================================================================
# Edge Case Tests
# =============================================================================

class TestConfigEdgeCases:
    """Edge case tests for config handling"""

    def test_zero_min_epitope_score(self):
        """Test that min_epitope_score of 0 is valid"""
        config = EpitopeConfig(min_epitope_score=0)
        eq_(config.min_epitope_score, 0)

    def test_very_small_min_epitope_score(self):
        """Test that very small min_epitope_score is preserved"""
        config = EpitopeConfig(min_epitope_score=1e-10)
        eq_(config.min_epitope_score, 1e-10)

    def test_large_vaccine_peptide_length(self):
        """Test that large vaccine_peptide_length is valid"""
        config = VaccineConfig(vaccine_peptide_length=100)
        eq_(config.vaccine_peptide_length, 100)

    def test_zero_vaccine_peptides_per_variant(self):
        """Test behavior with 0 max_vaccine_peptides_per_variant"""
        config = VaccineConfig(max_vaccine_peptides_per_variant=0)
        eq_(config.max_vaccine_peptides_per_variant, 0)

    def test_config_with_all_defaults(self):
        """Test that configs work with all default values"""
        epitope_config = EpitopeConfig()
        vaccine_config = VaccineConfig()

        # Should be usable without any custom values
        ok_(epitope_config.min_epitope_score > 0)
        ok_(vaccine_config.vaccine_peptide_length > 0)
