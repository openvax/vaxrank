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
Tests for YAML config system including EpitopeConfig, VaccineConfig,
and the config_from_args functions.
"""

import argparse
import os
import tempfile

import pytest
import msgspec

from vaxrank.epitope_config import (
    EpitopeConfig,
    DEFAULT_MIN_EPITOPE_SCORE,
    DEFAULT_BINDING_AFFINITY_CUTOFF,
)
from vaxrank.vaccine_config import VaccineConfig
from vaxrank.cli.epitope_config_args import (
    epitope_config_from_args,
    add_epitope_prediction_args,
)
from vaxrank.cli.vaccine_config_args import (
    vaccine_config_from_args,
    add_vaccine_peptide_args,
)

from .common import eq_


# =============================================================================
# EpitopeConfig Tests
# =============================================================================

class TestEpitopeConfig:
    """Tests for EpitopeConfig msgspec.Struct"""

    def test_default_values(self):
        """Test that EpitopeConfig has correct default values"""
        config = EpitopeConfig()
        eq_(config.logistic_epitope_score_midpoint, 350.0)
        eq_(config.logistic_epitope_score_width, 150.0)
        eq_(config.min_epitope_score, DEFAULT_MIN_EPITOPE_SCORE)
        eq_(config.binding_affinity_cutoff, DEFAULT_BINDING_AFFINITY_CUTOFF)

    def test_custom_values(self):
        """Test creating EpitopeConfig with custom values"""
        config = EpitopeConfig(
            logistic_epitope_score_midpoint=400.0,
            logistic_epitope_score_width=200.0,
            min_epitope_score=0.01,
            binding_affinity_cutoff=1000.0,
        )
        eq_(config.logistic_epitope_score_midpoint, 400.0)
        eq_(config.logistic_epitope_score_width, 200.0)
        eq_(config.min_epitope_score, 0.01)
        eq_(config.binding_affinity_cutoff, 1000.0)

    def test_partial_custom_values(self):
        """Test creating EpitopeConfig with only some custom values"""
        config = EpitopeConfig(min_epitope_score=0.05)
        eq_(config.min_epitope_score, 0.05)
        # Other values should be defaults
        eq_(config.logistic_epitope_score_midpoint, 350.0)
        eq_(config.logistic_epitope_score_width, 150.0)
        eq_(config.binding_affinity_cutoff, DEFAULT_BINDING_AFFINITY_CUTOFF)

    def test_msgspec_encode_decode(self):
        """Test that EpitopeConfig can be encoded/decoded with msgspec"""
        config = EpitopeConfig(
            logistic_epitope_score_midpoint=400.0,
            min_epitope_score=0.01,
        )
        # Encode to JSON
        json_bytes = msgspec.json.encode(config)
        # Decode back
        decoded = msgspec.json.decode(json_bytes, type=EpitopeConfig)
        eq_(decoded.logistic_epitope_score_midpoint, 400.0)
        eq_(decoded.min_epitope_score, 0.01)

    def test_yaml_decode(self):
        """Test decoding EpitopeConfig from YAML"""
        yaml_content = """
logistic_epitope_score_midpoint: 500.0
logistic_epitope_score_width: 100.0
min_epitope_score: 0.001
binding_affinity_cutoff: 2000.0
"""
        config = msgspec.yaml.decode(yaml_content, type=EpitopeConfig)
        eq_(config.logistic_epitope_score_midpoint, 500.0)
        eq_(config.logistic_epitope_score_width, 100.0)
        eq_(config.min_epitope_score, 0.001)
        eq_(config.binding_affinity_cutoff, 2000.0)

    def test_yaml_decode_partial(self):
        """Test decoding EpitopeConfig from YAML with partial values"""
        yaml_content = """
min_epitope_score: 0.05
"""
        config = msgspec.yaml.decode(yaml_content, type=EpitopeConfig)
        eq_(config.min_epitope_score, 0.05)
        # Defaults for unspecified values
        eq_(config.logistic_epitope_score_midpoint, 350.0)

    def test_convert_from_dict(self):
        """Test converting a dict to EpitopeConfig"""
        config_dict = {
            "logistic_epitope_score_midpoint": 450.0,
            "min_epitope_score": 0.02,
        }
        config = msgspec.convert(config_dict, EpitopeConfig)
        eq_(config.logistic_epitope_score_midpoint, 450.0)
        eq_(config.min_epitope_score, 0.02)
        # Defaults for unspecified
        eq_(config.logistic_epitope_score_width, 150.0)


# =============================================================================
# VaccineConfig Tests
# =============================================================================

class TestVaccineConfig:
    """Tests for VaccineConfig msgspec.Struct"""

    def test_default_values(self):
        """Test that VaccineConfig has correct default values"""
        config = VaccineConfig()
        eq_(config.vaccine_peptide_length, 25)
        eq_(config.padding_around_mutation, 5)
        eq_(config.max_vaccine_peptides_per_variant, 1)
        eq_(config.num_mutant_epitopes_to_keep, 1000)

    def test_custom_values(self):
        """Test creating VaccineConfig with custom values"""
        config = VaccineConfig(
            vaccine_peptide_length=30,
            padding_around_mutation=10,
            max_vaccine_peptides_per_variant=3,
            num_mutant_epitopes_to_keep=500,
        )
        eq_(config.vaccine_peptide_length, 30)
        eq_(config.padding_around_mutation, 10)
        eq_(config.max_vaccine_peptides_per_variant, 3)
        eq_(config.num_mutant_epitopes_to_keep, 500)

    def test_partial_custom_values(self):
        """Test creating VaccineConfig with only some custom values"""
        config = VaccineConfig(vaccine_peptide_length=35)
        eq_(config.vaccine_peptide_length, 35)
        # Other values should be defaults
        eq_(config.padding_around_mutation, 5)
        eq_(config.max_vaccine_peptides_per_variant, 1)

    def test_msgspec_encode_decode(self):
        """Test that VaccineConfig can be encoded/decoded with msgspec"""
        config = VaccineConfig(
            vaccine_peptide_length=30,
            max_vaccine_peptides_per_variant=5,
        )
        json_bytes = msgspec.json.encode(config)
        decoded = msgspec.json.decode(json_bytes, type=VaccineConfig)
        eq_(decoded.vaccine_peptide_length, 30)
        eq_(decoded.max_vaccine_peptides_per_variant, 5)

    def test_yaml_decode(self):
        """Test decoding VaccineConfig from YAML"""
        yaml_content = """
vaccine_peptide_length: 40
padding_around_mutation: 8
max_vaccine_peptides_per_variant: 2
num_mutant_epitopes_to_keep: 2000
"""
        config = msgspec.yaml.decode(yaml_content, type=VaccineConfig)
        eq_(config.vaccine_peptide_length, 40)
        eq_(config.padding_around_mutation, 8)
        eq_(config.max_vaccine_peptides_per_variant, 2)
        eq_(config.num_mutant_epitopes_to_keep, 2000)

    def test_yaml_decode_partial(self):
        """Test decoding VaccineConfig from YAML with partial values"""
        yaml_content = """
vaccine_peptide_length: 35
"""
        config = msgspec.yaml.decode(yaml_content, type=VaccineConfig)
        eq_(config.vaccine_peptide_length, 35)
        eq_(config.padding_around_mutation, 5)  # default

    def test_convert_from_dict(self):
        """Test converting a dict to VaccineConfig"""
        config_dict = {
            "vaccine_peptide_length": 28,
            "max_vaccine_peptides_per_variant": 4,
        }
        config = msgspec.convert(config_dict, VaccineConfig)
        eq_(config.vaccine_peptide_length, 28)
        eq_(config.max_vaccine_peptides_per_variant, 4)


# =============================================================================
# epitope_config_from_args Tests
# =============================================================================

class TestEpitopeConfigFromArgs:
    """Tests for epitope_config_from_args function"""

    def test_no_config_file_no_cli_args(self):
        """Test with no config file and no CLI overrides - should use defaults"""
        args = argparse.Namespace(config=None, min_epitope_score=None)
        config = epitope_config_from_args(args)
        eq_(config.min_epitope_score, DEFAULT_MIN_EPITOPE_SCORE)
        eq_(config.logistic_epitope_score_midpoint, 350.0)

    def test_cli_override_only(self):
        """Test CLI argument overrides without config file"""
        args = argparse.Namespace(config=None, min_epitope_score=0.05)
        config = epitope_config_from_args(args)
        eq_(config.min_epitope_score, 0.05)

    def test_yaml_config_file(self):
        """Test loading config from YAML file"""
        yaml_content = """
epitope_config:
  min_epitope_score: 0.01
  logistic_epitope_score_midpoint: 400.0
  binding_affinity_cutoff: 1000.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            args = argparse.Namespace(config=config_path, min_epitope_score=None)
            config = epitope_config_from_args(args)
            eq_(config.min_epitope_score, 0.01)
            eq_(config.logistic_epitope_score_midpoint, 400.0)
            eq_(config.binding_affinity_cutoff, 1000.0)
        finally:
            os.unlink(config_path)

    def test_cli_overrides_yaml(self):
        """Test that CLI arguments override YAML config values"""
        yaml_content = """
epitope_config:
  min_epitope_score: 0.01
  logistic_epitope_score_midpoint: 400.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            # CLI sets min_epitope_score to 0.05, should override YAML's 0.01
            args = argparse.Namespace(config=config_path, min_epitope_score=0.05)
            config = epitope_config_from_args(args)
            eq_(config.min_epitope_score, 0.05)  # CLI override
            eq_(config.logistic_epitope_score_midpoint, 400.0)  # from YAML
        finally:
            os.unlink(config_path)

    def test_yaml_without_epitope_config_section(self):
        """Test YAML file without epitope_config section uses defaults"""
        yaml_content = """
vaccine_config:
  vaccine_peptide_length: 30
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            args = argparse.Namespace(config=config_path, min_epitope_score=None)
            config = epitope_config_from_args(args)
            # Should use defaults since epitope_config section is missing
            eq_(config.min_epitope_score, DEFAULT_MIN_EPITOPE_SCORE)
        finally:
            os.unlink(config_path)

    def test_empty_yaml_file(self):
        """Test with empty YAML file uses defaults"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("")
            config_path = f.name

        try:
            args = argparse.Namespace(config=config_path, min_epitope_score=None)
            config = epitope_config_from_args(args)
            eq_(config.min_epitope_score, DEFAULT_MIN_EPITOPE_SCORE)
        finally:
            os.unlink(config_path)


# =============================================================================
# vaccine_config_from_args Tests
# =============================================================================

class TestVaccineConfigFromArgs:
    """Tests for vaccine_config_from_args function"""

    def test_no_config_file_no_cli_args(self):
        """Test with no config file and no CLI overrides - should use defaults"""
        args = argparse.Namespace(
            config=None,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_mutation=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        config = vaccine_config_from_args(args)
        eq_(config.vaccine_peptide_length, 25)
        eq_(config.padding_around_mutation, 5)
        eq_(config.max_vaccine_peptides_per_variant, 1)
        eq_(config.num_mutant_epitopes_to_keep, 1000)

    def test_cli_override_vaccine_peptide_length(self):
        """Test CLI override for vaccine_peptide_length"""
        args = argparse.Namespace(
            config=None,
            vaccine_peptide_length=35,
            padding_around_mutation=None,
            max_vaccine_peptides_per_mutation=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        config = vaccine_config_from_args(args)
        eq_(config.vaccine_peptide_length, 35)

    def test_cli_override_all_values(self):
        """Test CLI override for all vaccine config values"""
        args = argparse.Namespace(
            config=None,
            vaccine_peptide_length=30,
            padding_around_mutation=10,
            max_vaccine_peptides_per_mutation=5,
            num_epitopes_per_vaccine_peptide=500,
        )
        config = vaccine_config_from_args(args)
        eq_(config.vaccine_peptide_length, 30)
        eq_(config.padding_around_mutation, 10)
        eq_(config.max_vaccine_peptides_per_variant, 5)
        eq_(config.num_mutant_epitopes_to_keep, 500)

    def test_yaml_config_file(self):
        """Test loading vaccine config from YAML file"""
        yaml_content = """
vaccine_config:
  vaccine_peptide_length: 40
  padding_around_mutation: 8
  max_vaccine_peptides_per_variant: 3
  num_mutant_epitopes_to_keep: 2000
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            args = argparse.Namespace(
                config=config_path,
                vaccine_peptide_length=None,
                padding_around_mutation=None,
                max_vaccine_peptides_per_mutation=None,
                num_epitopes_per_vaccine_peptide=None,
            )
            config = vaccine_config_from_args(args)
            eq_(config.vaccine_peptide_length, 40)
            eq_(config.padding_around_mutation, 8)
            eq_(config.max_vaccine_peptides_per_variant, 3)
            eq_(config.num_mutant_epitopes_to_keep, 2000)
        finally:
            os.unlink(config_path)

    def test_cli_overrides_yaml(self):
        """Test that CLI arguments override YAML config values"""
        yaml_content = """
vaccine_config:
  vaccine_peptide_length: 40
  max_vaccine_peptides_per_variant: 3
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            args = argparse.Namespace(
                config=config_path,
                vaccine_peptide_length=50,  # Override YAML's 40
                padding_around_mutation=None,
                max_vaccine_peptides_per_mutation=None,
                num_epitopes_per_vaccine_peptide=None,
            )
            config = vaccine_config_from_args(args)
            eq_(config.vaccine_peptide_length, 50)  # CLI override
            eq_(config.max_vaccine_peptides_per_variant, 3)  # from YAML
        finally:
            os.unlink(config_path)

    def test_yaml_without_vaccine_config_section(self):
        """Test YAML file without vaccine_config section uses defaults"""
        yaml_content = """
epitope_config:
  min_epitope_score: 0.01
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            args = argparse.Namespace(
                config=config_path,
                vaccine_peptide_length=None,
                padding_around_mutation=None,
                max_vaccine_peptides_per_mutation=None,
                num_epitopes_per_vaccine_peptide=None,
            )
            config = vaccine_config_from_args(args)
            eq_(config.vaccine_peptide_length, 25)  # default
        finally:
            os.unlink(config_path)


# =============================================================================
# Combined Config Tests
# =============================================================================

class TestCombinedConfig:
    """Tests for loading both configs from a single YAML file"""

    def test_both_configs_from_single_yaml(self):
        """Test loading both epitope and vaccine configs from one YAML file"""
        yaml_content = """
epitope_config:
  min_epitope_score: 0.02
  logistic_epitope_score_midpoint: 400.0

vaccine_config:
  vaccine_peptide_length: 30
  max_vaccine_peptides_per_variant: 2
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            # Test epitope config
            epitope_args = argparse.Namespace(config=config_path, min_epitope_score=None)
            epitope_config = epitope_config_from_args(epitope_args)
            eq_(epitope_config.min_epitope_score, 0.02)
            eq_(epitope_config.logistic_epitope_score_midpoint, 400.0)

            # Test vaccine config
            vaccine_args = argparse.Namespace(
                config=config_path,
                vaccine_peptide_length=None,
                padding_around_mutation=None,
                max_vaccine_peptides_per_mutation=None,
                num_epitopes_per_vaccine_peptide=None,
            )
            vaccine_config = vaccine_config_from_args(vaccine_args)
            eq_(vaccine_config.vaccine_peptide_length, 30)
            eq_(vaccine_config.max_vaccine_peptides_per_variant, 2)
        finally:
            os.unlink(config_path)

    def test_mixed_cli_and_yaml_overrides(self):
        """Test mixing CLI and YAML overrides for both configs"""
        yaml_content = """
epitope_config:
  min_epitope_score: 0.01

vaccine_config:
  vaccine_peptide_length: 30
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            # Override epitope from CLI, keep vaccine from YAML
            epitope_args = argparse.Namespace(config=config_path, min_epitope_score=0.05)
            epitope_config = epitope_config_from_args(epitope_args)
            eq_(epitope_config.min_epitope_score, 0.05)  # CLI

            # Override vaccine from CLI, epitope from YAML is ignored
            vaccine_args = argparse.Namespace(
                config=config_path,
                vaccine_peptide_length=40,  # Override
                padding_around_mutation=None,
                max_vaccine_peptides_per_mutation=None,
                num_epitopes_per_vaccine_peptide=None,
            )
            vaccine_config = vaccine_config_from_args(vaccine_args)
            eq_(vaccine_config.vaccine_peptide_length, 40)  # CLI
        finally:
            os.unlink(config_path)


# =============================================================================
# CLI Argument Parser Tests
# =============================================================================

class TestCLIArgumentParsing:
    """Tests for CLI argument parsing with config options"""

    def test_add_epitope_prediction_args(self):
        """Test that add_epitope_prediction_args adds the expected arguments"""
        parser = argparse.ArgumentParser()
        add_epitope_prediction_args(parser)

        # Parse with min-epitope-score
        args = parser.parse_args(['--min-epitope-score', '0.05'])
        eq_(args.min_epitope_score, 0.05)

    def test_add_vaccine_peptide_args(self):
        """Test that add_vaccine_peptide_args adds the expected arguments"""
        parser = argparse.ArgumentParser()
        add_vaccine_peptide_args(parser)

        # Parse with various vaccine args
        args = parser.parse_args([
            '--vaccine-peptide-length', '30',
            '--padding-around-mutation', '10',
            '--max-vaccine-peptides-per-mutation', '5',
        ])
        eq_(args.vaccine_peptide_length, 30)
        eq_(args.padding_around_mutation, 10)
        eq_(args.max_vaccine_peptides_per_mutation, 5)

    def test_vaccine_peptide_args_defaults(self):
        """Test that vaccine peptide args default to None for config override"""
        parser = argparse.ArgumentParser()
        add_vaccine_peptide_args(parser)

        args = parser.parse_args([])
        eq_(args.vaccine_peptide_length, None)
        eq_(args.padding_around_mutation, None)
        eq_(args.max_vaccine_peptides_per_mutation, None)


# =============================================================================
# Error Handling Tests
# =============================================================================

class TestConfigErrorHandling:
    """Tests for error handling in config loading"""

    def test_nonexistent_config_file(self):
        """Test that nonexistent config file raises appropriate error"""
        args = argparse.Namespace(
            config='/nonexistent/path/config.yaml',
            min_epitope_score=None,
        )
        with pytest.raises(FileNotFoundError):
            epitope_config_from_args(args)

    def test_invalid_yaml_syntax(self):
        """Test that invalid YAML syntax raises appropriate error"""
        yaml_content = """
epitope_config:
  min_epitope_score: [invalid yaml
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            args = argparse.Namespace(config=config_path, min_epitope_score=None)
            with pytest.raises(msgspec.DecodeError):
                epitope_config_from_args(args)
        finally:
            os.unlink(config_path)

    def test_invalid_config_value_type(self):
        """Test that invalid value type raises appropriate error"""
        yaml_content = """
epitope_config:
  min_epitope_score: "not_a_number"
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            config_path = f.name

        try:
            args = argparse.Namespace(config=config_path, min_epitope_score=None)
            with pytest.raises((msgspec.ValidationError, TypeError)):
                epitope_config_from_args(args)
        finally:
            os.unlink(config_path)
