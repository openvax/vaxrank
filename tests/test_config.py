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

import msgspec
import pytest

from vaxrank.epitope_config import (
    EpitopeConfig,
    DEFAULT_MIN_EPITOPE_SCORE,
    DEFAULT_BINDING_AFFINITY_CUTOFF,
)
from vaxrank.config.defaults import (
    DEFAULT_PERCENTILE_RANK_CUTOFF,
    DEFAULT_VACCINE_PEPTIDE_SCORE_FRACTION_OF_BEST,
)
from vaxrank.config.loader import load_vaxrank_config
from vaxrank.vaccine_config import VaccineConfig
from vaxrank.cli.arg_parser import add_config_override_args
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

def test_epitope_config_default_values():
    """Test that EpitopeConfig has correct default values"""
    config = EpitopeConfig()
    eq_(config.logistic_epitope_score_midpoint, 350.0)
    eq_(config.logistic_epitope_score_width, 150.0)
    eq_(config.min_epitope_score, DEFAULT_MIN_EPITOPE_SCORE)
    eq_(config.binding_affinity_cutoff, DEFAULT_BINDING_AFFINITY_CUTOFF)
    eq_(config.percentile_rank_cutoff, DEFAULT_PERCENTILE_RANK_CUTOFF)


def test_epitope_config_custom_values():
    """Test creating EpitopeConfig with custom values"""
    config = EpitopeConfig(
        logistic_epitope_score_midpoint=400.0,
        logistic_epitope_score_width=200.0,
        min_epitope_score=0.01,
        binding_affinity_cutoff=1000.0,
        percentile_rank_cutoff=5.0,
    )
    eq_(config.logistic_epitope_score_midpoint, 400.0)
    eq_(config.logistic_epitope_score_width, 200.0)
    eq_(config.min_epitope_score, 0.01)
    eq_(config.binding_affinity_cutoff, 1000.0)
    eq_(config.percentile_rank_cutoff, 5.0)


def test_epitope_config_partial_custom_values():
    """Test creating EpitopeConfig with only some custom values"""
    config = EpitopeConfig(min_epitope_score=0.05)
    eq_(config.min_epitope_score, 0.05)
    # Other values should be defaults
    eq_(config.logistic_epitope_score_midpoint, 350.0)
    eq_(config.logistic_epitope_score_width, 150.0)
    eq_(config.binding_affinity_cutoff, DEFAULT_BINDING_AFFINITY_CUTOFF)


def test_epitope_config_msgspec_encode_decode():
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


def test_epitope_config_yaml_decode():
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


def test_epitope_config_yaml_decode_partial():
    """Test decoding EpitopeConfig from YAML with partial values"""
    yaml_content = """
min_epitope_score: 0.05
"""
    config = msgspec.yaml.decode(yaml_content, type=EpitopeConfig)
    eq_(config.min_epitope_score, 0.05)
    # Defaults for unspecified values
    eq_(config.logistic_epitope_score_midpoint, 350.0)


def test_epitope_config_convert_from_dict():
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

def test_vaccine_config_default_values():
    """Test that VaccineConfig has correct default values"""
    config = VaccineConfig()
    eq_(config.vaccine_peptide_length, 25)
    eq_(config.padding_around_mutation, 5)
    eq_(config.max_vaccine_peptides_per_variant, 1)
    eq_(config.num_mutant_epitopes_to_keep, 1000)
    eq_(config.score_fraction_of_best, DEFAULT_VACCINE_PEPTIDE_SCORE_FRACTION_OF_BEST)


def test_vaccine_config_custom_values():
    """Test creating VaccineConfig with custom values"""
    config = VaccineConfig(
        vaccine_peptide_length=30,
        padding_around_mutation=10,
        max_vaccine_peptides_per_variant=3,
        num_mutant_epitopes_to_keep=500,
        score_fraction_of_best=0.95,
    )
    eq_(config.vaccine_peptide_length, 30)
    eq_(config.padding_around_mutation, 10)
    eq_(config.max_vaccine_peptides_per_variant, 3)
    eq_(config.num_mutant_epitopes_to_keep, 500)
    eq_(config.score_fraction_of_best, 0.95)


def test_vaccine_config_partial_custom_values():
    """Test creating VaccineConfig with only some custom values"""
    config = VaccineConfig(vaccine_peptide_length=35)
    eq_(config.vaccine_peptide_length, 35)
    # Other values should be defaults
    eq_(config.padding_around_mutation, 5)
    eq_(config.max_vaccine_peptides_per_variant, 1)


def test_vaccine_config_msgspec_encode_decode():
    """Test that VaccineConfig can be encoded/decoded with msgspec"""
    config = VaccineConfig(
        vaccine_peptide_length=30,
        max_vaccine_peptides_per_variant=5,
    )
    json_bytes = msgspec.json.encode(config)
    decoded = msgspec.json.decode(json_bytes, type=VaccineConfig)
    eq_(decoded.vaccine_peptide_length, 30)
    eq_(decoded.max_vaccine_peptides_per_variant, 5)


def test_vaccine_config_yaml_decode():
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


def test_vaccine_config_yaml_decode_partial():
    """Test decoding VaccineConfig from YAML with partial values"""
    yaml_content = """
vaccine_peptide_length: 35
"""
    config = msgspec.yaml.decode(yaml_content, type=VaccineConfig)
    eq_(config.vaccine_peptide_length, 35)
    eq_(config.padding_around_mutation, 5)  # default


def test_vaccine_config_convert_from_dict():
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

def test_epitope_config_from_args_no_config_file_no_cli_args():
    """Test with no config file and no CLI overrides - should use defaults"""
    args = argparse.Namespace(config=None, min_epitope_score=None)
    config = epitope_config_from_args(args)
    eq_(config.min_epitope_score, DEFAULT_MIN_EPITOPE_SCORE)
    eq_(config.logistic_epitope_score_midpoint, 350.0)


def test_epitope_config_from_args_cli_override_only():
    """Test CLI argument overrides without config file"""
    args = argparse.Namespace(config=None, min_epitope_score=0.05)
    config = epitope_config_from_args(args)
    eq_(config.min_epitope_score, 0.05)


def test_epitope_config_from_args_yaml_config_file():
    """Test loading config from YAML file"""
    yaml_content = """
epitopes:
  filters:
    min_score: 0.01
  scoring:
    derived_fields:
      affinity_score:
        midpoint: 400.0
        cutoff: 1000.0
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


def test_epitope_config_from_args_cli_overrides_yaml():
    """Test that CLI arguments override YAML config values"""
    yaml_content = """
epitopes:
  filters:
    min_score: 0.01
  scoring:
    derived_fields:
      affinity_score:
        midpoint: 400.0
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


def test_epitope_config_from_args_yaml_without_epitopes_section():
    """Test YAML file without epitopes section uses defaults"""
    yaml_content = """
vaccine_peptides:
  generation:
    lengths: 30
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(config=config_path, min_epitope_score=None)
        config = epitope_config_from_args(args)
        # Should use defaults since epitopes section is missing
        eq_(config.min_epitope_score, DEFAULT_MIN_EPITOPE_SCORE)
    finally:
        os.unlink(config_path)


def test_epitope_config_from_args_empty_yaml_file():
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


def test_epitope_config_from_args_nested_yaml_config_file():
    yaml_content = """
epitopes:
  filters:
    min_score: 0.02
  scoring:
    mode: percentile_rank
    derived_fields:
      affinity_score:
        midpoint: 410.0
        width: 110.0
        cutoff: 1200.0
      percentile_score:
        worst: 7.5
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            min_epitope_score=None,
            scoring_mode=None,
            config_set_overrides=None,
            config_expr_overrides=None,
        )
        config = epitope_config_from_args(args)
        eq_(config.min_epitope_score, 0.02)
        eq_(config.scoring_mode, "percentile_rank")
        eq_(config.logistic_epitope_score_midpoint, 410.0)
        eq_(config.logistic_epitope_score_width, 110.0)
        eq_(config.binding_affinity_cutoff, 1200.0)
        eq_(config.percentile_rank_cutoff, 7.5)
    finally:
        os.unlink(config_path)


# =============================================================================
# vaccine_config_from_args Tests
# =============================================================================

def test_vaccine_config_from_args_no_config_file_no_cli_args():
    """Test with no config file and no CLI overrides - should use defaults"""
    args = argparse.Namespace(
        config=None,
        vaccine_peptide_length=None,
        padding_around_mutation=None,
        max_vaccine_peptides_per_variant=None,
        num_epitopes_per_vaccine_peptide=None,
    )
    config = vaccine_config_from_args(args)
    eq_(config.vaccine_peptide_length, 25)
    eq_(config.padding_around_mutation, 5)
    eq_(config.max_vaccine_peptides_per_variant, 1)
    eq_(config.num_mutant_epitopes_to_keep, 1000)


def test_vaccine_config_from_args_cli_override_vaccine_peptide_length():
    """Test CLI override for vaccine_peptide_length"""
    args = argparse.Namespace(
        config=None,
        vaccine_peptide_length=35,
        padding_around_mutation=None,
        max_vaccine_peptides_per_variant=None,
        num_epitopes_per_vaccine_peptide=None,
    )
    config = vaccine_config_from_args(args)
    eq_(config.vaccine_peptide_length, 35)


def test_vaccine_config_from_args_cli_override_all_values():
    """Test CLI override for all vaccine config values"""
    args = argparse.Namespace(
        config=None,
        vaccine_peptide_length=30,
        padding_around_mutation=10,
        max_vaccine_peptides_per_variant=5,
        num_epitopes_per_vaccine_peptide=500,
    )
    config = vaccine_config_from_args(args)
    eq_(config.vaccine_peptide_length, 30)
    eq_(config.padding_around_mutation, 10)
    eq_(config.max_vaccine_peptides_per_variant, 5)
    eq_(config.num_mutant_epitopes_to_keep, 500)


def test_vaccine_config_from_args_legacy_mutation_namespace_attr():
    """Backward compatibility for callers using old namespace attribute."""
    args = argparse.Namespace(
        config=None,
        vaccine_peptide_length=None,
        padding_around_mutation=None,
        max_vaccine_peptides_per_mutation=7,
        num_epitopes_per_vaccine_peptide=None,
    )
    config = vaccine_config_from_args(args)
    eq_(config.max_vaccine_peptides_per_variant, 7)


def test_vaccine_config_from_args_yaml_config_file():
    """Test loading vaccine config from YAML file"""
    yaml_content = """
vaccine_peptides:
  generation:
    lengths: 40
    padding_around_mutation: 8
  keep:
    per_mutation: 3
    max_epitopes_per_candidate: 2000
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        config = vaccine_config_from_args(args)
        eq_(config.vaccine_peptide_length, 40)
        eq_(config.padding_around_mutation, 8)
        eq_(config.max_vaccine_peptides_per_variant, 3)
        eq_(config.num_mutant_epitopes_to_keep, 2000)
    finally:
        os.unlink(config_path)


def test_vaccine_config_from_args_cli_overrides_yaml():
    """Test that CLI arguments override YAML config values"""
    yaml_content = """
vaccine_peptides:
  generation:
    lengths: 40
  keep:
    per_mutation: 3
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=50,  # Override YAML's 40
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        config = vaccine_config_from_args(args)
        eq_(config.vaccine_peptide_length, 50)  # CLI override
        eq_(config.max_vaccine_peptides_per_variant, 3)  # from YAML
    finally:
        os.unlink(config_path)


def test_vaccine_config_from_args_yaml_without_vaccine_peptides_section():
    """Test YAML file without vaccine_peptides section uses defaults"""
    yaml_content = """
epitopes:
  filters:
    min_score: 0.01
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        config = vaccine_config_from_args(args)
        eq_(config.vaccine_peptide_length, 25)  # default
    finally:
        os.unlink(config_path)


def test_vaccine_config_from_args_nested_yaml_config_file():
    yaml_content = """
vaccine_peptides:
  generation:
    lengths: [31]
    padding_around_mutation: 9
  keep:
    per_mutation: 2
    max_epitopes_per_candidate: 42
    score_fraction_of_best: 0.95
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
            config_set_overrides=None,
            config_expr_overrides=None,
        )
        config = vaccine_config_from_args(args)
        eq_(config.vaccine_peptide_length, 31)
        eq_(config.padding_around_mutation, 9)
        eq_(config.max_vaccine_peptides_per_variant, 2)
        eq_(config.num_mutant_epitopes_to_keep, 42)
        eq_(config.score_fraction_of_best, 0.95)
    finally:
        os.unlink(config_path)


def test_vaccine_config_from_args_nested_manufacturability():
    yaml_content = """
vaccine_peptides:
  manufacturability:
    max_c_terminal_hydropathy: 2.0
    max_kmer_hydropathy_high_priority: 3.0
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
            config_set_overrides=None,
            config_expr_overrides=None,
        )
        config = vaccine_config_from_args(args)
        eq_(config.max_c_terminal_hydropathy, 2.0)
        eq_(config.max_kmer_hydropathy_high_priority, 3.0)
        # Unset fields keep defaults
        eq_(config.min_kmer_hydropathy, 0.0)
        eq_(config.max_kmer_hydropathy_low_priority, 1.5)
    finally:
        os.unlink(config_path)


def test_vaccine_config_from_args_multiple_lengths_not_supported():
    yaml_content = """
vaccine_peptides:
  generation:
    lengths: [25, 31]
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
            config_set_overrides=None,
            config_expr_overrides=None,
        )
        with pytest.raises(ValueError):
            vaccine_config_from_args(args)
    finally:
        os.unlink(config_path)


# =============================================================================
# Combined Config Tests
# =============================================================================

def test_combined_config_both_configs_from_single_yaml():
    """Test loading both epitope and vaccine configs from one YAML file"""
    yaml_content = """
epitopes:
  filters:
    min_score: 0.02
  scoring:
    derived_fields:
      affinity_score:
        midpoint: 400.0

vaccine_peptides:
  generation:
    lengths: 30
  keep:
    per_mutation: 2
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
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        vaccine_config = vaccine_config_from_args(vaccine_args)
        eq_(vaccine_config.vaccine_peptide_length, 30)
        eq_(vaccine_config.max_vaccine_peptides_per_variant, 2)
    finally:
        os.unlink(config_path)


def test_combined_config_mixed_cli_and_yaml_overrides():
    """Test mixing CLI and YAML overrides for both configs"""
    yaml_content = """
epitopes:
  filters:
    min_score: 0.01

vaccine_peptides:
  generation:
    lengths: 30
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
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        vaccine_config = vaccine_config_from_args(vaccine_args)
        eq_(vaccine_config.vaccine_peptide_length, 40)  # CLI
    finally:
        os.unlink(config_path)


# =============================================================================
# CLI Argument Parser Tests
# =============================================================================

def test_cli_add_epitope_prediction_args():
    """Test that add_epitope_prediction_args adds the expected arguments"""
    parser = argparse.ArgumentParser()
    add_epitope_prediction_args(parser)

    # Parse with min-epitope-score
    args = parser.parse_args(['--min-epitope-score', '0.05'])
    eq_(args.min_epitope_score, 0.05)


def test_cli_add_vaccine_peptide_args():
    """Test that add_vaccine_peptide_args adds the expected arguments"""
    parser = argparse.ArgumentParser()
    add_vaccine_peptide_args(parser)

    # Parse with various vaccine args
    args = parser.parse_args([
        '--vaccine-peptide-length', '30',
        '--padding-around-mutation', '10',
        '--max-vaccine-peptides-per-variant', '5',
    ])
    eq_(args.vaccine_peptide_length, 30)
    eq_(args.padding_around_mutation, 10)
    eq_(args.max_vaccine_peptides_per_variant, 5)


def test_cli_add_config_override_args():
    parser = argparse.ArgumentParser()
    add_config_override_args(parser)

    args = parser.parse_args(
        [
            "--config-value", "epitope_config.min_epitope_score=0.05",
            "--config-text", "vaccine_peptides.ranking.score=sqrt(n_alt_reads)",
        ]
    )
    eq_(args.config_set_overrides, ["epitope_config.min_epitope_score=0.05"])
    eq_(args.config_expr_overrides, ["vaccine_peptides.ranking.score=sqrt(n_alt_reads)"])


def test_load_vaxrank_config_preserves_mixed_override_order():
    parser = argparse.ArgumentParser()
    add_config_override_args(parser)

    args_expr_then_set = parser.parse_args([
        "--config-text", "epitopes.scoring.mode=percentile_rank",
        "--config-value", "epitopes.scoring.mode=affinity",
    ])
    args_set_then_expr = parser.parse_args([
        "--config-value", "epitopes.scoring.mode=affinity",
        "--config-text", "epitopes.scoring.mode=percentile_rank",
    ])

    config_expr_then_set = load_vaxrank_config(args_expr_then_set)
    config_set_then_expr = load_vaxrank_config(args_set_then_expr)

    eq_(config_expr_then_set["epitopes"]["scoring"]["mode"], "affinity")
    eq_(config_set_then_expr["epitopes"]["scoring"]["mode"], "percentile_rank")


def test_cli_add_vaccine_peptide_args_legacy_mutation_alias():
    """Legacy mutation flag maps to variant destination for compatibility."""
    parser = argparse.ArgumentParser()
    add_vaccine_peptide_args(parser)

    args = parser.parse_args(['--max-vaccine-peptides-per-mutation', '5'])
    eq_(args.max_vaccine_peptides_per_variant, 5)


def test_cli_vaccine_peptide_args_defaults():
    """Test that vaccine peptide args default to None for config override"""
    parser = argparse.ArgumentParser()
    add_vaccine_peptide_args(parser)

    args = parser.parse_args([])
    eq_(args.vaccine_peptide_length, None)
    eq_(args.padding_around_mutation, None)
    eq_(args.max_vaccine_peptides_per_variant, None)


# =============================================================================
# CLI Version Tests
# =============================================================================

def test_cli_version_cached_parser_has_version_flag():
    """Test --version works in cached mode"""
    from vaxrank.cli.arg_parser import choose_arg_parser

    parser = choose_arg_parser(["--input-json-file", "dummy.json"])
    with pytest.raises(SystemExit) as excinfo:
        parser.parse_args(["--version"])
    eq_(excinfo.value.code, 0)


def test_choose_arg_parser_detects_input_json_file_equals_form():
    """--input-json-file=<path> should select the cached-mode parser."""
    from vaxrank.cli.arg_parser import choose_arg_parser

    parser = choose_arg_parser(["--input-json-file=dummy.json"])
    args = parser.parse_args(["--input-json-file=dummy.json", "--output-ascii-report", "out.txt"])
    eq_(args.input_json_file, "dummy.json")


# =============================================================================
# Error Handling Tests
# =============================================================================

def test_config_error_nonexistent_config_file():
    """Test that nonexistent config file raises appropriate error"""
    args = argparse.Namespace(
        config='/nonexistent/path/config.yaml',
        min_epitope_score=None,
    )
    with pytest.raises(FileNotFoundError):
        epitope_config_from_args(args)


def test_config_error_invalid_yaml_syntax():
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


def test_config_error_invalid_config_value_type():
    """Test that invalid value type raises appropriate error"""
    yaml_content = """
epitopes:
  filters:
    min_score: "not_a_number"
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


def test_load_vaxrank_config_rejects_schema_version():
    yaml_content = """
schema_version: 2
epitope_config:
  min_epitope_score: 0.01
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        with pytest.raises(ValueError, match="schema_version is not supported"):
            load_vaxrank_config(config_path=config_path)
    finally:
        os.unlink(config_path)


def test_load_vaxrank_config_rejects_unsupported_top_level_section():
    yaml_content = """
self_proteome:
  exclude_gene_ids: ["ENSG000001"]
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        with pytest.raises(msgspec.ValidationError, match="self_proteome"):
            load_vaxrank_config(config_path=config_path)
    finally:
        os.unlink(config_path)


def test_load_vaxrank_config_rejects_unsupported_nested_section():
    yaml_content = """
vaccine_peptides:
  ranking:
    score: sqrt(n_alt_reads)
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        with pytest.raises(msgspec.ValidationError, match="ranking"):
            load_vaxrank_config(config_path=config_path)
    finally:
        os.unlink(config_path)


def test_load_vaxrank_config_set_override():
    args = argparse.Namespace(
        config=None,
        config_set_overrides=[
            "epitopes.scoring.derived_fields.affinity_score.midpoint=425.0",
            "vaccine_peptides.generation.lengths=[31]",
        ],
        config_expr_overrides=None,
    )
    config = load_vaxrank_config(args)
    eq_(config["epitopes"]["scoring"]["derived_fields"]["affinity_score"]["midpoint"], 425.0)
    eq_(config["vaccine_peptides"]["generation"]["lengths"], [31])


def test_load_vaxrank_config_expr_override():
    args = argparse.Namespace(
        config=None,
        config_set_overrides=None,
        config_expr_overrides=[
            "epitopes.scoring.mode=percentile_rank"
        ],
    )
    config = load_vaxrank_config(args)
    eq_(config["epitopes"]["scoring"]["mode"], "percentile_rank")


# =============================================================================
# Legacy Config Key Rejection Tests
# =============================================================================

def test_load_vaxrank_config_rejects_legacy_epitope_config_key():
    """Legacy epitope_config top-level key is now rejected."""
    yaml_content = """
epitope_config:
  min_epitope_score: 0.01
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        with pytest.raises(msgspec.ValidationError, match="epitope_config"):
            load_vaxrank_config(config_path=config_path)
    finally:
        os.unlink(config_path)


def test_load_vaxrank_config_rejects_legacy_vaccine_config_key():
    """Legacy vaccine_config top-level key is now rejected."""
    yaml_content = """
vaccine_config:
  vaccine_peptide_length: 30
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        with pytest.raises(msgspec.ValidationError, match="vaccine_config"):
            load_vaxrank_config(config_path=config_path)
    finally:
        os.unlink(config_path)


# =============================================================================
# Range Validation Tests
# =============================================================================

def test_epitope_config_rejects_negative_midpoint():
    with pytest.raises(ValueError, match="logistic_epitope_score_midpoint"):
        EpitopeConfig(logistic_epitope_score_midpoint=-1)


def test_epitope_config_rejects_zero_midpoint():
    with pytest.raises(ValueError, match="logistic_epitope_score_midpoint"):
        EpitopeConfig(logistic_epitope_score_midpoint=0)


def test_epitope_config_rejects_zero_width():
    with pytest.raises(ValueError, match="logistic_epitope_score_width"):
        EpitopeConfig(logistic_epitope_score_width=0)


def test_epitope_config_rejects_negative_width():
    with pytest.raises(ValueError, match="logistic_epitope_score_width"):
        EpitopeConfig(logistic_epitope_score_width=-10)


def test_epitope_config_rejects_negative_min_score():
    with pytest.raises(ValueError, match="min_epitope_score"):
        EpitopeConfig(min_epitope_score=-0.1)


def test_epitope_config_rejects_negative_binding_affinity_cutoff():
    with pytest.raises(ValueError, match="binding_affinity_cutoff"):
        EpitopeConfig(binding_affinity_cutoff=-100)


def test_epitope_config_rejects_zero_binding_affinity_cutoff():
    with pytest.raises(ValueError, match="binding_affinity_cutoff"):
        EpitopeConfig(binding_affinity_cutoff=0)


def test_epitope_config_rejects_invalid_scoring_mode():
    with pytest.raises(ValueError, match="scoring_mode"):
        EpitopeConfig(scoring_mode="invalid")


def test_epitope_config_rejects_percentile_rank_cutoff_zero():
    with pytest.raises(ValueError, match="percentile_rank_cutoff"):
        EpitopeConfig(percentile_rank_cutoff=0)


def test_epitope_config_rejects_percentile_rank_cutoff_over_100():
    with pytest.raises(ValueError, match="percentile_rank_cutoff"):
        EpitopeConfig(percentile_rank_cutoff=101)


def test_vaccine_config_rejects_zero_peptide_length():
    with pytest.raises(ValueError, match="vaccine_peptide_length"):
        VaccineConfig(vaccine_peptide_length=0)


def test_vaccine_config_rejects_negative_peptide_length():
    with pytest.raises(ValueError, match="vaccine_peptide_length"):
        VaccineConfig(vaccine_peptide_length=-5)


def test_vaccine_config_rejects_negative_padding():
    with pytest.raises(ValueError, match="padding_around_mutation"):
        VaccineConfig(padding_around_mutation=-1)


def test_vaccine_config_rejects_negative_max_peptides():
    with pytest.raises(ValueError, match="max_vaccine_peptides_per_variant"):
        VaccineConfig(max_vaccine_peptides_per_variant=-1)


def test_vaccine_config_rejects_negative_epitopes_to_keep():
    with pytest.raises(ValueError, match="num_mutant_epitopes_to_keep"):
        VaccineConfig(num_mutant_epitopes_to_keep=-1)


def test_vaccine_config_rejects_score_fraction_zero():
    with pytest.raises(ValueError, match="score_fraction_of_best"):
        VaccineConfig(score_fraction_of_best=0)


def test_vaccine_config_rejects_score_fraction_over_1():
    with pytest.raises(ValueError, match="score_fraction_of_best"):
        VaccineConfig(score_fraction_of_best=1.5)


def test_vaccine_config_allows_zero_max_peptides_per_variant():
    """0 means disabled — should be allowed."""
    config = VaccineConfig(max_vaccine_peptides_per_variant=0)
    eq_(config.max_vaccine_peptides_per_variant, 0)


def test_vaccine_config_allows_zero_epitopes_to_keep():
    """0 means keep all — should be allowed."""
    config = VaccineConfig(num_mutant_epitopes_to_keep=0)
    eq_(config.num_mutant_epitopes_to_keep, 0)


def test_epitope_config_allows_zero_min_score():
    """0 means no filtering — should be allowed."""
    config = EpitopeConfig(min_epitope_score=0)
    eq_(config.min_epitope_score, 0)


def test_validation_through_yaml_loading():
    """Invalid values via YAML should raise during config construction."""
    yaml_content = """
vaccine_peptides:
  generation:
    lengths: -5
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
        )
        # msgspec.convert wraps ValueError in ValidationError
        with pytest.raises((ValueError, msgspec.ValidationError), match="vaccine_peptide_length"):
            vaccine_config_from_args(args)
    finally:
        os.unlink(config_path)


# =============================================================================
# Conflict Detection Tests
# =============================================================================

def test_extract_vaccine_config_kwargs_errors_on_conflicting_epitope_keep():
    """Setting num_mutant_epitopes_to_keep via both paths should error."""
    from vaxrank.config.loader import extract_vaccine_config_kwargs

    config = {
        "epitopes": {"keep": {"top_n_per_candidate": 10}},
        "vaccine_peptides": {"keep": {"max_epitopes_per_candidate": 20}},
    }
    with pytest.raises(ValueError, match="Cannot set both"):
        extract_vaccine_config_kwargs(config)


def test_extract_vaccine_config_kwargs_accepts_single_epitope_keep_path():
    """Only one path set should work fine."""
    from vaxrank.config.loader import extract_vaccine_config_kwargs

    config_a = {
        "epitopes": {"keep": {"top_n_per_candidate": 10}},
    }
    kwargs_a = extract_vaccine_config_kwargs(config_a)
    eq_(kwargs_a["num_mutant_epitopes_to_keep"], 10)

    config_b = {
        "vaccine_peptides": {"keep": {"max_epitopes_per_candidate": 20}},
    }
    kwargs_b = extract_vaccine_config_kwargs(config_b)
    eq_(kwargs_b["num_mutant_epitopes_to_keep"], 20)
