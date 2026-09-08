"""Exercise the actual Isovar factory through Vaxrank's public RNA entry point."""

from importlib.resources import files
from types import SimpleNamespace
from unittest.mock import Mock

import msgspec
import pytest
from isovar.cli import protein_sequence_creator_from_args

from vaxrank.cli import make_vaxrank_arg_parser, run_vaxrank_from_parsed_args
from vaxrank.cli.isovar_config_args import resolve_isovar_args
from vaxrank.vaccine_config import VaccineConfig


@pytest.fixture
def run_config(monkeypatch, tmp_path):
    from vaxrank.cli import entry_point

    # Stub only I/O and expensive downstream ranking, not config or the creator.
    monkeypatch.setattr(entry_point, "variant_collection_from_args", lambda args: [])
    monkeypatch.setattr(entry_point, "filter_unannotatable_variants", lambda v: v)
    monkeypatch.setattr(entry_point, "alignment_file_from_args", lambda args: None)
    monkeypatch.setattr(entry_point, "read_collector_from_args", lambda args: None)
    monkeypatch.setattr(entry_point, "predictors_from_args", lambda args: [Mock()])
    monkeypatch.setattr(entry_point, "run_vaxrank", lambda **kwargs: kwargs)
    captured = {}

    def capture(**kwargs):
        captured["creator"] = kwargs["protein_sequence_creator"]
        return []

    monkeypatch.setattr(entry_point, "run_isovar", capture)

    def run(flags=(), yaml=None):
        options = ["--vcf", "unused.vcf", "--bam", "unused.bam",
                   "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01"]
        if yaml is not None:
            config = tmp_path / "config.yaml"
            config.write_text(yaml)
            options += ["--config", str(config)]
        args = make_vaxrank_arg_parser().parse_args(options + list(flags))
        result = run_vaxrank_from_parsed_args(args)
        return args, captured["creator"], result

    return run


@pytest.mark.parametrize("peptide,context", [(15, 29), (25, 49), (30, 59)])
@pytest.mark.parametrize("source", ["cli", "yaml"])
def test_adaptive_context_uses_resolved_peptide_size(run_config, peptide, context, source):
    flags = ["--vaccine-peptide-length", str(peptide)] if source == "cli" else []
    yaml = f"vaccine_peptides:\n  preferred_length: {peptide}\n" if source == "yaml" else None
    args, creator, result = run_config(flags, yaml)
    assert creator.protein_sequence_length == context
    assert args.protein_sequence_length == context  # effective run provenance
    assert args.protein_context_peptide_length == peptide
    assert creator.protein_sequence_preference == args.protein_sequence_preference == "balanced"
    assert creator.min_protein_sequence_support_fraction == args.min_protein_sequence_support_fraction == .85
    assert creator.min_variant_sequence_coverage == args.min_variant_sequence_coverage == 2
    assert result["vaccine_config"].padding_around_mutation == 5  # DNA fallback unchanged


def test_printed_default_yaml_does_not_reenable_old_padding(run_config):
    yaml = files("vaxrank.config").joinpath("default.yaml").read_text()
    args, creator, _ = run_config(yaml=yaml)
    assert args.protein_sequence_length == creator.protein_sequence_length == 49


@pytest.mark.parametrize("padding", [0, 5, 12])
@pytest.mark.parametrize("source", ["cli", "yaml"])
def test_explicit_legacy_padding_is_preserved(run_config, padding, source):
    flags = ["--padding-around-mutation", str(padding)] if source == "cli" else []
    yaml = f"vaccine_peptides:\n  padding_around_mutation: {padding}\n" if source == "yaml" else None
    args, creator, _ = run_config(flags, yaml)
    assert creator.protein_sequence_length == args.protein_sequence_length == 25 + 2 * padding


def test_explicit_context_length_wins_over_yaml_and_padding(run_config):
    args, creator, _ = run_config(
        ["--protein-sequence-length", "40", "--padding-around-mutation", "12"],
        "isovar:\n  protein_sequence_length: 60\nvaccine_peptides:\n  padding_around_mutation: 7\n")
    assert creator.protein_sequence_length == args.protein_sequence_length == 40


def test_yaml_context_length_wins_over_legacy_padding(run_config):
    _, creator, _ = run_config(
        ["--padding-around-mutation", "12"], "isovar:\n  protein_sequence_length: 40\n")
    assert creator.protein_sequence_length == 40


def test_isovar_peptide_override_is_separate_from_vaccine_size(run_config):
    args, creator, result = run_config(
        ["--protein-context-peptide-length", "30"],
        "vaccine_peptides:\n  preferred_length: 15\nisovar:\n  protein_context_peptide_length: 20\n")
    assert creator.protein_sequence_length == 59
    assert args.protein_context_peptide_length == 30
    assert result["vaccine_config"].preferred_peptide_length == 15


def test_yaml_controls_and_cli_explicit_defaults(run_config):
    yaml = ("isovar:\n  protein_sequence_preference: context\n"
            "  min_protein_sequence_support_fraction: 0.95\n"
            "  min_variant_sequence_coverage: 5\n")
    _, creator, _ = run_config(yaml=yaml)
    assert (creator.protein_sequence_preference, creator.min_protein_sequence_support_fraction,
            creator.min_variant_sequence_coverage) == ("context", .95, 5)
    _, creator, _ = run_config(
        ["--protein-sequence-preference", "balanced", "--min-protein-sequence-support-fraction", ".85",
         "--min-variant-sequence-coverage", "2"], yaml)
    assert (creator.protein_sequence_preference, creator.min_protein_sequence_support_fraction,
            creator.min_variant_sequence_coverage) == ("balanced", .85, 2)


def test_set_overrides_and_explicit_zero_coverage(run_config):
    _, creator, _ = run_config([
        "--config-value", "isovar.protein_context_peptide_length=30",
        "--config-value", "isovar.min_variant_sequence_coverage=0",
        "--config-value", "isovar.min_protein_sequence_support_fraction=0.95",
    ])
    assert creator.protein_sequence_length == 59
    assert creator.min_variant_sequence_coverage == 0
    assert creator.min_protein_sequence_support_fraction == .95


def test_historical_policy_remains_explicit(run_config):
    args, creator, _ = run_config([
        "--protein-sequence-length", "20", "--protein-sequence-preference", "support"])
    assert creator.protein_sequence_preference == args.protein_sequence_preference == "support"
    assert creator.candidate_context_lengths() == [20]


@pytest.mark.parametrize("name,value", [
    ("protein_sequence_preference", "invalid"),
    ("min_protein_sequence_support_fraction", "1.01"),
    ("min_protein_sequence_support_fraction", ".nan"),
    ("min_variant_sequence_coverage", "-1"),
    ("min_variant_sequence_coverage", "true"),
    ("protein_context_peptide_length", "0"),
    ("protein_sequence_length", "-1"),
    ("min_variant_sequence_coverge", "2"),
])
def test_invalid_yaml_fails_before_io(run_config, name, value):
    with pytest.raises((ValueError, msgspec.ValidationError)):
        run_config(yaml=f"isovar:\n  {name}: {value}\n")


def test_null_yaml_and_null_python_namespace_use_defaults(run_config):
    yaml = "isovar:\n" + "".join(f"  {name}: null\n" for name in (
        "protein_sequence_length", "protein_context_peptide_length", "protein_sequence_preference",
        "min_protein_sequence_support_fraction", "min_variant_sequence_coverage"))
    _, creator, _ = run_config(yaml=yaml)
    assert creator.protein_sequence_length == 49
    args = SimpleNamespace(protein_sequence_length=None, protein_context_peptide_length=None,
                           protein_sequence_preference=None, min_protein_sequence_support_fraction=None,
                           min_variant_sequence_coverage=None)
    resolve_isovar_args(args, VaccineConfig(preferred_peptide_length=30, max_peptide_length=30), {})
    assert args.protein_sequence_length is None  # Isovar, not Vaxrank, derives this
    assert args.protein_context_peptide_length == 30
    assert args.min_variant_sequence_coverage == 2


def test_python_callers_explicit_requests_survive_resolution():
    args = make_vaxrank_arg_parser().parse_args([
        "--vcf", "unused", "--bam", "unused", "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01"])
    args.protein_sequence_length = 40
    args.min_variant_sequence_coverage = 5
    resolve_isovar_args(args, VaccineConfig(preferred_peptide_length=30, max_peptide_length=30), {})
    creator = protein_sequence_creator_from_args(args)
    assert creator.protein_sequence_length == 40
    assert creator.min_variant_sequence_coverage == 5
