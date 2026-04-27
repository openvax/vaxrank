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

"""Tests for configurable manufacturability rules and combined_score on
``VaccinePeptide``. Separate from ``test_epitope_prediction.py`` because
these don't require an MHC predictor."""

import argparse
import os
import tempfile
from types import SimpleNamespace

import numpy as np
import pytest

from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
from vaxrank.manufacturability import DEFAULT_MANUFACTURABILITY_RULES
from vaxrank.vaccine_peptide import VaccinePeptide


def _make_fragment(seq="A" * 25, n_alt_reads=9):
    """Minimal duck-typed MutantProteinFragment for VaccinePeptide tests."""
    return SimpleNamespace(
        amino_acids=seq,
        n_alt_reads=n_alt_reads,
        n_alt_reads_supporting_protein_sequence=n_alt_reads,
        n_mutant_amino_acids=1,
        mutation_distance_from_edge=5,
    )


def _make_mutant_epitope(ic50=10.0):
    """Duck-typed EpitopePrediction that overlaps mutation, not in reference."""
    return SimpleNamespace(
        ic50=ic50,
        percentile_rank=0.5,
        overlaps_mutation=True,
        occurs_in_reference=False,
        logistic_epitope_score=lambda **kw: 0.9,
    )


def test_default_manufacturability_tuple_length():
    """Default sort tuple length = 10 (legacy behavior, unchanged)."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
    )
    tup = vp.peptide_synthesis_difficulty_score_tuple()
    assert len(tup) == 10
    assert len(DEFAULT_MANUFACTURABILITY_RULES) == 10


def test_custom_rules_shorten_tuple():
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
        manufacturability_rules=["cysteine_count", "cterm_hydropathy"],
    )
    tup = vp.peptide_synthesis_difficulty_score_tuple(
        rules=vp.manufacturability_rules,
    )
    assert len(tup) == 2


def test_combined_score_default_mode_is_legacy():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
    )
    # Legacy: sqrt(9) * mutant_epitope_score
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(np.sqrt(9) * vp.mutant_epitope_score)


def test_combined_score_reads_times_epitope_mode():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
        combined_score_mode="reads_times_epitope",
    )
    # expression_score is an honest metric; mode only affects combined_score
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(9.0 * vp.mutant_epitope_score)


def test_combined_score_epitope_only_mode():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
        combined_score_mode="epitope_only",
    )
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(vp.mutant_epitope_score)


def test_vaccine_peptide_rejects_unknown_combined_score_mode():
    with pytest.raises(ValueError, match="combined_score_mode"):
        VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitope_predictions=[_make_mutant_epitope()],
            combined_score_mode="garbage",
        )


def test_end_to_end_yaml_rules_through_vaccine_config_to_peptide():
    """YAML with manufacturability.rules loads through load_vaxrank_config,
    then vaccine_config_from_args, producing a tuple on VaccineConfig.
    A VaccinePeptide built with those rules gets a sort tuple of the
    expected length."""
    yaml_content = """
manufacturability:
  rules:
    - cysteine_count
    - cterm_hydropathy
    - aspartate_proline
vaccine_peptides:
  combined_score_mode: epitope_only
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False
    ) as f:
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
        vc = vaccine_config_from_args(args)
        # Type contract: VaccineConfig declares Optional[tuple[str, ...]]
        assert vc.manufacturability_rules == (
            "cysteine_count", "cterm_hydropathy", "aspartate_proline",
        )
        assert isinstance(vc.manufacturability_rules, tuple)
        assert vc.combined_score_mode == "epitope_only"

        vp = VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitope_predictions=[_make_mutant_epitope()],
            manufacturability_rules=vc.manufacturability_rules,
            combined_score_mode=vc.combined_score_mode,
        )
        tup = vp.peptide_synthesis_difficulty_score_tuple(
            rules=vp.manufacturability_rules,
        )
        assert len(tup) == 3
        # combined_score under epitope_only is just the epitope score
        assert vp.combined_score == pytest.approx(vp.mutant_epitope_score)
    finally:
        os.unlink(config_path)


def test_vaccine_peptide_roundtrip_preserves_custom_config():
    """to_dict / from_dict preserves non-default rules + score mode."""
    fragment = _make_fragment()
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
        manufacturability_rules=["cysteine_count"],
        combined_score_mode="epitope_only",
    )
    d = vp.to_dict()
    assert d["manufacturability_rules"] == ["cysteine_count"]
    assert d["combined_score_mode"] == "epitope_only"


def test_vaccine_peptide_roundtrip_omits_defaults():
    """to_dict doesn't include None/default rule + mode fields."""
    fragment = _make_fragment()
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
    )
    d = vp.to_dict()
    assert "manufacturability_rules" not in d
    assert "combined_score_mode" not in d
    assert "ranking_rules" not in d


def test_vaccine_peptide_ranking_rules_tuple_coercion():
    """List of ranking rule names coerced to tuple (same contract as
    manufacturability_rules)."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
        ranking_rules=["mutant_epitope_score", "n_alt_reads"],
    )
    assert vp.ranking_rules == ("mutant_epitope_score", "n_alt_reads")
    assert isinstance(vp.ranking_rules, tuple)
    # Sort tuple length reflects the custom 2-rule list.
    assert len(vp.lexicographic_sort_key()) == 2


def test_vaccine_peptide_custom_ranking_rules_in_to_dict():
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
        ranking_rules=["mutant_epitope_score"],
    )
    d = vp.to_dict()
    assert d["ranking_rules"] == ["mutant_epitope_score"]


def test_end_to_end_yaml_ranking_rules():
    """YAML with vaccine_peptides.ranking_rules loads through the config
    pipeline to a VaccineConfig tuple and threads into VaccinePeptide."""
    yaml_content = """
vaccine_peptides:
  ranking_rules:
    - mutant_epitope_score
    - n_alt_reads
    - manufacturability
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
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
        vc = vaccine_config_from_args(args)
        assert vc.ranking_rules == (
            "mutant_epitope_score", "n_alt_reads", "manufacturability",
        )
        assert isinstance(vc.ranking_rules, tuple)
        # Custom rules drop the extra_score_tuple tiers (n_alt_reads_supporting,
        # wildtype_epitope_score, n_mutant_amino_acids, mutation_distance_from_edge).
        vp = VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitope_predictions=[_make_mutant_epitope()],
            ranking_rules=vc.ranking_rules,
        )
        # 2 essential + 10 manufacturability defaults = 12 entries
        assert len(vp.lexicographic_sort_key()) == 12
    finally:
        os.unlink(config_path)


def test_vaccine_config_rejects_unknown_ranking_rule():
    from vaxrank.vaccine_config import VaccineConfig

    with pytest.raises(ValueError, match="Unknown ranking rule"):
        VaccineConfig(ranking_rules=("not_a_real_rule",))


# ── require_mutant_epitopes_in_variant ───────────────────────────────────────

def test_vaccine_config_default_requires_mutant_epitopes():
    """Legacy behavior preserved: default is to drop variants with no
    mutant-overlapping epitopes."""
    from vaxrank.vaccine_config import VaccineConfig
    assert VaccineConfig().require_mutant_epitopes_in_variant is True


def test_require_mutant_epitopes_yaml_override(tmp_path):
    """vaccine_peptides.require_mutant_epitopes_in_variant: false in YAML
    propagates into the config struct."""
    import yaml
    from vaxrank.config.loader import (
        extract_vaccine_config_kwargs, load_vaxrank_config)
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump({
        "vaccine_peptides": {"require_mutant_epitopes_in_variant": False}
    }))
    import argparse
    args = argparse.Namespace(config=[str(cfg_path)],
                              config_set_overrides=None,
                              config_expr_overrides=None)
    merged = load_vaxrank_config(args)
    kwargs = extract_vaccine_config_kwargs(merged)
    assert kwargs["require_mutant_epitopes_in_variant"] is False


def test_create_vaccine_peptides_dict_respects_require_mutant_epitopes():
    """When require_mutant_epitopes_in_variant=False, variants with no
    mutant-overlapping vaccine peptides should still appear in the
    output dict (instead of being silently dropped)."""
    from unittest.mock import MagicMock, patch
    from vaxrank.core_logic import create_vaccine_peptides_dict
    from vaxrank.vaccine_config import VaccineConfig

    # Fake a single isovar result + a non-mutant-covering vaccine peptide.
    iso = MagicMock()
    iso.variant = "v1"
    no_epi_peptide = MagicMock()
    no_epi_peptide.contains_mutant_epitopes.return_value = False

    with patch("vaxrank.core_logic.vaccine_peptides_for_variant",
               return_value=[no_epi_peptide]):
        # Default (require=True): variant dropped.
        out_strict = create_vaccine_peptides_dict(
            isovar_results=[iso], mhc_predictor=MagicMock(),
            vaccine_config=VaccineConfig())
        assert out_strict == {}

        # require=False: variant kept.
        out_lax = create_vaccine_peptides_dict(
            isovar_results=[iso], mhc_predictor=MagicMock(),
            vaccine_config=VaccineConfig(
                require_mutant_epitopes_in_variant=False))
        assert "v1" in out_lax


# ── Bundled default.yaml ─────────────────────────────────────────────────────

def test_bundled_default_yaml_is_present_and_parses():
    """The shipped vaxrank/config/default.yaml must be installed and
    parse as valid YAML."""
    from importlib.resources import files
    import yaml
    text = files("vaxrank.config").joinpath("default.yaml").read_text()
    parsed = yaml.safe_load(text)
    # Top-level sections present
    assert set(parsed) >= {"epitopes", "vaccine_peptides", "manufacturability"}


def test_bundled_default_yaml_round_trips_to_default_configs(tmp_path):
    """Loading the shipped default.yaml produces configs whose effective
    behavior matches plain EpitopeConfig() / VaccineConfig() defaults.

    Equality isn't quite ``==`` because the YAML lists ranking_rules /
    manufacturability_rules explicitly (so users can see the names) while
    the dataclasses default these to ``None`` (which the rule executor
    interprets as DEFAULT_RANKING_RULES / DEFAULT_MANUFACTURABILITY_RULES).
    Compare resolved tuples instead of raw struct equality.
    """
    from importlib.resources import files
    import argparse
    from vaxrank.cli.epitope_config_args import epitope_config_from_args
    from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.manufacturability import DEFAULT_MANUFACTURABILITY_RULES
    from vaxrank.ranking import DEFAULT_RANKING_RULES
    from vaxrank.vaccine_config import VaccineConfig

    src = files("vaxrank.config").joinpath("default.yaml").read_text()
    path = tmp_path / "default.yaml"
    path.write_text(src)

    args = argparse.Namespace(
        config=[str(path)], config_set_overrides=None,
        config_expr_overrides=None,
        min_epitope_score=None, scoring_mode=None,
        default_affinity_predictor=None,
        default_stability_predictor=None,
        default_presentation_predictor=None,
        vaccine_peptide_length=None,
        padding_around_mutation=None,
        max_vaccine_peptides_per_mutation=None,
        num_epitopes_per_vaccine_peptide=None,
    )
    ec = epitope_config_from_args(args)
    vc = vaccine_config_from_args(args)
    # Epitope side has no list-typed fields, so == works.
    assert ec == EpitopeConfig()
    # Vaccine side: effective rule tuples should match the defaults.
    assert vc.ranking_rules == DEFAULT_RANKING_RULES
    assert vc.manufacturability_rules == DEFAULT_MANUFACTURABILITY_RULES
    # All scalar fields match defaults too.
    defaults = VaccineConfig()
    for field in (
        "preferred_peptide_length", "min_peptide_length",
        "max_peptide_length", "padding_around_mutation",
        "max_vaccine_peptides_per_variant", "num_mutant_epitopes_to_keep",
        "score_fraction_of_best", "combined_score_mode",
        "max_c_terminal_hydropathy", "min_kmer_hydropathy",
        "max_kmer_hydropathy_low_priority", "max_kmer_hydropathy_high_priority",
        "require_mutant_epitopes_in_variant",
    ):
        assert getattr(vc, field) == getattr(defaults, field), (
            f"default.yaml drift: {field}")


def test_default_yaml_uncommented_keys_match_schema():
    """Every uncommented key under each top-level section in default.yaml
    must be a recognized schema field. Catches drift like keys placed in
    the wrong section or renamed without updating the YAML."""
    from importlib.resources import files
    import yaml
    from vaxrank.config.schema import (
        EpitopesConfigSchema,
        ManufacturabilityConfigSchema,
        VaccinePeptidesConfigSchema,
    )

    text = files("vaxrank.config").joinpath("default.yaml").read_text()
    parsed = yaml.safe_load(text)

    schema_for_section = {
        "epitopes": EpitopesConfigSchema,
        "vaccine_peptides": VaccinePeptidesConfigSchema,
        "manufacturability": ManufacturabilityConfigSchema,
    }
    for section, schema_cls in schema_for_section.items():
        assert section in parsed, f"Missing section in default.yaml: {section}"
        section_dict = parsed[section] or {}
        # Drop the nested manufacturability sub-object inside vaccine_peptides
        # (handled separately above) so it doesn't trigger forbid_unknown_fields
        # when it's actually a different schema.
        if section == "vaccine_peptides":
            section_dict = {k: v for k, v in section_dict.items()
                            if k != "manufacturability"}
        # forbid_unknown_fields=True on each schema means msgspec will raise
        # if any key in section_dict is not a declared schema field.
        import msgspec
        msgspec.convert(section_dict, schema_cls)


def test_default_yaml_commented_examples_are_valid():
    """Commented-out example lines in default.yaml are advertised as 'just
    uncomment to enable' starters. Walk every commented YAML key in the
    shipped file and verify it's a recognized schema field in its
    section. Catches the kind of bug where a future edit puts an
    `epitopes`-only key under `vaccine_peptides`.
    """
    from importlib.resources import files
    import re
    from vaxrank.config.schema import (
        EpitopesConfigSchema,
        ManufacturabilityConfigSchema,
        VaccinePeptidesConfigSchema,
    )

    text = files("vaxrank.config").joinpath("default.yaml").read_text()

    schema_for_section = {
        "epitopes": EpitopesConfigSchema,
        "vaccine_peptides": VaccinePeptidesConfigSchema,
        "manufacturability": ManufacturabilityConfigSchema,
    }
    section_fields = {
        s: set(c.__struct_fields__) for s, c in schema_for_section.items()
    }

    section = None
    # `# - foo` rule-list comments, `# foo: ...` field comments — capture
    # both. Skip lines under a nested `manufacturability:` block inside
    # vaccine_peptides (handled separately under its own section).
    inside_nested_manuf = False
    for raw in text.splitlines():
        line = raw.rstrip()
        # Top-level section header
        m = re.match(r"^([a-z_]+):\s*$", line)
        if m and m.group(1) in schema_for_section:
            section = m.group(1)
            inside_nested_manuf = False
            continue
        if section is None:
            continue
        # Nested manufacturability inside vaccine_peptides
        if section == "vaccine_peptides" and re.match(
                r"^\s+manufacturability:\s*$", line):
            inside_nested_manuf = True
            continue
        if inside_nested_manuf:
            continue
        # Match commented field lines: `  # field: value`
        m = re.match(r"^\s*#\s*([a-z_]+):\s", line)
        if not m:
            continue
        key = m.group(1)
        # Skip rule-list entries (those are values, not field names).
        # Heuristic: rule-list lines start with `# - ` not `# key:`.
        # Field comments always have a colon AFTER the name and at least
        # one char before the # (the indent).
        assert key in section_fields[section], (
            f"Commented example `{key}:` in section `{section}` is not a "
            f"valid field (allowed: {sorted(section_fields[section])}). "
            f"Either move the line to the right section or fix the key name."
        )


def test_print_default_config_cli_dumps_yaml(capsys):
    """`vaxrank --print-default-config` should write the bundled YAML and
    exit with code 0 before doing any real work."""
    from vaxrank.cli.entry_point import main
    with pytest.raises(SystemExit) as excinfo:
        main(["--print-default-config"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert captured.out.startswith("# Vaxrank default configuration")
    assert "epitopes:" in captured.out
    assert "vaccine_peptides:" in captured.out
    assert "manufacturability:" in captured.out
