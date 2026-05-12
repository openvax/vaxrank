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
from mhctools.pred import Prediction

from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
from vaxrank.manufacturability import DEFAULT_MANUFACTURABILITY_RULES
from vaxrank.candidate_epitope import CandidateEpitope, Peptide
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


def _make_mutant_epitope(ic50=10.0, peptide="A" * 9, source="A" * 25,
                         offset=0):
    """One ``CandidateEpitope`` with a single pMHC_affinity prediction, marked
    as overlapping the mutation."""
    pred = Prediction(
        kind='pMHC_affinity',
        predictor_name='test',
        predictor_version='',
        allele='HLA-A*02:01',
        peptide=peptide,
        value=ic50,
        score=0.0,
        percentile_rank=0.5,
    )
    return CandidateEpitope(
        mutant=Peptide(
            sequence=peptide,
            source=source,
            offset=offset,
            predictions=(pred,),
        ),
        overlaps_mutation=True,
        occurs_in_reference=False,
    )


def test_default_manufacturability_tuple_length():
    """Default sort tuple length = 10 (legacy behavior, unchanged)."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitopes=[_make_mutant_epitope()],
    )
    tup = vp.peptide_synthesis_difficulty_score_tuple()
    assert len(tup) == 10
    assert len(DEFAULT_MANUFACTURABILITY_RULES) == 10


def test_custom_rules_shorten_tuple():
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitopes=[_make_mutant_epitope()],
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
        epitopes=[_make_mutant_epitope()],
    )
    # Legacy: sqrt(9) * mutant_epitope_score
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(np.sqrt(9) * vp.mutant_epitope_score)


def test_combined_score_reads_times_epitope_mode():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=[_make_mutant_epitope()],
        combined_score_mode="reads_times_epitope",
    )
    # expression_score is an honest metric; mode only affects combined_score
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(9.0 * vp.mutant_epitope_score)


def test_combined_score_epitope_only_mode():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=[_make_mutant_epitope()],
        combined_score_mode="epitope_only",
    )
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(vp.mutant_epitope_score)


def test_vaccine_peptide_rejects_unknown_combined_score_mode():
    with pytest.raises(ValueError, match="combined_score_mode"):
        VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitopes=[_make_mutant_epitope()],
            combined_score_mode="garbage",
        )


def test_end_to_end_yaml_rules_through_manufacturability_config_to_peptide():
    """YAML with peptide.manufacturability.rules loads through
    load_vaxrank_config, then manufacturability_config_from_args,
    producing a tuple on ManufacturabilityConfig. A VaccinePeptide
    built with those rules gets a sort tuple of the expected
    length. Post-2.20: manufacturability lives on its own struct,
    not on VaccineConfig."""
    from vaxrank.cli.vaccine_config_args import (
        manufacturability_config_from_args)
    yaml_content = """
peptide:
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
        mfg = manufacturability_config_from_args(args)
        vc = vaccine_config_from_args(args)
        assert mfg.rules == (
            "cysteine_count", "cterm_hydropathy", "aspartate_proline",
        )
        assert isinstance(mfg.rules, tuple)
        assert vc.combined_score_mode == "epitope_only"

        vp = VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitopes=[_make_mutant_epitope()],
            manufacturability_rules=mfg.rules,
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
        epitopes=[_make_mutant_epitope()],
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
        epitopes=[_make_mutant_epitope()],
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
        epitopes=[_make_mutant_epitope()],
        ranking_rules=["mutant_epitope_score", "n_alt_reads"],
    )
    assert vp.ranking_rules == ("mutant_epitope_score", "n_alt_reads")
    assert isinstance(vp.ranking_rules, tuple)
    # Sort tuple length reflects the custom 2-rule list.
    assert len(vp.lexicographic_sort_key()) == 2


def test_vaccine_peptide_custom_ranking_rules_in_to_dict():
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitopes=[_make_mutant_epitope()],
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
            epitopes=[_make_mutant_epitope()],
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
    # Top-level sections present (post-2.19 layout: peptide /
    # mrna at top level; manufacturability nests under peptide).
    assert set(parsed) >= {"epitopes", "vaccine_peptides", "peptide", "mrna"}
    assert "manufacturability" in parsed["peptide"]


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
    from vaxrank.cli.vaccine_config_args import (
        manufacturability_config_from_args, vaccine_config_from_args)
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.manufacturability import DEFAULT_MANUFACTURABILITY_RULES
    from vaxrank.manufacturability_config import ManufacturabilityConfig
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
    mfg = manufacturability_config_from_args(args)
    # CandidateEpitope side has no list-typed fields, so == works.
    assert ec == EpitopeConfig()
    # Vaccine side: ranking_rules tuple matches the default.
    assert vc.ranking_rules == DEFAULT_RANKING_RULES
    # Manufacturability rules now live on ManufacturabilityConfig.
    assert mfg.rules == DEFAULT_MANUFACTURABILITY_RULES
    # Vaccine scalar fields match defaults.
    vc_defaults = VaccineConfig()
    for field in (
        "preferred_peptide_length", "min_peptide_length",
        "max_peptide_length", "padding_around_mutation",
        "max_vaccine_peptides_per_variant", "num_mutant_epitopes_to_keep",
        "score_fraction_of_best", "combined_score_mode",
        "require_mutant_epitopes_in_variant",
    ):
        assert getattr(vc, field) == getattr(vc_defaults, field), (
            f"default.yaml VaccineConfig drift: {field}")
    # Manufacturability scalar fields match defaults.
    mfg_defaults = ManufacturabilityConfig()
    for field in (
        "max_c_terminal_hydropathy", "min_kmer_hydropathy",
        "max_kmer_hydropathy_low_priority",
        "max_kmer_hydropathy_high_priority",
    ):
        assert getattr(mfg, field) == getattr(mfg_defaults, field), (
            f"default.yaml ManufacturabilityConfig drift: {field}")


def test_default_yaml_uncommented_keys_match_schema():
    """The shipped default.yaml must validate against
    ``VaxrankConfigSchema`` end-to-end. ``forbid_unknown_fields=True``
    on each sub-schema catches drift like a key placed in the wrong
    section or renamed without updating the YAML."""
    from importlib.resources import files
    import yaml
    import msgspec
    from vaxrank.config.schema import VaxrankConfigSchema

    text = files("vaxrank.config").joinpath("default.yaml").read_text()
    parsed = yaml.safe_load(text)
    msgspec.convert(parsed, VaxrankConfigSchema)


def test_default_yaml_commented_examples_are_valid():
    """Commented-out example lines under ``epitopes:`` in default.yaml
    are advertised as 'just uncomment to enable' starters. Walk every
    commented field name under that section and verify it's a
    recognized schema field. Catches drift like a key renamed in the
    schema but left dangling in the YAML.

    Scope is intentionally narrow: only ``epitopes:`` carries
    commented-out scalar examples in default.yaml; the
    ``vaccine_peptides:`` / ``peptide:`` / ``mrna:`` sections ship
    fully uncommented and are validated by
    ``test_default_yaml_uncommented_keys_match_schema`` above.
    """
    from importlib.resources import files
    import re
    from vaxrank.config.schema import EpitopesConfigSchema

    text = files("vaxrank.config").joinpath("default.yaml").read_text()
    epitopes_fields = set(EpitopesConfigSchema.__struct_fields__)

    section = None
    for raw in text.splitlines():
        line = raw.rstrip()
        m = re.match(r"^([a-z_]+):\s*$", line)
        if m:
            section = m.group(1)
            continue
        if section != "epitopes":
            continue
        # Match top-level commented field lines: `  # field: value`.
        # Require exactly one space between `#` and the field name so
        # that deeper-nested commented values
        # (`#   pMHC_affinity: mhcflurry` under a parent
        # `# default_methods:`) are skipped — they're values of the
        # parent dict, not section-level fields.
        m = re.match(r"^\s*# ([a-z_]+):\s", line)
        if not m:
            continue
        key = m.group(1)
        assert key in epitopes_fields, (
            f"Commented example `{key}:` in section `epitopes` is not a "
            f"valid field (allowed: {sorted(epitopes_fields)}). "
            f"Either move the line to the right section or fix the key name."
        )


def test_default_yaml_default_methods_block_uncomments_cleanly():
    """``epitopes.default_methods`` is the only nested-dict example in
    default.yaml that ships commented out. Strip the leading ``# `` from
    that block and confirm the result parses + validates against the
    schema, so the documented "uncomment to enable" pattern actually
    works for that block."""
    from importlib.resources import files
    import re
    import yaml
    import msgspec
    from vaxrank.config.schema import VaxrankConfigSchema

    text = files("vaxrank.config").joinpath("default.yaml").read_text()
    activated = re.sub(
        r"^(\s*)# (default_methods:|  pMHC_\w+:.*)$",
        r"\1\2",
        text, flags=re.MULTILINE,
    )
    parsed = yaml.safe_load(activated)
    msgspec.convert(parsed, VaxrankConfigSchema)
    assert (parsed["epitopes"]["default_methods"]["pMHC_affinity"]
            == "mhcflurry")
    assert (parsed["epitopes"]["default_methods"]["pMHC_stability"]
            == "netmhcstabpan")
    assert (parsed["epitopes"]["default_methods"]["pMHC_presentation"]
            == "mhcflurry")


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
