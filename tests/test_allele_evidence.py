# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for frame-level allele attribution."""

import pandas as pd
import pytest

from vaxrank.allele_evidence import (
    ALLELE_SOURCE_BROADCAST,
    ALLELE_SOURCE_FROM_INPUT,
    ALLELE_SOURCE_SELECTED,
    AllelePolicy,
    AllelePolicyError,
    attribute_frame,
    is_peptide_level_kind,
)
from vaxrank.epitope_dsl import (
    PIPELINE_GROUP_COLUMNS,
    PREDICTION_GROUP_COLUMNS,
)

A = "HLA-A*02:01"
B = "HLA-B*07:02"


def _frame(identity_column, *, processing=True):
    """One peptide: strong affinity for A, weak for B, plus a processing row."""
    rows = [
        {identity_column: "s", "peptide": "SIINFEKL", "peptide_offset": 0,
         "allele": allele, "kind": "pMHC_affinity",
         "prediction_method_name": "mhcflurry", "predictor_version": "2.1.1",
         "value": affinity, "score": None, "percentile_rank": None}
        for allele, affinity in ((A, 50.0), (B, 900.0))]
    if processing:
        rows.append(
            {identity_column: "s", "peptide": "SIINFEKL", "peptide_offset": 0,
             "allele": "", "kind": "antigen_processing",
             "prediction_method_name": "mhcflurry",
             "predictor_version": "2.1.1", "value": None, "score": 0.8,
             "percentile_rank": None})
    return pd.DataFrame(rows)


FRAME_SHAPES = [
    pytest.param("prediction_id", list(PREDICTION_GROUP_COLUMNS), id="external"),
    pytest.param(
        "source_sequence_name", list(PIPELINE_GROUP_COLUMNS), id="native"),
]


@pytest.mark.parametrize("identity_column,group_columns", FRAME_SHAPES)
@pytest.mark.parametrize("policy_text,expected", [
    ("all", [A, B]),
    ("selected", [A]),
    ("top:2", [A, B]),
    ("top:1", [A]),
    ("from_input", [A, B]),
])
def test_both_input_paths_attribute_identically(
        identity_column, group_columns, policy_text, expected):
    """The same policy must produce the same answer on either frame shape.

    This is the point of attributing at the frame layer. The previous
    implementation ran on CandidateEpitope objects, which only the external
    report path builds before scoring, so the setting applied to LENS runs
    and silently did nothing on the native pipeline — the path most runs
    take (openvax/vaxrank#349).
    """
    policy = AllelePolicy.parse(policy_text)
    attributions = attribute_frame(
        _frame(identity_column), policy, group_columns=group_columns)

    [attributed] = attributions.values()
    assert [a.allele for a in attributed] == expected


@pytest.mark.parametrize("identity_column,group_columns", FRAME_SHAPES)
def test_a_selection_records_how_it_was_made(identity_column, group_columns):
    """A selected allele carries the evidence that selected it.

    Without this the attribution is an assertion the reader cannot check.
    """
    attributions = attribute_frame(
        _frame(identity_column), AllelePolicy.parse("selected"),
        group_columns=group_columns)

    [[attributed]] = attributions.values()
    assert attributed.allele == A
    assert attributed.source == ALLELE_SOURCE_SELECTED
    assert attributed.rank_kind == "pMHC_affinity"
    assert attributed.rank_predictor == "mhcflurry"
    assert attributed.rank_value == 50.0
    assert attributed.rank_position == 1
    # Plain builtins, not numpy scalars: these are serialized, and a
    # numpy.int64 fails with "Cannot convert 1 : numpy.int64".
    assert type(attributed.rank_value) is float
    assert type(attributed.rank_position) is int


@pytest.mark.parametrize("identity_column,group_columns", FRAME_SHAPES)
def test_a_peptide_with_nothing_to_rank_keeps_its_whole_genotype(
        identity_column, group_columns):
    """No allele-scoped evidence means no winner, not an invented one.

    "We do not know which allele presents this" is not the same claim as
    "this allele does", so the fallback is the genotype rather than a pick.
    """
    frame = pd.DataFrame([
        {identity_column: "s", "peptide": "SIINFEKL", "peptide_offset": 0,
         "allele": "", "kind": "antigen_processing",
         "prediction_method_name": "mhcflurry", "predictor_version": "2.1.1",
         "value": None, "score": 0.8, "percentile_rank": None}])

    attributions = attribute_frame(
        frame, AllelePolicy.parse("selected"),
        genotype_for=lambda key: (A, B), group_columns=group_columns)

    [attributed] = attributions.values()
    assert [a.allele for a in attributed] == [A, B]
    assert {a.source for a in attributed} == {ALLELE_SOURCE_BROADCAST}


@pytest.mark.parametrize("identity_column,group_columns", FRAME_SHAPES)
def test_a_peptide_with_no_allele_free_evidence_is_not_attributed(
        identity_column, group_columns):
    """Attribution says where allele-free evidence counts.

    A peptide that has none has nothing to attribute, and inventing an entry
    would put a claim in the audit record that no evidence supports.
    """
    attributions = attribute_frame(
        _frame(identity_column, processing=False),
        AllelePolicy.parse("selected"), group_columns=group_columns)

    assert attributions == {}


@pytest.mark.parametrize("identity_column,group_columns", FRAME_SHAPES)
def test_an_explicit_genotype_can_narrow_as_well_as_widen(
        identity_column, group_columns):
    """The patient's typing wins over what the input file happened to pair.

    Unioning the input's alleles in would mean an explicit genotype could
    only ever add alleles, so it could not correct an input produced with
    the wrong or a wider panel.
    """
    attributions = attribute_frame(
        _frame(identity_column), AllelePolicy.parse("all"),
        genotype_for=lambda key: (A,), group_columns=group_columns)

    [attributed] = attributions.values()
    assert [a.allele for a in attributed] == [A]
    assert attributed[0].source == ALLELE_SOURCE_FROM_INPUT


def test_the_peptide_level_partition_comes_from_topiary():
    """No local copy of the MHC-dependence table (openvax/vaxrank#357)."""
    import topiary

    for kind, dependence in topiary.KIND_MHC_DEPENDENCE.items():
        assert is_peptide_level_kind(kind) is (dependence == "none")
    # An unclassified kind is never peptide-level, so a topiary release
    # adding one cannot make vaxrank credit an allele with it.
    assert is_peptide_level_kind("a_kind_topiary_has_not_defined") is False


@pytest.mark.parametrize("bad", ["nonsense", "top:0", "top:x", "top:-1"])
def test_an_unusable_policy_fails_at_config_time(bad):
    """A typo must fail when the config is read, not mid-run."""
    from vaxrank.epitope_config import EpitopeConfig

    with pytest.raises(AllelePolicyError):
        EpitopeConfig(allele_free_evidence=bad)


def test_the_policy_is_reachable_from_yaml_and_the_cli():
    """A knob nobody can set is not a knob.

    Review of the first attempt found the field settable only by
    constructing EpitopeConfig in Python — neither the YAML schema nor the
    CLI exposed it.
    """
    import argparse

    from vaxrank.cli.epitope_config_args import (
        add_epitope_prediction_args, epitope_config_from_args)
    from vaxrank.config.loader import _EPITOPE_CONFIG_MAPPING

    parser = argparse.ArgumentParser()
    add_epitope_prediction_args(parser)
    args = parser.parse_args(["--allele-free-evidence", "top:2"])
    assert epitope_config_from_args(args).allele_policy.limit == 2

    assert ("epitopes.allele_free_evidence",
            "allele_free_evidence") in _EPITOPE_CONFIG_MAPPING


def test_attributions_are_recorded_on_the_candidate(tmp_path):
    """The computed attribution has to land somewhere or it is waste.

    Review of the first cut found it computed on every scoring pass and
    used only to decide whether to warn — a full pass over the frame whose
    result was discarded. It is recorded on the candidate now, which is
    both what makes it observable and what the report layer reads.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_lens_report

    path = tmp_path / "processing.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t0.8\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t900\t0.8\n")

    [selected] = read_lens_report(
        path,
        epitope_config=EpitopeConfig(allele_free_evidence="selected")).epitopes
    [every] = read_lens_report(
        path,
        epitope_config=EpitopeConfig(allele_free_evidence="all")).epitopes

    assert [a.allele for a in selected.allele_attributions] == [A]
    assert selected.allele_attributions[0].source == ALLELE_SOURCE_SELECTED
    assert [a.allele for a in every.allele_attributions] == [A, B]

    # Per-allele scores are deliberately unchanged between the two: topiary
    # projects a peptide-level value across every allele group keyed on the
    # kind (openvax/topiary#232), so a narrowing policy records the credited
    # allele without yet narrowing the score. The warning says so; this pins
    # that the two really do still agree, so the claim is not stale.
    assert selected.per_allele_scores == every.per_allele_scores


def test_a_frame_with_no_allele_free_evidence_is_not_walked():
    """The common frame must not pay for a feature it cannot use.

    Turning a large frame into dicts to discover there is nothing to
    attribute cost real time on every scoring pass — 0.4s on a 200k-row
    frame, for an empty result.
    """
    frame = _frame("prediction_id", processing=False)
    calls = []

    class CountingFrame(type(frame)):
        def to_dict(self, *args, **kwargs):
            calls.append(1)
            return super().to_dict(*args, **kwargs)

    counting = CountingFrame(frame)
    assert attribute_frame(
        counting, AllelePolicy.parse("selected"),
        group_columns=list(PREDICTION_GROUP_COLUMNS)) == {}
    assert calls == []


@pytest.mark.parametrize("policy_text", ["all", "selected", "top:2", "from_input"])
def test_attributions_survive_a_native_round_trip(tmp_path, policy_text):
    """A saved run must reload with the attribution it actually made.

    The attribution depends on the policy in force when the candidate was
    scored, and ``load_predictions`` takes no config — so recomputing on
    load would silently answer a different question and present it as what
    the run produced. Recording it is what makes a finished ranking
    explainable later.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import (
        load_predictions, read_lens_report, save_predictions)

    source = tmp_path / "processing.tsv"
    source.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t0.8\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t900\t0.8\n")

    [scored] = read_lens_report(
        source,
        epitope_config=EpitopeConfig(
            allele_free_evidence=policy_text)).epitopes
    assert scored.allele_attributions  # the test would pass vacuously otherwise

    path = tmp_path / "native.csv"
    save_predictions([scored], str(path))
    [reloaded] = load_predictions(str(path))

    assert reloaded.allele_attributions == scored.allele_attributions


def test_a_reload_does_not_re_derive_under_the_default_policy(tmp_path):
    """The recorded answer wins over what the loader would compute.

    A run scored under 'selected' credits one allele. Reloading it must not
    quietly produce the two-allele answer the default policy would give,
    which is what recomputation on load would do.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import (
        load_predictions, read_lens_report, save_predictions)

    source = tmp_path / "processing.tsv"
    source.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t0.8\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t900\t0.8\n")

    [scored] = read_lens_report(
        source,
        epitope_config=EpitopeConfig(allele_free_evidence="selected")).epitopes
    path = tmp_path / "native.csv"
    save_predictions([scored], str(path))

    [reloaded] = load_predictions(str(path))

    assert [a.allele for a in reloaded.allele_attributions] == [A]
    assert reloaded.allele_attributions[0].source == ALLELE_SOURCE_SELECTED
    # And the ranking provenance is intact, not just the allele name.
    assert reloaded.allele_attributions[0].rank_predictor == "mhcflurry"
    assert reloaded.allele_attributions[0].rank_value == 50.0
