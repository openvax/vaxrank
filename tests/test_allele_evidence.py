# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Attribution of peptide-level evidence to MHC alleles."""

import pytest
from mhctools.pred import Prediction

from vaxrank.allele_evidence import (
    ALLELE_SOURCE_BROADCAST,
    ALLELE_SOURCE_FROM_INPUT,
    ALLELE_SOURCE_SELECTED,
    AllelePolicy,
    AllelePolicyError,
    attribute_alleles,
)
from vaxrank.candidate_epitope import CandidateEpitope


def _affinity(allele, ic50, predictor="mhcflurry"):
    return Prediction(
        kind="pMHC_affinity", predictor_name=predictor, predictor_version="1",
        allele=allele, peptide="SIINFEKL", value=ic50, score=0.0)


def _presentation(allele, score, predictor="mhcflurry"):
    return Prediction(
        kind="pMHC_presentation", predictor_name=predictor,
        predictor_version="1", allele=allele, peptide="SIINFEKL",
        value=None, score=score)


def _processing(score=0.7):
    return Prediction(
        kind="antigen_processing", predictor_name="mhcflurry",
        predictor_version="1", allele="", peptide="SIINFEKL",
        value=None, score=score)


def _epitope(*predictions, patient_alleles=()):
    return CandidateEpitope(
        sequence="SIINFEKL", source_sequence="SIINFEKL", offset=0,
        prediction_id="p", patient_alleles=tuple(patient_alleles),
        predictions=tuple(predictions))


def test_selection_ranks_within_one_predictor():
    """A recorded (predictor, value) pair must be one the predictor emitted.

    Taking a per-allele best across predictors compares two independent
    scales *and* records a value the named predictor never produced — which
    makes the attribution unreconstructable, the one thing it exists for.
    """
    epitope = _epitope(
        _affinity("HLA-A*02:01", 9000.0, "mhcflurry"),
        _affinity("HLA-A*02:01", 5.0, "netmhcpan"),
        _affinity("HLA-B*07:02", 100.0, "mhcflurry"),
        _affinity("HLA-B*07:02", 8000.0, "netmhcpan"),
        _processing(),
        patient_alleles=("HLA-A*02:01", "HLA-B*07:02"))
    genotype = ["HLA-A*02:01", "HLA-B*07:02"]

    [chosen] = attribute_alleles(
        epitope, genotype, AllelePolicy.parse("selected"))
    # mhcflurry is the canonical default, so its numbers rank: B (100) beats
    # A (9000). The cross-predictor minimum would have picked A on 5.0 and
    # attributed that value to mhcflurry, which never predicted it.
    assert chosen.allele == "HLA-B*07:02"
    assert chosen.rank_predictor == "mhcflurry"
    assert chosen.rank_value == pytest.approx(100.0)

    # The value recorded is one the named predictor actually emitted.
    emitted = {
        (p.predictor_name, p.allele): p.value
        for p in epitope.predictions_flat() if p.kind == "pMHC_affinity"}
    assert emitted[(chosen.rank_predictor, chosen.allele)] == pytest.approx(
        chosen.rank_value)


def test_selection_honors_default_methods():
    """The config's chosen model decides, not the most optimistic one."""
    epitope = _epitope(
        _affinity("HLA-A*02:01", 9000.0, "mhcflurry"),
        _affinity("HLA-A*02:01", 5.0, "netmhcpan"),
        _affinity("HLA-B*07:02", 100.0, "mhcflurry"),
        _affinity("HLA-B*07:02", 8000.0, "netmhcpan"),
        _processing(),
        patient_alleles=("HLA-A*02:01", "HLA-B*07:02"))

    [chosen] = attribute_alleles(
        epitope, ["HLA-A*02:01", "HLA-B*07:02"],
        AllelePolicy.parse("selected"),
        default_methods={"pMHC_affinity": "netmhcpan"})
    assert (chosen.allele, chosen.rank_predictor) == (
        "HLA-A*02:01", "netmhcpan")
    assert chosen.rank_value == pytest.approx(5.0)


def test_auto_axis_picks_the_axis_that_covers_the_alleles():
    """Partial presentation columns must not drop the unranked alleles.

    A file carrying presentation for one allele and affinity for all of them
    should rank on affinity. Short-circuiting on the first non-empty axis
    silently reduced a ``top:3`` to a single allele.
    """
    epitope = _epitope(
        _presentation("HLA-A*02:01", 0.4),
        _affinity("HLA-A*02:01", 900.0, "netmhcpan"),
        _affinity("HLA-B*07:02", 10.0, "netmhcpan"),
        _affinity("HLA-C*07:01", 20.0, "netmhcpan"),
        _processing(),
        patient_alleles=("HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"))

    got = attribute_alleles(
        epitope, ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"],
        AllelePolicy.parse("top:3"))
    assert [a.allele for a in got] == [
        "HLA-B*07:02", "HLA-C*07:01", "HLA-A*02:01"]
    assert {a.rank_kind for a in got} == {"pMHC_affinity"}


def test_from_input_reflects_the_input_not_the_growing_allele_set():
    """Attribution must be idempotent across repeated scoring passes.

    ``CandidateEpitope.__post_init__`` folds every scored allele back into
    ``patient_alleles``, so deriving "the input paired these" from that field
    made a broadcast allele claim, on a second pass, that the source file had
    named it.
    """
    epitope = _epitope(
        _affinity("HLA-A*02:01", 50.0), _processing(),
        patient_alleles=("HLA-A*02:01",))
    genotype = ["HLA-A*02:01", "HLA-C*07:01"]

    first = attribute_alleles(epitope, genotype, AllelePolicy.parse("all"))
    by_allele = {a.allele: a.source for a in first}
    assert by_allele == {
        "HLA-A*02:01": ALLELE_SOURCE_FROM_INPUT,
        "HLA-C*07:01": ALLELE_SOURCE_BROADCAST,
    }

    # Simulate the widened patient_alleles a scoring pass produces.
    widened = _epitope(
        _affinity("HLA-A*02:01", 50.0), _processing(),
        patient_alleles=("HLA-A*02:01", "HLA-C*07:01"))
    again = attribute_alleles(widened, genotype, AllelePolicy.parse("all"))
    assert {a.allele: a.source for a in again} == by_allele


def test_explicit_genotype_can_narrow_not_only_widen():
    """A patient's real typing must be able to correct a wider input panel."""
    epitope = _epitope(
        _affinity("HLA-A*02:01", 50.0),
        _affinity("HLA-B*07:02", 60.0),
        _processing(),
        patient_alleles=("HLA-A*02:01", "HLA-B*07:02"))

    got = attribute_alleles(
        epitope, ["HLA-A*02:01"], AllelePolicy.parse("all"))
    assert [a.allele for a in got] == ["HLA-A*02:01"]


def test_only_peptide_level_kinds_are_attributed():
    """A malformed blank-allele row must not become per-allele evidence."""
    from vaxrank.epitope_dsl import epitopes_to_topiary_df

    # An affinity prediction that lost its allele is malformed, not
    # peptide-level. It is neither attributed nor carried through: see
    # test_malformed_allele_scoped_predictions_never_reach_the_frame for why
    # leaving it in the frame is unsafe.
    epitope = _epitope(
        Prediction(kind="pMHC_affinity", predictor_name="mhcflurry",
                   predictor_version="1", allele="", peptide="SIINFEKL",
                   value=50.0, score=0.0),
        patient_alleles=("HLA-A*02:01", "HLA-B*07:02"))
    frame = epitopes_to_topiary_df(
        [epitope], policy=AllelePolicy.parse("all"))
    assert frame.empty


def test_malformed_policies_are_rejected_where_they_are_written():
    """Typos fail at config time, not silently mid-run."""
    with pytest.raises(AllelePolicyError, match="Unknown allele policy"):
        AllelePolicy.parse("typo")
    with pytest.raises(AllelePolicyError, match="at least one allele"):
        AllelePolicy.parse("top:0")
    with pytest.raises(AllelePolicyError, match="integer N"):
        AllelePolicy.parse("top:many")
    # The axis is validated too — it used to escape until a selection ran,
    # and for "all" / "from_input" it never ran at all.
    with pytest.raises(AllelePolicyError, match="selection axis"):
        AllelePolicy.parse("selected", axis="affinity")
    with pytest.raises(AllelePolicyError, match="selection axis"):
        AllelePolicy.parse("all", axis="nonsense")

    # A hand-built policy that never went through parse is refused rather
    # than silently treated as top-N.
    epitope = _epitope(_processing(), patient_alleles=("HLA-A*02:01",))
    with pytest.raises(AllelePolicyError, match="Unknown allele policy name"):
        attribute_alleles(epitope, ["HLA-A*02:01"], AllelePolicy(name="typo"))


def test_attribution_requires_reconstructable_provenance():
    """A selection without its ranking axis could not be reproduced."""
    from vaxrank.allele_evidence import AlleleAttribution

    with pytest.raises(ValueError, match="must record the axis"):
        AlleleAttribution(allele="HLA-A*02:01", source=ALLELE_SOURCE_SELECTED)
    with pytest.raises(ValueError, match="Only a selected allele"):
        AlleleAttribution(
            allele="HLA-A*02:01", source=ALLELE_SOURCE_BROADCAST,
            rank_kind="pMHC_affinity")


def test_predictor_is_chosen_by_coverage_unless_the_config_names_one():
    """A predictor covering one allele must not silently reduce a top:3.

    Preference order alone picked the canonical name even when it could rank
    a single allele out of six — the same silent reduction as choosing the
    axis by presence rather than coverage, arriving through a second door.
    """
    genotype = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"]
    epitope = _epitope(
        _affinity("HLA-A*02:01", 50.0, "mhcflurry"),      # covers 1
        _affinity("HLA-A*02:01", 900.0, "netmhcpan"),     # covers 3
        _affinity("HLA-B*07:02", 10.0, "netmhcpan"),
        _affinity("HLA-C*07:01", 20.0, "netmhcpan"),
        _processing(),
        patient_alleles=tuple(genotype))

    got = attribute_alleles(epitope, genotype, AllelePolicy.parse("top:3"))
    assert [a.allele for a in got] == [
        "HLA-B*07:02", "HLA-C*07:01", "HLA-A*02:01"]
    assert {a.rank_predictor for a in got} == {"netmhcpan"}

    # An explicit choice is the user's, and still wins over coverage.
    named = attribute_alleles(
        epitope, genotype, AllelePolicy.parse("top:3"),
        default_methods={"pMHC_affinity": "mhcflurry"})
    assert [(a.allele, a.rank_predictor) for a in named] == [
        ("HLA-A*02:01", "mhcflurry")]


def test_recorded_values_survive_serialization():
    """Predictions from a pandas frame carry numpy scalars.

    A recorded value has to round-trip to stay reconstructable, and the
    serializer refuses numpy types outright, so the coercion happens where
    the value is captured.
    """
    import numpy as np

    from vaxrank.allele_evidence import AlleleAttribution

    epitope = _epitope(
        _affinity("HLA-A*02:01", np.float64(50.0)), _processing(),
        patient_alleles=("HLA-A*02:01",))
    [attribution] = attribute_alleles(
        epitope, ["HLA-A*02:01"], AllelePolicy.parse("selected"))

    assert type(attribution.rank_value) is float
    assert AlleleAttribution.from_json(attribution.to_json()) == attribution


def test_mhc_dependence_table_covers_every_topiary_kind():
    """A new prediction kind must be classified, not silently un-attributable.

    "Is this prediction allele-scoped?" used to be re-derived inline at four
    call sites, two of them testing only for a blank allele — which is also
    what a malformed row looks like. The table is the single answer, and it
    is checked against topiary at import.
    """
    from topiary.ranking import KIND_ALIASES

    from vaxrank.allele_evidence import (
        ALLELE_SCOPED_KINDS, PEPTIDE_LEVEL_KINDS,
    )

    classified = PEPTIDE_LEVEL_KINDS | ALLELE_SCOPED_KINDS
    unclassified = sorted(set(KIND_ALIASES.values()) - classified)
    # CI is where this drift should surface. At runtime an unclassified kind
    # is treated as not peptide-level (never projected, so nothing is
    # fabricated) and warned about, rather than refusing to import.
    assert not unclassified, (
        "topiary knows prediction kind(s) %s that vaxrank.allele_evidence "
        "does not classify" % ", ".join(unclassified))
    # The two sets are a partition, not overlapping guesses.
    assert not PEPTIDE_LEVEL_KINDS & ALLELE_SCOPED_KINDS


def test_peptide_level_requires_both_the_kind_and_a_blank_allele():
    """Neither half of the predicate is sufficient alone."""
    from vaxrank.allele_evidence import is_allele_scoped, is_peptide_level

    # Blank allele, but an allele-scoped kind: malformed, not peptide-level.
    malformed = Prediction(
        kind="pMHC_affinity", predictor_name="mhcflurry",
        predictor_version="1", allele="", peptide="SIINFEKL",
        value=50.0, score=0.0)
    assert not is_peptide_level(malformed)
    assert not is_allele_scoped(malformed)

    # Peptide-level kind that already carries an allele: already attributed.
    attributed = Prediction(
        kind="antigen_processing", predictor_name="mhcflurry",
        predictor_version="1", allele="HLA-A*02:01", peptide="SIINFEKL",
        value=None, score=0.7)
    assert not is_peptide_level(attributed)

    assert is_peptide_level(_processing())
    assert is_allele_scoped(_affinity("HLA-A*02:01", 50.0))


def test_malformed_allele_scoped_predictions_never_reach_the_frame(caplog):
    """A blank-allele affinity row must not become per-allele evidence.

    Topiary classifies a kind's MHC dependence by scanning the rows it is
    given, so a candidate whose only affinity leaves carry no allele reads as
    peptide-level and has its value projected across every allele group —
    inventing a binding prediction for alleles the model never scored.
    Verified live against topiary 5.21.0; openvax/topiary#197 fixes it
    upstream, but vaxrank supports the whole >=5.17.1 range.
    """
    import logging

    from topiary.ranking import EvalContext, parse

    from vaxrank.epitope_dsl import (
        PREDICTION_GROUP_COLUMNS, epitopes_to_topiary_df,
    )

    def stability(allele, value):
        return Prediction(
            kind="pMHC_stability", predictor_name="netmhcstabpan",
            predictor_version="1", allele=allele, peptide="SIINFEKL",
            value=value, score=0.0)

    malformed = Prediction(
        kind="pMHC_affinity", predictor_name="mhcflurry", predictor_version="1",
        allele="", peptide="SIINFEKL", value=7.0, score=0.0)
    epitope = _epitope(
        malformed,
        stability("HLA-A*02:01", 100.0), stability("HLA-B*07:02", 200.0),
        patient_alleles=("HLA-A*02:01", "HLA-B*07:02"))

    with caplog.at_level(logging.WARNING):
        frame = epitopes_to_topiary_df([epitope])

    # The malformed row is absent, and the drop is reported with enough to
    # find it — a silent drop would be its own defect.
    assert set(frame["kind"]) == {"pMHC_stability"}
    warning = "\n".join(
        record.getMessage() for record in caplog.records
        if record.levelno >= logging.WARNING)
    assert "no allele" in warning
    assert "SIINFEKL" in warning and "pMHC_affinity" in warning

    # And nothing invents an affinity for alleles that were never scored.
    context = EvalContext(frame, group_keys=list(PREDICTION_GROUP_COLUMNS))
    values = parse("affinity.value").eval(context).reindex(
        context.group_index).tolist()
    assert all(value != value for value in values)  # all NaN, not 7.0
