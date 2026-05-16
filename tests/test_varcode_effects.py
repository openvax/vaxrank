from varcode import MultiOutcomeEffect
from varcode.effect_candidates import EffectCandidate
from varcode.effects.effect_classes import MutationEffect

from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.report import TemplateDataCreator
from vaxrank.varcode_effects import (
    OUTCOME_SELECTION_HIGHEST_PRIORITY,
    OUTCOME_SELECTION_MOST_LIKELY,
    OUTCOME_SELECTION_MULTI_OUTCOME,
    select_varcode_effect_outcome,
    summarize_varcode_effect_outcomes,
)


class _LikelyEffect(MutationEffect):
    @property
    def short_description(self):
        return "likely"

    @property
    def transcript_name(self):
        return "TX-LIKELY"

    @property
    def transcript_id(self):
        return "ENST-LIKELY"


class _PriorityEffect(MutationEffect):
    @property
    def short_description(self):
        return "priority"

    @property
    def transcript_name(self):
        return "TX-PRIORITY"

    @property
    def transcript_id(self):
        return "ENST-PRIORITY"


class _FakeOutcomeSet(MultiOutcomeEffect):
    def __init__(self, likely_effect, priority_effect):
        MultiOutcomeEffect.__init__(self, variant=None)
        self._likely_effect = likely_effect
        self._priority_effect = priority_effect
        self._candidates = (
            EffectCandidate(likely_effect, source="producer-order"),
            EffectCandidate(priority_effect, source="priority-order"),
        )

    @property
    def candidates(self):
        return self._candidates

    @property
    def most_likely_effect(self):
        return self._likely_effect

    @property
    def highest_priority_effect(self):
        return self._priority_effect


class _Transcript:
    id = "ENST-TEST"


class _Variant:
    def __init__(self, effect):
        self.effect = effect

    def __str__(self):
        return "VariantWithOutcomeSet"

    def effect_on_transcript(self, transcript):
        return self.effect


def _outcome_set():
    likely = _LikelyEffect(variant=None)
    priority = _PriorityEffect(variant=None)
    return _FakeOutcomeSet(likely, priority), likely, priority


def test_select_varcode_effect_outcome_distinguishes_selection_modes():
    outcomes, likely, priority = _outcome_set()

    assert select_varcode_effect_outcome(
        outcomes,
        OUTCOME_SELECTION_MOST_LIKELY) is likely
    assert select_varcode_effect_outcome(
        outcomes,
        OUTCOME_SELECTION_HIGHEST_PRIORITY) is priority
    assert select_varcode_effect_outcome(
        outcomes,
        OUTCOME_SELECTION_MULTI_OUTCOME) is outcomes
    assert select_varcode_effect_outcome(likely) is likely


def test_varcode_effect_outcome_summary_includes_candidates():
    outcomes, _likely, _priority = _outcome_set()

    summary = summarize_varcode_effect_outcomes(outcomes)

    assert summary["Outcome set type"] == "_FakeOutcomeSet"
    assert summary["Most likely effect"] == "_LikelyEffect: likely"
    assert summary["Highest priority effect"] == "_PriorityEffect: priority"
    assert "_LikelyEffect: likely (producer-order)" in summary[
        "Candidate effects"]
    assert "_PriorityEffect: priority (priority-order)" in summary[
        "Candidate effects"]


def test_report_effect_data_uses_selected_effect_and_shows_outcomes():
    outcomes, _likely, priority = _outcome_set()
    creator = object.__new__(TemplateDataCreator)

    effect_data = creator._effect_data(outcomes, selected_effect=priority)

    assert effect_data["Effect type"] == "_PriorityEffect"
    assert effect_data["Transcript name"] == "TX-PRIORITY"
    assert effect_data["Transcript ID"] == "ENST-PRIORITY"
    assert effect_data["Outcome set type"] == "_FakeOutcomeSet"
    assert effect_data["Most likely effect"] == "_LikelyEffect: likely"
    assert effect_data["Highest priority effect"] == "_PriorityEffect: priority"


def test_mutant_protein_fragment_predicted_effect_selects_outcome():
    outcomes, likely, priority = _outcome_set()
    frag = MutantProteinFragment(
        variant=_Variant(outcomes),
        gene_name="TEST",
        amino_acids="MTEST",
        mutant_amino_acid_start_offset=0,
        mutant_amino_acid_end_offset=1,
        supporting_reference_transcripts=[_Transcript()],
        n_overlapping_reads=0,
        n_alt_reads=0,
        n_ref_reads=0,
        n_alt_reads_supporting_protein_sequence=0,
    )

    assert frag.predicted_effect() is priority
    assert frag.predicted_effect(
        outcome_selection=OUTCOME_SELECTION_MOST_LIKELY) is likely
    assert frag.predicted_effect(
        outcome_selection=OUTCOME_SELECTION_MULTI_OUTCOME) is outcomes
