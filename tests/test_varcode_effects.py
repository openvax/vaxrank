import pytest

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
    assert summary["Producer-preferred effect"] == "_LikelyEffect: likely"
    assert summary["Highest priority effect"] == "_PriorityEffect: priority"
    assert "_LikelyEffect: likely (producer-order)" in summary[
        "Candidate effects"]
    assert "_PriorityEffect: priority (priority-order)" in summary[
        "Candidate effects"]


def test_report_effect_data_uses_selected_effect_and_shows_outcomes():
    outcomes, _likely, priority = _outcome_set()
    creator = object.__new__(TemplateDataCreator)

    effect_data = creator.effect_data(outcomes, selected_effect=priority)

    assert effect_data["Effect type"] == "_PriorityEffect"
    assert effect_data["Transcript name"] == "TX-PRIORITY"
    assert effect_data["Transcript ID"] == "ENST-PRIORITY"
    assert effect_data["Outcome set type"] == "_FakeOutcomeSet"
    assert effect_data["Producer-preferred effect"] == "_LikelyEffect: likely"
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


@pytest.mark.parametrize("position,ref,alt,preferred_type", [
    (117171168, "G", "T", "Substitution"),
    (117171171, "A", "C", "Intronic"),
])
def test_real_splice_outcomes_survive_dna_fallback_and_reporting(
        human_genome_grch37, position, ref, alt, preferred_type):
    """Exercise Varcode 11's ordinary candidates through both annotation doors.

    The last base of CFTR exon 4 is G in Ensembl 75 (GRCh37). Its G>T
    substitution has a coding consequence if normal splicing persists, and
    a distinct exon-skipping candidate. This checks API composition, not
    which predicted splice mechanism occurs in a patient.
    """
    from varcode import Variant
    from varcode.effects import Deletion
    from vaxrank.gene_pathway_check import GenePathwayCheck

    variant = Variant("7", position, ref, alt, human_genome_grch37)
    fragment = MutantProteinFragment.from_variant_dna(variant, 35, effects=variant.effects())
    assert fragment is not None
    assert fragment.gene_name == "CFTR"

    outcomes = fragment.predicted_effect(OUTCOME_SELECTION_MULTI_OUTCOME)
    assert type(outcomes).__name__ == preferred_type
    likely = fragment.predicted_effect(OUTCOME_SELECTION_MOST_LIKELY)
    priority = fragment.predicted_effect(OUTCOME_SELECTION_HIGHEST_PRIORITY)
    assert type(likely).__name__ == preferred_type
    assert isinstance(priority, Deletion)
    assert likely.mutant_protein_sequence != priority.mutant_protein_sequence
    assert fragment.amino_acids in priority.mutant_protein_sequence
    assert fragment.amino_acids not in likely.mutant_protein_sequence
    # Exon skipping retains a zero-width target at the deletion junction.
    assert fragment.mutant_amino_acid_start_offset == 17
    assert fragment.mutant_amino_acid_end_offset == 17

    summary = summarize_varcode_effect_outcomes(outcomes)
    assert summary["Outcome set type"] == preferred_type
    assert preferred_type in summary["Candidate effects"]
    assert "Deletion" in summary["Candidate effects"]
    assert "Unresolved" in summary["Candidate effects"]
    assert "not calibrated probabilities" in summary["Outcome preference"]
    # Preserve the source products and hypotheses in existing native fragment
    # provenance. There is no reclassification of unknown outcomes as safe.
    path = fragment.dna_effect_selection["candidate_path"]
    evidence = path[-1]["evidence"]
    assert path[-1]["source"] == "varcode.transcript_model"
    assert evidence["probability"] is None
    assert evidence["outcomes"][0]["baseline"]["cdna_sequence"]
    assert evidence["outcomes"][0]["mutant"]["protein_sequence"] == priority.mutant_protein_sequence
    assert evidence["outcomes"][0]["hypothesis"]["splice_plan"]
    restored = MutantProteinFragment.from_json(fragment.to_json())
    assert restored.dna_effect_selection == fragment.dna_effect_selection
    assert restored.predicted_effect().mutant_protein_sequence == priority.mutant_protein_sequence
    from vaxrank.native_serialization import from_native_json, to_native_json
    archived = from_native_json(to_native_json(fragment), MutantProteinFragment)
    assert archived.dna_effect_selection == fragment.dna_effect_selection
    assert archived.predicted_effect().mutant_protein_sequence == priority.mutant_protein_sequence
    # Driver annotation uses Variant.effects(), while report metadata uses
    # effect_on_transcript(). Neither may require the removed keyword.
    assert GenePathwayCheck().make_variant_dict(variant)


def test_template_data_include_manufacturability_override():
    """The split-report layout (#17) toggles manufacturability per
    report via the include_manufacturability override: forced off for
    the modality-agnostic core + mRNA reports, following the run's
    --manufacturability resolution (None) for the peptide report."""
    from types import SimpleNamespace
    from vaxrank.report import TemplateDataCreator

    args = {'manufacturability': True, 'wt_epitopes': True,
            'vaccine_type': ['peptide', 'mrna']}
    patient = SimpleNamespace(patient_id='p', mhc_alleles=[])

    def make(**kw):
        return TemplateDataCreator(
            candidate_details=True,
            ranked_variants_with_vaccine_peptides=[], patient_info=patient,
            final_review='', reviewers='', args_for_report=args,
            input_json_file=None, **kw).template_data

    assert make()['include_manufacturability'] is True          # follows args
    assert make(include_manufacturability=None)[
        'include_manufacturability'] is True                    # explicit follow
    assert make(include_manufacturability=False)[
        'include_manufacturability'] is False                   # forced off (core/mrna)


def test_shared_candidate_paths_preserve_nested_provenance_and_unknown_cycles():
    from varcode import iter_candidate_paths
    from topiary import varcode_candidate_provenance
    outer, likely, priority = _outcome_set()
    inner = _FakeOutcomeSet(likely, priority)
    outer._candidates = (EffectCandidate(inner, source="outer", evidence={"sample": "x"}),)
    inner._candidates = (
        EffectCandidate(likely, source="first"),
        EffectCandidate(priority, source="second"),
        EffectCandidate(outer, source="cycle"),
    )
    paths = list(iter_candidate_paths(outer))
    assert [effect for effect, _path in paths] == [likely, priority, None]
    assert [[c.source for c in path] for _effect, path in paths] == [
        ["outer", "first"], ["outer", "second"], ["outer", "cycle"]]
    assert varcode_candidate_provenance(paths[0][1][0]) == dict(
        source="outer", evidence={"sample": "x"})
