"""Cleavage evidence belongs to a source occurrence, not just a sequence."""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from vaxrank.processing import annotate_processing
from vaxrank.report import TemplateDataCreator

from .test_processing import StubPepsickle, _ep


@pytest.mark.parametrize("reverse", [False, True])
def test_repeated_ligands_keep_distinct_processing_and_report_scores(reverse):
    source = "AAAAAKLMNPVGGGKLMNPVCCCC"
    peptide = "KLMNPV"
    first, second = (_ep(peptide, source, offset) for offset in (5, 14))
    probabilities = [0.1] * len(source)
    probabilities[10] = 0.2
    probabilities[16] = 0.7
    probabilities[19] = 0.9
    epitopes = [second, first] if reverse else [first, second]
    predictor = StubPepsickle({source: probabilities})

    count, records = annotate_processing(epitopes, predictor=predictor)

    assert count == 2
    assert len(records) == 2
    assert predictor.call_count == 1
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = records
    for epitope, expected_c_term, expected_internal in (
            (first, 0.2, 0.1), (second, 0.9, 0.7)):
        key = (peptide, source, epitope.offset, "pepsickle")
        assert records[key].peptide_offset == epitope.offset
        assert records[key].c_term_cleavage_prob == expected_c_term
        assert records[key].max_internal_cut_prob == expected_internal
        row = creator.epitope_data(
            epitope, epitope.best_affinity(), include_processing=True)
        assert row["Processing: C-term"] == f"{expected_c_term:.2f}"
        assert row["Processing: max internal"] == f"{expected_internal:.2f}"


@pytest.mark.parametrize("offset", [0, 9, 8])
def test_resolved_occurrence_is_shared_between_annotation_and_report(offset):
    source = "KLMNPVGGGKLMNPV"
    epitope = _ep("KLMNPV", source, offset)
    probabilities = [0.1] * len(source)
    probabilities[5] = 0.2
    probabilities[14] = 0.9
    count, records = annotate_processing(
        [epitope], predictor=StubPepsickle({source: probabilities}))
    resolved = 0 if offset == 0 else 9
    assert count == 1
    record, = records.values()
    assert record.peptide_offset == resolved
    assert epitope.offset == offset  # Annotation must not mutate the input.
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = records
    row = creator.epitope_data(
        epitope, epitope.best_affinity(), include_processing=True)
    assert row["Processing: C-term"] == ("0.20" if resolved == 0 else "0.90")


def test_one_occurrence_can_share_processing_across_alleles():
    source = "KLMNPVGGG"
    first = _ep("KLMNPV", source, 0)
    second = _ep("KLMNPV", source, 0, allele="HLA-B*07:02")
    count, records = annotate_processing(
        [first, second], predictor=StubPepsickle({source: [0.1] * len(source)}))
    assert count == 2
    assert len(records) == 1


def test_unscored_occurrence_does_not_borrow_another_occurrences_score():
    source = "KLMNPVGGGKLMNPV"
    first = _ep("KLMNPV", source, 0)
    second = _ep("KLMNPV", source, 9)
    _, records = annotate_processing(
        [first], predictor=StubPepsickle({source: [0.1] * len(source)}))
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = records
    row = creator.epitope_data(
        second, second.best_affinity(), include_processing=True)
    assert row["Processing: C-term"] == "—"


@pytest.mark.parametrize("source_class", ["self", "virus"])
def test_processing_occurrence_is_not_mutation_specific(source_class):
    source = "KLMNPVGGG"
    epitope = replace(
        _ep("KLMNPV", source, 0), source_class=source_class, overlaps_mutation=False)
    count, records = annotate_processing(
        [epitope], predictor=StubPepsickle({source: [0.1] * len(source)}))
    assert count == 1
    assert (epitope.sequence, source, 0, "pepsickle") in records


def test_cli_annotation_preserves_occurrences_across_input_paths(monkeypatch):
    from vaxrank.cli.entry_point import annotate_predictions_with_processing

    source = "KLMNPVGGGKLMNPV"
    first = _ep("KLMNPV", source, 0)
    second = _ep("KLMNPV", source, 9)
    predictor = StubPepsickle({source: [0.1] * len(source)})
    monkeypatch.setattr(
        "vaxrank.processing.load_default_processing_predictor",
        lambda **kwargs: predictor)
    ranked = [(None, [SimpleNamespace(target_epitopes=[first, second])])]
    count, records = annotate_predictions_with_processing(ranked, [replace(second)])
    assert count == 2
    assert len(records) == 2
    assert predictor.call_count == 1
