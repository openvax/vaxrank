"""Sequence-end sentinels are not cleavage probabilities (#434)."""

from dataclasses import asdict, replace

import pytest

from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.processing import annotate_processing, processing_component_probabilities
from vaxrank.processing_prediction import ProcessingPrediction
from vaxrank.report import TemplateDataCreator, make_ascii_report, make_html_report

from .test_processing import StubPepsickle, _ep


@pytest.mark.parametrize("sentinel", [0.0, .95])
def test_endpoint_retains_internal_evidence_but_no_cleavage_or_composite(sentinel, tmp_path):
    source, peptide = "GKRFHATISFDTDTGLKQALETKK", "LKQALETKK"
    scores = [.2] * len(source)
    scores[17], scores[-1] = .9177, sentinel
    epitope = _ep(peptide, source, 15)
    n, records = annotate_processing([epitope], predictor=StubPepsickle({source: scores}))
    record, = records.values()
    assert n == 1
    assert record.c_boundary_status == "sequence_endpoint"
    assert record.c_term_cleavage_prob is None and record.processing_score is None
    assert record.max_internal_cut_prob == .9177
    restored = from_native_json(to_native_json(record), ProcessingPrediction)
    assert restored == record and restored.c_boundary_status == "sequence_endpoint"
    assert asdict(restored)["processing_score"] is None
    assert asdict(restored)["c_boundary_status"] == "sequence_endpoint"
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = records
    row = creator.epitope_data(epitope, epitope.best_affinity(), include_processing=True)
    assert row["Processing: C-term"] == "sequence endpoint"
    assert row["Processing: max internal"] == "0.92"
    assert row["Processing: combined"] == "—"
    # Exercise the actual report templates, not just the Python row mapping.
    template_data = dict(patient_info={}, package_versions={}, args=[], vaccine_constructions={},
        variants=[dict(num=1, short_description="endpoint regression", variant_data={}, effect_data={},
            peptides=[dict(header_display_data={}, peptide_data={}, epitopes=[row],
                           ascii_epitopes=" | ".join(map(str, row.values())))])])
    for write, suffix in ((make_ascii_report, "txt"), (make_html_report, "html")):
        path = tmp_path / ("endpoint." + suffix)
        write(template_data, path)
        assert "sequence endpoint" in path.read_text()


def test_actual_internal_zero_remains_zero_and_single_residue_endpoint_is_missing():
    assert processing_component_probabilities([.1, 0., .9, 0.], 0, 2) == (0., .1)
    assert processing_component_probabilities([0.], 0, 1) == (None, 0.)


@pytest.mark.parametrize("scores", [[.1] * 8, [.1] * 10, [float("nan")] * 9, [1.1] * 9])
def test_invalid_length_or_probabilities_cannot_turn_endpoint_into_a_real_bond(scores, caplog):
    source = "KLMNPVGGG"
    assert annotate_processing([_ep("KLMNPV", source, 0)], predictor=StubPepsickle({source: scores})) == (0, {})
    assert "Invalid cleavage output" in caplog.text


@pytest.mark.parametrize("change", [{"peptide_offset": -1}, {"peptide_sequence": "AAA"},
    {"c_term_cleavage_prob": 0.}, {"processing_score": 0.}, {"c_boundary_status": "observed"},
    {"max_internal_cut_prob": float("nan")}, {"max_internal_cut_prob": 2.}])
def test_native_endpoint_record_cannot_silently_assert_inconsistent_evidence(change):
    record = ProcessingPrediction("KLM", "GKLM", "pepsickle", 1, max_internal_cut_prob=.2)
    with pytest.raises(ValueError):
        replace(record, **change)
