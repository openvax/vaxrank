"""Malformed evidence fails at input, not while rendering other results."""

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from mhctools import Prediction
import pandas as pd
import pytest
from topiary import KIND_MHC_DEPENDENCE

from vaxrank.allele_validation import validate_prediction_allele
from vaxrank.candidate_epitope import CandidateEpitope, Peptide
from vaxrank.epitope_io import (
    NATIVE_EPITOPE_COLUMN, load_predictions, predictions_to_dataframe,
    read_lens_report, read_pvacseq_report, save_predictions,
)
from vaxrank.report import JINJA_ENVIRONMENT, TemplateDataCreator, epitope_report_row_inputs


def candidate(allele="HLA-A*02:01", kind="pMHC_affinity"):
    return CandidateEpitope(
        sequence="SIINFEKL", source_sequence="SIINFEKL", offset=0,
        predictions=(Prediction(kind=kind, predictor_name="example", predictor_version="1",
                                peptide="SIINFEKL", allele=allele, score=.3, value=100.),))


@pytest.mark.parametrize("kind", [k for k, scope in KIND_MHC_DEPENDENCE.items() if scope == "single_allele"])
@pytest.mark.parametrize("allele", [None, "", " ", "nan", "NA", "."])
def test_all_allele_scoped_kinds_require_a_stated_allele(kind, allele):
    with pytest.raises(ValueError, match="input.tsv, row 7.*requires a patient allele"):
        validate_prediction_allele(kind, allele, "input.tsv, row 7")


@pytest.mark.parametrize("kind", [k for k, scope in KIND_MHC_DEPENDENCE.items() if scope == "none"])
def test_allele_free_kinds_remain_valid_without_a_patient_allele(kind):
    validate_prediction_allele(kind, "", "input.tsv, row 7")


@pytest.mark.parametrize("allele", ["HLA-A*02:01", "H-2-Kb", "Mamu-A*001:01",
                                  "HLA-DPA1*01:03-DPB1*04:01", "HLA-DRB1*08:01"])
def test_real_alleles_and_class_ii_pairs_are_accepted(allele):
    validate_prediction_allele("pMHC_affinity", allele, "input.tsv, row 7")


def test_garbage_allele_has_actionable_error():
    with pytest.raises(ValueError, match="input.tsv, row 7.*invalid allele.*not-an-allele"):
        validate_prediction_allele("pMHC_affinity", "not-an-allele", "input.tsv, row 7")


@pytest.mark.parametrize("native_payload", [True, False])
def test_native_file_rejects_missing_allele_with_filename_row_and_kind(tmp_path, native_payload):
    path = tmp_path / "broken.csv"
    frame = predictions_to_dataframe([candidate(allele="")])
    if not native_payload:
        frame = frame.drop(columns=[NATIVE_EPITOPE_COLUMN])
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="broken.csv.*native row 2.*pMHC_affinity"):
        load_predictions(path)


def test_non_wt_comparator_cannot_hide_a_missing_allele(tmp_path):
    original = replace(candidate(), comparators={
        "nearest_non_CTA": Peptide(sequence="SIINFEKL", predictions=candidate(allele="").predictions)})
    path = tmp_path / "comparator.tsv"
    save_predictions([original], path)
    with pytest.raises(ValueError, match="comparator.tsv.*row 2.*nearest_non_CTA.*requires"):
        load_predictions(path)
    with pytest.raises(ValueError, match="Native peptide JSON.*nearest_non_CTA.*requires"):
        CandidateEpitope.from_json(original.to_json())


def test_canonical_allele_free_native_file_loads(tmp_path):
    original = candidate(allele="", kind="antigen_processing")
    path = tmp_path / "processing.csv"
    save_predictions([original], path)
    assert load_predictions(path) == [original]


@pytest.mark.parametrize("reader,fixture,column", [
    (read_lens_report, "lens_example.tsv", "allele"),
    (read_pvacseq_report, "pvacseq_example.tsv", "Allele"),
])
def test_external_loader_rejects_missing_allele_before_scoring(tmp_path, reader, fixture, column):
    source = Path(__file__).parent / "data" / "epitope_fixtures" / fixture
    frame = pd.read_csv(source, sep="\t")
    assert column in frame
    frame.loc[0, column] = ""
    path = tmp_path / fixture
    frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match=fixture + ".*row.*requires a patient allele"):
        reader(path)


def test_report_json_reload_validates_before_output_overrides(tmp_path):
    from vaxrank.cli.entry_point import ranked_vaccine_peptides_with_metadata_from_parsed_args

    path = tmp_path / "prior-report.json"
    path.write_text("{}")
    data = {"variants": [(None, [SimpleNamespace(epitopes=[candidate(allele="")])])], "args": {}}
    args = SimpleNamespace(input_json_file=str(path), max_mutations_in_report=None)
    with patch("vaxrank.cli.entry_point.serializable.from_json", return_value=data):
        with pytest.raises(ValueError, match="prior-report.json.*variant 1.*vaccine peptide 1.*requires"):
            ranked_vaccine_peptides_with_metadata_from_parsed_args(args)


def test_lens_processing_without_typing_is_not_malformed_mhc_evidence(tmp_path):
    path = tmp_path / "untyped-processing.tsv"
    pd.DataFrame([{"allele": "", "peptide": "SIINFEKL", "pep_context": "SSIINFEKL",
                   "mhcflurry_2.1.1.proc_score": .8}]).to_csv(path, sep="\t", index=False)
    report = read_lens_report(path)
    assert len(report.epitopes) == 1
    epitope = report.epitopes[0]
    assert epitope.patient_alleles == ()
    assert all(p.allele == "" for p in epitope.predictions_flat())
    rows = epitope_report_row_inputs(epitope)
    assert len(rows) == 1 and rows[0].allele == ""


@pytest.mark.parametrize("kind", ["pMHC_affinity", "pMHC_presentation", "antigen_processing"])
def test_direct_legacy_report_objects_render_unknown_allele_without_crashing(kind):
    epitope = candidate(allele="", kind=kind)
    inputs = epitope_report_row_inputs(epitope)
    assert len(inputs) == 1 and inputs[0].allele == ""
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    row = creator.epitope_data(epitope, inputs[0].prediction, allele=inputs[0].allele)
    assert row["Allele"] == "" and row["Score"] is None
    scored = replace(candidate(), per_allele_scores={"HLA-A*02:01": .25})
    scored_row = creator.epitope_data(scored, scored.predictions_flat()[0])
    assert pd.api.types.is_numeric_dtype(pd.DataFrame([row, scored_row])["Score"])
    template = JINJA_ENVIRONMENT.from_string("|{{ score|display_epitope_value }}|")
    assert template.render(score=row["Score"]) == "||"
    assert template.render(score=scored_row["Score"]) == "|0.25|"


@pytest.mark.parametrize("bad_output", ["source", "WT"])
def test_predictor_output_validation_is_not_swallowed_as_a_backend_failure(bad_output):
    from topiary import TopiaryPredictor
    from varcode import Variant
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_logic import predict_epitopes
    from vaxrank.mutant_protein_fragment import MutantProteinFragment

    class Predictor(TopiaryPredictor):
        def __init__(self):
            pass

        def frame(self, named, output):
            name, sequence = next(iter(named.items()))
            return pd.DataFrame([{
                "source_sequence_name": name, "peptide": sequence, "peptide_offset": 0,
                "peptide_length": len(sequence), "kind": "pMHC_affinity",
                "allele": "" if output == bad_output else "HLA-A*02:01",
                "prediction_method_name": "example", "predictor_version": "1",
                "affinity": 10., "value": 10., "score": .8, "percentile_rank": .1}])

        def predict_from_named_sequences(self, named):
            return self.frame(named, "source")

        def predict_from_named_peptides(self, named):
            return self.frame(named, "WT")

    fragment = MutantProteinFragment(
        variant=Variant("1", 100, "A", "T"), gene_name="GENE",
        amino_acids="SIINFEKL", mutant_amino_acid_start_offset=3,
        mutant_amino_acid_end_offset=4, supporting_reference_transcripts=[],
        n_overlapping_reads=10, n_alt_reads=5, n_ref_reads=5,
        n_alt_reads_supporting_protein_sequence=5)
    with patch("vaxrank.epitope_logic.ReferenceProteome"), patch.object(
            MutantProteinFragment, "predicted_effect",
            return_value=SimpleNamespace(original_protein_sequence="SIISFEKL")), patch.object(
                MutantProteinFragment, "global_start_pos", return_value=0):
        prefix = "WT MHC" if bad_output == "WT" else "MHC"
        with pytest.raises(ValueError, match=prefix + " output.*normalized row 2.*requires"):
            predict_epitopes(Predictor(), fragment, EpitopeConfig(min_epitope_score=0))


@pytest.mark.parametrize("score", [None, float("nan"), float("inf"), -float("inf")])
def test_nonfinite_report_score_is_missing_not_zero(score):
    epitope = replace(candidate(), per_allele_scores={"HLA-A*02:01": score})
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    row = creator.epitope_data(epitope, epitope.predictions_flat()[0])
    assert row["Score"] is None


def test_template_pipeline_keeps_unknown_and_known_allele_rows(tmp_path):
    """Exercise real table assembly and packaged HTML/ASCII templates."""
    from vaxrank.patient_info import PatientInfo
    from vaxrank.report import make_ascii_report, make_html_report

    unknown = replace(candidate(allele=""), overlaps_mutation=True)
    known = replace(candidate(), overlaps_mutation=True,
                    per_allele_scores={"HLA-A*02:01": .25})
    fragment = SimpleNamespace(gene_name="GENE", predicted_effect=lambda **kw: None)
    vaccine = SimpleNamespace(
        mutant_protein_fragment=fragment, target_epitopes=[unknown, known],
        self_epitopes=[], contains_target_epitopes=lambda: True)
    creator = TemplateDataCreator(
        [(SimpleNamespace(short_description="test variant"), [vaccine])],
        PatientInfo("TEST"), final_review="", reviewers="", input_json_file=None,
        args_for_report={"manufacturability": False, "wt_epitopes": False})
    # Variant annotation is orthogonal to epitope table rendering.
    for method in ("_variant_data", "effect_data", "_databases", "_peptide_data",
                   "_manufacturability_data"):
        setattr(creator, method, lambda *a, **kw: {})
    creator._peptide_header_display_data = lambda *a: {
        "num": 1, "aa_before_mutation": "SII", "aa_mutant": "N", "aa_after_mutation": "FEKL"}
    data = creator.compute_template_data()
    peptide = data["variants"][0]["peptides"][0]
    assert len(peptide["epitopes"]) == 2
    assert peptide["epitopes"][0]["Score"] is None
    assert peptide["epitopes"][1]["Score"] == .25
    assert "None" not in peptide["ascii_epitopes"]
    for suffix, writer in (("html", make_html_report), ("txt", make_ascii_report)):
        path = tmp_path / ("report." + suffix)
        writer(data, path)
        rendered = path.read_text()
        assert rendered.count("SIINFEKL") >= 2
        assert "0.25" in rendered and "A*02:01" in rendered
