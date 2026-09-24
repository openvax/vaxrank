"""Configured construct scores must reach selection, assembly and reports."""

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from varcode import Variant

from vaxrank.candidate_epitope import CandidateEpitope, Peptide
from vaxrank.core_logic import ranked_vaccine_peptides, vaccine_peptides_from_epitopes
from vaxrank.external_input import load_external_ranked
from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.patient_info import PatientInfo
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.report import TemplateDataCreator, make_ascii_report
from vaxrank.vaccine_config import VaccineConfig
from vaxrank.vaccine_library import iter_named_antigens

DATA = Path(__file__).parent / "data" / "epitope_fixtures"
LENS = "lens=" + str(DATA / "lens_example.tsv")
PVACSEQ = "pvacseq=" + str(DATA / "pvacseq_example.tsv")
INVERSE_SCORE = "1 / (1 + target_epitope_score)"


def _load(inputs, expression, **config):
    return load_external_ranked(
        SimpleNamespace(external_input=inputs),
        vaccine_config=VaccineConfig(combined_score_expr=expression, **config),
    )[0]


def _assert_outputs_choose_lowest_epitope_score(ranked, tmp_path):
    """Exercise both capped assemblers and the data used by rendered reports."""
    scores = [peptides[0].target_epitope_score for _, peptides in ranked]
    assert len(scores) >= 2 and scores[0] < scores[-1]
    assert scores == sorted(scores)
    expected_name = next(iter_named_antigens(ranked))[0]
    peptide_constructs = assemble_peptide_constructs(
        ranked, options=PeptideConstructConfig(
            min_antigen_length_aa=5, max_constructs=1))
    assert len(peptide_constructs) == 1
    assert peptide_constructs[0].antigen_names == [expected_name]
    assert peptide_constructs[0].sequence in ranked[0][1][0].amino_acids
    mrna_constructs = assemble_mrna_constructs(
        ranked, options=RNAConstructConfig(
            signal_peptide=None, include_mitd=False,
            min_antigen_length_aa=5, antigens_per_construct=1,
            max_constructs=1, optimize_linkers=False))
    assert len(mrna_constructs) == 1
    assert mrna_constructs[0].antigen_names == [expected_name]
    data = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked,
        patient_info=PatientInfo(patient_id="construct-ranking"),
        final_review=None, reviewers=None,
        args_for_report={"manufacturability": False, "wt_epitopes": False},
        input_json_file=None,
    ).compute_template_data()
    reported = [float(row["variant_data"]["Top score"]) for row in data["variants"]]
    assert reported == pytest.approx([
        peptides[0].combined_score for _, peptides in ranked], abs=1e-5)
    assert reported == sorted(reported, reverse=True)
    report_path = tmp_path / "ranking.txt"
    make_ascii_report(data, str(report_path))
    report = report_path.read_text()
    positions = [report.index(row["short_description"]) for row in data["variants"]]
    assert positions == sorted(positions)


@pytest.mark.parametrize("inputs", [[LENS], [PVACSEQ], [LENS, PVACSEQ]])
def test_external_combined_score_controls_final_constructs_and_reports(inputs, tmp_path):
    baseline = _load(inputs, "target_epitope_score")
    ranked = _load(inputs, INVERSE_SCORE)
    assert [str(source) for source, _ in ranked] != [str(source) for source, _ in baseline]
    assert {str(source) for source, _ in ranked} == {str(source) for source, _ in baseline}
    _assert_outputs_choose_lowest_epitope_score(ranked, tmp_path)


def test_direct_selected_windows_follow_the_same_combined_score_policy(tmp_path):
    selected = {}
    for index, (sequence, score) in enumerate([
        ("AAAASIINFEKLAAA", 0.9), ("AAAAYLLPAIVHIAA", 0.2),
    ]):
        variant = Variant("1", 100 + index, "A", "T", genome=None)
        fragment = MutantProteinFragment(
            variant=variant, gene_name="GENE%d" % index, amino_acids=sequence,
            mutant_amino_acid_start_offset=6, mutant_amino_acid_end_offset=7,
            supporting_reference_transcripts=[], n_overlapping_reads=20,
            n_alt_reads=10, n_ref_reads=10,
            n_alt_reads_supporting_protein_sequence=10,
            n_alt_fragments=5, rna_evidence_subject="fragments",
            rna_evidence_method="rna_alignment",
        )
        epitope = CandidateEpitope.from_peptide(
            Peptide(sequence=sequence[4:13], source_sequence=sequence, offset=4),
            overlaps_mutation=True, occurs_in_reference=False,
            per_allele_scores={"HLA-A*02:01": score},
        )
        selected[variant] = vaccine_peptides_from_epitopes(
            variant, fragment, [epitope], vaccine_config=VaccineConfig(
                preferred_peptide_length=15, min_peptide_length=15,
                max_peptide_length=15, combined_score_expr=INVERSE_SCORE))
    assert all(selected.values())
    ranked = ranked_vaccine_peptides(selected)
    assert ranked[0][0].start == 101
    _assert_outputs_choose_lowest_epitope_score(ranked, tmp_path)


@pytest.mark.parametrize("reverse_inputs", [False, True])
def test_repeated_sources_choose_by_combined_score_without_combining_evidence(
        tmp_path, reverse_inputs):
    frame = pd.read_csv(DATA / "lens_example.tsv", sep="\t")
    for metric in ("netmhcpan_4.1b.aff_nm", "mhcflurry_2.1.1.aff"):
        frame[metric] = 450
    second = tmp_path / "second.tsv"
    frame.to_csv(second, sep="\t", index=False)
    second_input = "lens=" + str(second)
    original = _load([LENS], INVERSE_SCORE)
    lower_binding = _load([second_input], INVERSE_SCORE)
    assert max(p[0].target_epitope_score for _, p in lower_binding) < min(
        p[0].target_epitope_score for _, p in original)
    inputs = [LENS, second_input]
    if reverse_inputs:
        inputs.reverse()
    ranked = _load(inputs, INVERSE_SCORE)
    assert len(ranked) == len(original)
    expected = {str(s): p[0] for s, p in lower_binding}
    for source, peptides in ranked:
        actual = peptides[0]
        assert actual.combined_score == expected[str(source)].combined_score
        assert [e.prediction_id for e in actual.epitopes] == [
            e.prediction_id for e in expected[str(source)].epitopes]


def test_source_agnostic_antigens_rank_without_mutation_rna(tmp_path):
    path = DATA / "real_lens_subsets" / "lens_v1.9_real_subset.tsv"
    ranked = _load(["lens=" + str(path)], INVERSE_SCORE,
                   included_antigen_sources=("FUSION", "SPLICE", "CTA/SELF", "ERV"))
    assert {p[0].antigen.kind for _, p in ranked} == {"fusion", "splice", "CTA", "ERV"}
    assert all(p[0].mutant_protein_fragment is None for _, p in ranked)
    _assert_outputs_choose_lowest_epitope_score(ranked, tmp_path)


def test_default_source_agnostic_scores_and_missing_rna_remain_explicit():
    path = DATA / "real_lens_subsets" / "lens_v1.9_real_subset.tsv"
    ranked, *_ = load_external_ranked(
        SimpleNamespace(external_input=[LENS, "lens=" + str(path)]))
    antigens = [p[0] for _, p in ranked if p[0].mutant_protein_fragment is None]
    assert antigens
    assert all(p.combined_score_expr == "target_epitope_score" for p in antigens)
    assert all(p.combined_score == p.target_epitope_score for p in antigens)
    # The simple LENS fixture supplies no RNA counts or derivation. Ordering
    # must not label its legacy zero placeholders as measured support.
    missing = _load([LENS], INVERSE_SCORE)
    assert all(not p[0].mutant_protein_fragment.rna_evidence_method for _, p in missing)
    assert all(p[0].mutant_protein_fragment.n_alt_fragments is None for _, p in missing)


def test_occurrence_selection_stays_separate_from_final_construct_ranking(tmp_path):
    path = tmp_path / "contexts.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tAAAASIINFEKLAAA\tSNV\tchr1:100\tA\tT\tG\t450\n"
        "YLLPAIVHI\tHLA-A02:01\tAAAAYLLPAIVHIAA\tSNV\tchr1:100\tA\tT\tG\t10\n")
    # Within one report the best eligible epitope still picks its own context.
    # The final inverse score must not change which occurrence was selected.
    ranked = _load(["lens=" + str(path)], INVERSE_SCORE)
    assert len(ranked) == 1
    peptide = ranked[0][1][0]
    assert peptide.amino_acids == "AAAAYLLPAIVHIAA"
    assert [e.sequence for e in peptide.epitopes] == ["YLLPAIVHI"]


@pytest.mark.parametrize("source_format", ["lens", "pvacseq"])
def test_external_target_admission_preserves_audit_and_honors_opt_out(source_format):
    from vaxrank.epitope_io import read_lens_report, read_pvacseq_report
    from vaxrank.external_input import (
        ExternalConstructOptions, lens_ranking_result, pvacseq_ranking_result,
    )

    reader, ranker = {
        "lens": (read_lens_report, lens_ranking_result),
        "pvacseq": (read_pvacseq_report, pvacseq_ranking_result),
    }[source_format]
    report = reader(str(DATA / (source_format + "_example.tsv")))
    self_only = [replace(e, occurs_in_reference=True) for e in report.epitopes]
    options = ExternalConstructOptions.from_configs(VaccineConfig(
        require_target_epitopes_in_variant=False, combined_score_expr=INVERSE_SCORE))
    admitted = ranker(report, self_only, options=options)
    filtered = ranker(report, self_only)
    assert len(admitted.ranked) == 3
    assert all(not p[0].contains_target_epitopes() for _, p in admitted.ranked)
    assert not filtered.ranked
    assert len(filtered.entries) == len(admitted.entries) == 3
    assert all(entry.vaccine_peptide is None for entry in filtered.entries)
    assert filtered.input_summary == admitted.input_summary


def test_unified_admission_excludes_self_only_without_losing_input_evidence(tmp_path):
    # A self-only observation must not replace the admitted construct from
    # another report just because its inverse combined score is larger.
    frame = pd.read_csv(DATA / "pvacseq_example.tsv", sep="\t")
    frame["Ref Match"] = True
    path = tmp_path / "self-only.tsv"
    frame.to_csv(path, sep="\t", index=False)
    inputs = ["pvacseq=" + str(path), PVACSEQ]
    ranked, report, predictions, patient, _ = load_external_ranked(
        SimpleNamespace(external_input=inputs),
        vaccine_config=VaccineConfig(combined_score_expr=INVERSE_SCORE))
    assert len(ranked) == 2
    assert {p[0].antigen.gene_name for _, p in ranked} == {"TP53", "BRAF"}
    assert all(p[0].contains_target_epitopes() for _, p in ranked)
    assert report["Input source"].nunique() == 2
    assert len(predictions) == 6
    assert patient.num_somatic_variants == 3
    assert patient.num_variants_with_vaccine_peptides == 2
    opt_out = _load(inputs, INVERSE_SCORE, require_target_epitopes_in_variant=False)
    assert len(opt_out) == 3
    assert all(not p[0].contains_target_epitopes() for _, p in opt_out)
