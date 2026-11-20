"""Tests for longitudinal library support matrix (#818) and matched-normal evidence (#819)."""

from dataclasses import asdict
import json

from mhctools import Prediction
from pyensembl import EnsemblRelease
from varcode import Variant

from vaxrank.antigen_input import AntigenDataset
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.cta_admission import CTAAdmissionPolicy
from vaxrank.cta_expression import (
    CTAExpressionDecision,
    CTAExpressionInput,
    CTAExpressionResult,
    CTASequenceAssessment,
    CTATargetExclusionPolicy,
)
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.input_scope import InputProvenance, InputScope
from vaxrank.isovar_antigens import antigens_from_isovar
from vaxrank.isovar_input import load_isovar_input
from vaxrank.longitudinal_support import build_longitudinal_support_view
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.native_serialization import to_native_json
from vaxrank.normal_observations import (
    MatchedNormalObservation,
    group_normal_observations_by_library,
    summarize_normal_assessment,
)
from vaxrank.patient_info import PatientInfo
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.report import TemplateDataCreator, make_ascii_report, make_html_report
from vaxrank.vaccine_antigen import (
    AminoAcidInterval,
    TargetableMask,
    TumorSpecificityAttestation,
    VaccineAntigen,
)
from vaxrank.vaccine_peptide import VaccinePeptide


def _pred(peptide):
    return Prediction(
        peptide=peptide,
        allele="HLA-A*02:01",
        kind="pMHC_affinity",
        value=18.0,
        score=0.92,
        percentile_rank=0.15,
        predictor_name="mock",
        predictor_version="1.0",
    )


def _ep(seq, source_seq, pid, prov=None, offset=0):
    return CandidateEpitope(
        sequence=seq,
        source_sequence=source_seq,
        offset=offset,
        prediction_id=pid,
        occurs_in_reference=False,
        overlaps_targetable=True,
        predictions=[_pred(seq)],
        per_allele_scores={"HLA-A*02:01": 1.0},
        input_provenance=prov,
    )


def _cta_antigen(gene_id, gene_name, transcript_id, seq, admitted=True):
    return VaccineAntigen(
        kind="CTA",
        amino_acids=seq,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(seq)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status="admitted" if admitted else "held_out",
            evidence_kind="cta_tumor_expression",
            evidence_source="rsem",
            patient_specific=True,
            rationale_code="" if admitted else "below_expression_threshold",
        ),
        gene_name=gene_name,
        gene_id=gene_id,
        transcript_ids=(transcript_id,),
        protein_ids=("ENSP_" + gene_name,),
        species="homo_sapiens",
        source_identifier=gene_id,
        source_protein_completeness="complete",
    )


def _cta_result(sample_id, sha_prefix, decisions):
    return CTAExpressionResult(
        input_contract=CTAExpressionInput(
            measurement_level="gene",
            id_column="gene_id",
            value_column="TPM",
            expression_unit="TPM",
            sample_id=sample_id,
            evidence_source="RSEM",
            evidence_version="1.3.3",
            assay="RNA-seq",
        ),
        admission_policy=CTAAdmissionPolicy(2.0, "TPM"),
        exclusion_policy=CTATargetExclusionPolicy(gene_patterns=(), gene_exceptions=()),
        input_sha256=(sha_prefix * 64)[:64],
        reference_assembly="GRCh38",
        annotation_name="ensembl",
        annotation_version="112",
        decisions=tuple(decisions),
    )


def test_vanilla_single_bam_vcf_omits_longitudinal_matrix():
    genome = EnsemblRelease(75)
    variant = Variant("12", 25398284, "C", "A", genome=genome)
    frag = MutantProteinFragment(
        variant=variant,
        gene_name="KRAS",
        amino_acids="MTEYKLVVVGADGVGKSALTI",
        mutant_amino_acid_start_offset=11,
        mutant_amino_acid_end_offset=12,
        supporting_reference_transcripts=[],
        n_overlapping_reads=40,
        n_alt_reads=15,
        n_ref_reads=25,
        n_alt_reads_supporting_protein_sequence=15,
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=frag,
        epitopes=[_ep("GADGVGKSA", frag.amino_acids, "kras-1")],
    )
    pinfo = PatientInfo(
        patient_id="OpenVax1-Patient",
        vcf_paths=["somatic.vcf"],
        bam_path="tumor.bam",
        mhc_alleles=["HLA-A*02:01"],
        num_somatic_variants=1,
        num_coding_effect_variants=1,
        num_variants_with_rna_support=1,
        num_variants_with_vaccine_peptides=1,
    )
    assert build_longitudinal_support_view([(variant, [vp])], {}, {}, pinfo) is None


def test_longitudinal_support_matrix_cta_below_threshold_zero_missing_and_replay(tmp_path):
    cta_seq = "SIINFEKLAKK"
    mut_seq_t0 = "GILGFVFTLKK"
    mut_seq_t2_other = "KLGFVFTLKKQ"

    prov_t0_cta = InputProvenance(
        input_id="in-t0-cta",
        source_format="cta_expression",
        path="t0_cta.tsv",
        content_sha256="a" * 64,
        scope=InputScope(
            patient_id="P1",
            reference_assembly="GRCh38",
            mhc_alleles=("HLA-A*02:01",),
            sample_id="T0-tumor",
            timepoint="T0",
            library_id="lib-T0",
        ),
    )
    # Dependent second pipeline on the EXACT same (sample_id='T0-tumor', library_id='lib-T0')
    prov_t0_mut = InputProvenance(
        input_id="in-t0-mut",
        source_format="isovar",
        path="t0_mut.json",
        content_sha256="b" * 64,
        scope=InputScope(
            patient_id="P1",
            reference_assembly="GRCh38",
            mhc_alleles=("HLA-A*02:01",),
            sample_id="T0-tumor",
            timepoint="T0",
            library_id="lib-T0",
        ),
    )
    prov_t1_targeted = InputProvenance(
        input_id="in-t1-targeted",
        source_format="isovar",
        path="t1_targeted.json",
        content_sha256="c" * 64,
        scope=InputScope(
            patient_id="P1",
            reference_assembly="GRCh38",
            mhc_alleles=("HLA-A*02:01",),
            sample_id="T1-targeted",
            timepoint="T1",
            library_id="lib-T1",
        ),
    )
    prov_t2_cta = InputProvenance(
        input_id="in-t2-cta",
        source_format="cta_expression",
        path="t2_cta.tsv",
        content_sha256="d" * 64,
        scope=InputScope(
            patient_id="P1",
            reference_assembly="GRCh38",
            mhc_alleles=("HLA-A*02:01",),
            sample_id="T2-tumor",
            timepoint="T2",
            library_id="lib-T2",
        ),
    )

    ag_ctag1b = _cta_antigen("ENSG00000184033", "CTAG1B", "ENST00000328435", cta_seq, admitted=True)
    cta_t0 = _cta_result(
        "T0-tumor",
        "1",
        [
            CTAExpressionDecision(
                input_identifier="ENSG00000184033",
                value=18.5,
                gene_id="ENSG00000184033",
                gene_name="CTAG1B",
                transcript_id="ENST00000328435",
                status="admitted",
                assessment=CTASequenceAssessment(
                    antigen=ag_ctag1b,
                    bundle_id="bundle-1",
                    record_id="rec-1",
                    occurrence={"gene_id": "ENSG00000184033"},
                ),
            ),
            CTAExpressionDecision(
                input_identifier="ENSG00000221867",
                value=0.4,
                gene_id="ENSG00000221867",
                gene_name="MAGEA3",
                transcript_id="ENST00000370239",
                status="held_out",
            ),
            CTAExpressionDecision(
                input_identifier="ENSG00000111640",
                value=512.0,
                gene_id="ENSG00000111640",
                gene_name="GAPDH",
                transcript_id="ENST00000229239",
                status="not_cta",
            ),
        ],
    )
    # At T2: CTAG1B drops to measured 0.0 TPM (< 2.0), MAGEA3 is None (missing!)
    cta_t2 = _cta_result(
        "T2-tumor",
        "2",
        [
            CTAExpressionDecision(
                input_identifier="ENSG00000184033",
                value=0.0,
                gene_id="ENSG00000184033",
                gene_name="CTAG1B",
                transcript_id="ENST00000328435",
                status="held_out",
            ),
            CTAExpressionDecision(
                input_identifier="ENSG00000221867",
                value=None,
                gene_id="ENSG00000221867",
                gene_name="MAGEA3",
                transcript_id="ENST00000370239",
                status="unknown_expression",
            ),
        ],
    )

    genome = EnsemblRelease(75)
    variant = Variant("12", 25398284, "C", "A", genome=genome)
    frag_t0 = MutantProteinFragment(
        variant=variant,
        gene_name="KRAS",
        amino_acids=mut_seq_t0,
        mutant_amino_acid_start_offset=0,
        mutant_amino_acid_end_offset=1,
        supporting_reference_transcripts=[],
        n_overlapping_reads=30,
        n_alt_reads=12,
        n_ref_reads=18,
        n_alt_reads_supporting_protein_sequence=12,
        n_alt_fragments=10,
    )
    frag_t2_alt = MutantProteinFragment(
        variant=variant,
        gene_name="KRAS",
        amino_acids=mut_seq_t2_other,
        mutant_amino_acid_start_offset=0,
        mutant_amino_acid_end_offset=1,
        supporting_reference_transcripts=[],
        n_overlapping_reads=14,
        n_alt_reads=5,
        n_ref_reads=9,
        n_alt_reads_supporting_protein_sequence=5,
        n_alt_fragments=4,
    )

    ep_cta = _ep(cta_seq[:9], cta_seq, "pred-cta", prov=prov_t0_cta)
    ep_mut_t0 = _ep(mut_seq_t0[:9], mut_seq_t0, "pred-mut-t0", prov=prov_t0_mut)
    ep_mut_t2 = _ep(mut_seq_t2_other[:9], mut_seq_t2_other, "pred-mut-t2", prov=prov_t2_cta)

    vp_cta = VaccinePeptide(antigen=ag_ctag1b, epitopes=[ep_cta])
    vp_mut_t0 = VaccinePeptide(mutant_protein_fragment=frag_t0, epitopes=[ep_mut_t0])
    vp_mut_t2 = VaccinePeptide(mutant_protein_fragment=frag_t2_alt, epitopes=[ep_mut_t2])

    ranked = [(ag_ctag1b, [vp_cta]), (variant, [vp_mut_t0, vp_mut_t2])]
    products = {
        "peptide": assemble_peptide_constructs(
            ranked,
            PeptideConstructConfig(min_antigen_length_aa=11, max_antigen_length_aa=11),
            target_alleles=["HLA-A*02:01"],
        )
    }

    selection = {
        "source_evidence": {
            "in-t0-cta": {"cta_expression_admission": to_native_json(cta_t0)},
            "in-t2-cta": {"cta_expression_admission": to_native_json(cta_t2)},
            "iso-t1": {
                "family": "variant_orfs",
                "decisions": [],
                "isovar_result": {
                    "sample_id": "T1-targeted",
                    "variant": {
                        "reference_name": variant.genome.reference_name,
                        "contig": str(variant.contig),
                        "start": int(variant.start),
                        "ref": variant.ref,
                        "alt": variant.alt,
                    },
                    "protein_candidates": [{"rna_support": {"fragments": 2}}],
                }
            },
        }
    }

    # Verify persistence through EpitopeDataset save -> load (replay without sidecars!)
    ds = EpitopeDataset.from_predictions(
        [ep_cta, ep_mut_t0, ep_mut_t2],
        antigens={"pred-cta": ag_ctag1b, "pred-mut-t0": vp_mut_t0.antigen, "pred-mut-t2": vp_mut_t2.antigen},
        mutation_fragments={"pred-mut-t0": frag_t0, "pred-mut-t2": frag_t2_alt},
        config=EpitopeConfig(score_expr="1.0"),
    )
    ds.provenance = (prov_t0_cta, prov_t0_mut, prov_t1_targeted, prov_t2_cta)
    ds.selection.update(selection)
    saved_tsv = tmp_path / "longitudinal.tsv"
    ds.save(saved_tsv)
    reloaded = EpitopeDataset.load(saved_tsv)

    pinfo = PatientInfo(
        patient_id="P1",
        vcf_paths=[],
        bam_path=None,
        mhc_alleles=["HLA-A*02:01"],
        num_somatic_variants=2,
        num_coding_effect_variants=2,
        num_variants_with_rna_support=2,
        num_variants_with_vaccine_peptides=2,
        input_coverage={
            "sources": [],
            "original_inputs": [asdict(p) for p in reloaded.provenance],
            "current_inputs": [],
        },
    )

    view = build_longitudinal_support_view(ranked, products, reloaded.selection, pinfo)
    assert view is not None
    assert view["design_scope"] == "longitudinal_union"
    assert view["latest_timepoint"] == "T2"
    # Dependent pipelines on ("T0-tumor", "lib-T0") share ONE library column (3 total: T0, T1, T2)
    assert [lib["sample_id"] for lib in view["libraries"]] == ["T0-tumor", "T1-targeted", "T2-tumor"]
    assert set(view["libraries"][0]["pipelines"]) >= {"cta_expression", "isovar"}

    by_gene_win = {(r["gene"], r["window"]): r for r in view["rows"]}
    assert "GAPDH" not in {r["gene"] for r in view["rows"]}
    ctag_row = by_gene_win[("CTAG1B", cta_seq)]
    assert ctag_row["selected"] is True
    assert ctag_row["qualifying_timepoints"] == ["T0"]
    assert ctag_row["qualified_only_earlier"] is True
    assert ctag_row["below_threshold_at_latest"] is True
    assert ctag_row["cells"][0]["state"] == "qualified"
    assert "18.5 TPM" in ctag_row["cells"][0]["summary"]
    assert ctag_row["cells"][1]["state"] == "unmeasured"
    assert ctag_row["cells"][2]["state"] == "below_threshold"
    assert "[measured zero]" in ctag_row["cells"][2]["summary"]

    mage_row = by_gene_win[("MAGEA3", "")]
    assert mage_row["selected"] is False
    assert mage_row["cells"][0]["state"] == "below_threshold"
    assert "0.4 TPM" in mage_row["cells"][0]["summary"]
    assert mage_row["cells"][2]["state"] == "missing"
    assert "[measured zero]" not in mage_row["cells"][2]["summary"]

    kras_t0_row = by_gene_win[("KRAS", mut_seq_t0)]
    assert kras_t0_row["cells"][0]["state"] == "exact_window"
    assert "12 alt reads, 10 alt fragments" in kras_t0_row["cells"][0]["summary"]
    assert kras_t0_row["cells"][1]["state"] == "no_ranked_window"
    assert "source RNA observed (max 2 fragments)" in kras_t0_row["cells"][1]["summary"]
    assert kras_t0_row["cells"][2]["state"] == "different_window_only"
    assert "does not imply absent RNA" in kras_t0_row["cells"][2]["summary"]

    creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked,
        patient_info=pinfo,
        final_review="",
        reviewers="",
        args_for_report={"manufacturability": True, "wt_epitopes": False},
        input_json_file=None,
        products=products,
        selection=reloaded.selection,
        source_evidence=reloaded.selection["source_evidence"],
    )
    tdata = creator.compute_template_data()
    html_path = tmp_path / "longitudinal.html"
    txt_path = tmp_path / "longitudinal.txt"
    make_html_report(tdata, str(html_path))
    make_ascii_report(tdata, str(txt_path))
    for text in (html_path.read_text(), txt_path.read_text()):
        assert "Longitudinal and library support" in text
        assert "[measured zero]" in text
        assert "Qualified only through earlier timepoints" in text or "QUALIFIED ONLY THROUGH EARLIER TIMEPOINTS" in text


def test_matched_normal_observations_conflict_dependent_pipelines_and_replay(tmp_path):
    # Two dependent pipelines on ONE normal DNA library + one blood RNA library
    obs_dna_dragen = MatchedNormalObservation(
        assay="DNA",
        library_id="wgs-blood-normal-1",
        pipeline="DRAGEN",
        pipeline_version="4.2",
        tissue="blood",
        depth=12,
        alt_count=2,
        min_adequate_depth=20,
    )
    obs_dna_onco = MatchedNormalObservation(
        assay="DNA",
        library_id="wgs-blood-normal-1",
        pipeline="oncoanalyser",
        pipeline_version="2.0",
        tissue="blood",
        depth=11,
        alt_count=2,
        min_adequate_depth=20,
    )
    obs_rna_blood = MatchedNormalObservation(
        assay="RNA",
        library_id="rna-blood-1",
        pipeline="STAR-pileup",
        tissue="blood",
        depth=28,
        alt_count=1,
        min_adequate_depth=10,
    )

    grouped = group_normal_observations_by_library([obs_dna_dragen, obs_dna_onco, obs_rna_blood])
    assert len(grouped) == 2
    dna_lib = next(g for g in grouped if g["assay"] == "DNA")
    # Dependent pipelines on one library are NEVER summed (max_observed_alt == 2, not 4!)
    assert dna_lib["dependent_pipeline_count"] == 2
    assert dna_lib["max_observed_alt"] == 2
    assert dna_lib["low_coverage_present"] is True

    isovar_doc = {
        "schema": "isovar.protein_hypotheses.v2",
        "sample_id": "T0-tumor",
        "source": "tumor.bam",
        "events": [
            {
                "event_id": "ev-1",
                "variant": {
                    "reference_name": "GRCh38",
                    "contig": "7",
                    "start": 140753336,
                    "ref": "A",
                    "alt": "T",
                },
                "filters": {"passes_all_filters": True},
                "protein_hypotheses": [
                    {
                        "hypothesis_id": "hyp-1",
                        "amino_acids": "MSIINFEKLAKK",
                        "ends_with_stop_codon": True,
                        "n_terminus": "annotated_start_codon",
                        "contains_mutation": True,
                        "mutation_interval": [1, 10],
                        "isovar_rank": 1,
                        "rna_support": {"reads": 8, "fragments": 6},
                        "gene_names": ["BRAF"],
                        "gene_ids": ["ENSG00000157764"],
                        "transcript_ids": ["ENST00000288602"],
                        "translations": [
                            {
                                "nucleotide_sequence": "ATGAGTATAATCAACTTTGAAAAACTGGCCAAAAAGTAA",
                                "translated_interval": [0, 39],
                                "reference_context": {"translation_table": 1},
                            }
                        ],
                    }
                ],
            }
        ],
    }

    # 1) Default hold_for_review holds out somatic-attested candidate when normal alt reads are observed
    held_ds = antigens_from_isovar(
        isovar_doc,
        somatic=True,
        normal_observations=[obs_dna_dragen, obs_dna_onco, obs_rna_blood],
    )
    (held_ag,) = held_ds.antigens.values()
    assert held_ag.tumor_specificity.status == "held_out"
    assert held_ag.tumor_specificity.requires_review is True
    assert "observed_normal_dna_alt_support" in held_ag.tumor_specificity.rationale_code
    assert "observed_normal_rna_alt_support" in held_ag.tumor_specificity.rationale_code
    assert "observed_normal_alt_at_low_coverage" in held_ag.tumor_specificity.rationale_code
    summary = summarize_normal_assessment(held_ag.tumor_specificity, held_ag.normal_observations)
    assert summary["conflict"] is True
    assert "do not establish germline origin" in summary["note"]

    # 2) JSON envelope load_isovar_input with flag_for_review retains candidate for review & survives replay
    env_path = tmp_path / "isovar_with_normal.json"
    env_path.write_text(json.dumps({
        "schema": "vaxrank.isovar_input.v1",
        "somatic": True,
        "normal_conflict_policy": "flag_for_review",
        "normal_observations": [o.to_dict() for o in (obs_dna_dragen, obs_dna_onco, obs_rna_blood)],
        "result": isovar_doc,
    }))
    flagged_ds = load_isovar_input(env_path)
    ag_json_path = tmp_path / "antigens.json"
    flagged_ds.save(ag_json_path)
    reloaded_ag_ds = AntigenDataset.load(ag_json_path)
    (flagged_ag,) = reloaded_ag_ds.antigens.values()
    assert flagged_ag.tumor_specificity.status == "admitted"
    assert flagged_ag.tumor_specificity.requires_review is True
    assert len(flagged_ag.normal_observations) == 3
    assert flagged_ag.sliced(1, 12).normal_observations == flagged_ag.normal_observations

    # Roundtrip through EpitopeDataset + selected product reports
    ep = _ep("SIINFEKLA", flagged_ag.amino_acids, "pred-braf", offset=1)
    ep_ds = EpitopeDataset.from_predictions(
        [ep], antigens={"pred-braf": flagged_ag}, config=EpitopeConfig(score_expr="1.0")
    )
    ep_path = tmp_path / "braf_epitopes.tsv"
    ep_ds.save(ep_path)
    ep_reloaded = EpitopeDataset.load(ep_path)
    restored_ag = ep_reloaded.antigens["pred-braf"]
    assert restored_ag.normal_observations == flagged_ag.normal_observations

    vp = VaccinePeptide(antigen=restored_ag.sliced(1, 12), epitopes=[ep.sliced(1, 12)])
    ranked = [(restored_ag, [vp])]
    products = {
        "peptide": assemble_peptide_constructs(
            ranked,
            PeptideConstructConfig(min_antigen_length_aa=11, max_antigen_length_aa=11),
            target_alleles=["HLA-A*02:01"],
        )
    }
    creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked,
        patient_info=PatientInfo(
            patient_id="P-NORMAL",
            vcf_paths=[],
            bam_path=None,
            mhc_alleles=["HLA-A*02:01"],
            num_somatic_variants=1,
            num_coding_effect_variants=1,
            num_variants_with_rna_support=1,
            num_variants_with_vaccine_peptides=1,
        ),
        final_review="",
        reviewers="",
        args_for_report={"manufacturability": True, "wt_epitopes": False},
        input_json_file=None,
        products=products,
        candidate_details=True,
    )
    tdata = creator.compute_template_data()
    html_out = tmp_path / "normal_report.html"
    txt_out = tmp_path / "normal_report.txt"
    make_html_report(tdata, str(html_out))
    make_ascii_report(tdata, str(txt_out))
    for rendered in (html_out.read_text(), txt_out.read_text()):
        assert "wgs-blood-normal-1" in rendered
        assert "DRAGEN" in rendered and "oncoanalyser" in rendered
        assert "rna-blood-1" in rendered
        assert "do not establish germline origin" in rendered

    # 3) RNA-only normal observation leaves normal DNA explicitly missing and requires review
    rna_only_obs = MatchedNormalObservation(
        assay="RNA",
        library_id="rna-blood-clean",
        pipeline="STAR-pileup",
        tissue="blood",
        depth=30,
        alt_count=0,
        min_adequate_depth=10,
    )
    rna_only_ds = antigens_from_isovar(
        isovar_doc,
        somatic=True,
        normal_observations=[rna_only_obs],
    )
    (rna_only_ag,) = rna_only_ds.antigens.values()
    assert rna_only_ag.tumor_specificity.status == "admitted"
    assert rna_only_ag.tumor_specificity.requires_review is True
    assert "matched_normal_dna_missing" in rna_only_ag.tumor_specificity.rationale_code


def test_source_outcomes_omits_not_cta_background_and_keeps_concise_pipeline_headers(tmp_path):
    """Regression test for #887: not_cta background rows are summarized, not dumped verbatim."""
    from vaxrank.input_coverage import coverage_view

    long_source_uri = (
        "rna-seq/reprocessed/BG003082/BG003082.rsem.genes.results:"
        "osteosarc:cbbf4b9d834a8a25faf55897b31a8b21a3b48fe670586f50332e16a8b34d93fc"
    )
    ag_cta = _cta_antigen(
        "ENSG00000241476",
        "SSX2",
        "ENST00000337502",
        "KMKASEKIFYVYMKRKYEAMTKLGF",
    )
    cta_res = CTAExpressionResult(
        input_contract=CTAExpressionInput(
            measurement_level="gene",
            id_column="gene_id",
            value_column="TPM",
            expression_unit="TPM",
            sample_id="T0-BG",
            evidence_source=long_source_uri,
            evidence_version="1.0",
            assay="bulk_rna_seq",
        ),
        admission_policy=CTAAdmissionPolicy(5.0, "TPM"),
        exclusion_policy=CTATargetExclusionPolicy(gene_patterns=(), gene_exceptions=()),
        input_sha256="b" * 64,
        reference_assembly="GRCh38",
        annotation_name="ensembl",
        annotation_version="112",
        decisions=(
            CTAExpressionDecision(
                input_identifier="ENSG00000241476",
                gene_name="SSX2",
                gene_id="ENSG00000241476",
                transcript_id="ENST00000337502",
                value=24.5,
                status="admitted",
                assessment=CTASequenceAssessment(
                    antigen=ag_cta,
                    bundle_id="bundle-1",
                    record_id="rec-1",
                    occurrence={"gene_id": "ENSG00000241476"},
                ),
            ),
            CTAExpressionDecision(
                input_identifier="ENSG00000185686",
                gene_name="PRAME",
                gene_id="ENSG00000185686",
                transcript_id="ENST00000398743",
                value=1.2,
                status="held_out",
            ),
            *(
                CTAExpressionDecision(
                    input_identifier=f"ENSG00000900{i:03d}",
                    gene_name=f"HOUSEKEEPING_{i}",
                    gene_id=f"ENSG00000900{i:03d}",
                    transcript_id=f"ENST00000900{i:03d}",
                    value=10.0,
                    status="not_cta",
                )
                for i in range(50)
            ),
        ),
    )

    source_evidence = {"cta:T0-BG": {"cta_expression_admission": to_native_json(cta_res)}}
    cov = coverage_view({}, source_evidence, [], {})
    assert len(cov["expression"]) == 1
    expr_entry = cov["expression"][0]
    assert expr_entry["omitted_not_cta_count"] == 50
    assert [d["gene_name"] for d in expr_entry["decisions"]] == ["SSX2", "PRAME"]

    ep = _ep("KMKASEKIF", ag_cta.amino_acids, "pred-ssx2", offset=0)
    vp = VaccinePeptide(antigen=ag_cta, epitopes=[ep])
    ranked = [(ag_cta, [vp])]
    pinfo = PatientInfo(
        patient_id="P-887",
        vcf_paths=[],
        bam_path=None,
        mhc_alleles=["HLA-A*02:01"],
        num_somatic_variants=0,
        num_coding_effect_variants=0,
        num_variants_with_rna_support=0,
        num_variants_with_vaccine_peptides=0,
        input_coverage={"sources": [], "original_inputs": [], "current_inputs": []},
    )
    creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked,
        patient_info=pinfo,
        final_review="",
        reviewers="",
        args_for_report={"manufacturability": True, "wt_epitopes": False},
        input_json_file=None,
        source_evidence=source_evidence,
    )
    tdata = creator.compute_template_data()
    html_out = tmp_path / "report_887.html"
    txt_out = tmp_path / "report_887.txt"
    make_html_report(tdata, str(html_out))
    make_ascii_report(tdata, str(txt_out))
    html_text = html_out.read_text()
    txt_text = txt_out.read_text()
    assert "50 non-CTA background features omitted" in html_text
    assert "Omitted 50 non-CTA background features" in txt_text
    assert "HOUSEKEEPING_0" not in html_text
    assert "HOUSEKEEPING_0" not in txt_text
    assert "PRAME" in html_text and "held_out" in html_text
    assert long_source_uri not in html_text
    assert "bulk_rna_seq" in html_text
    assert "<h3>CTA expression inputs</h3>" in html_text
    assert "CTA expression input:" not in html_text


def test_many_normal_libraries_render_concise_headline_and_collapsible_audit(tmp_path):
    """Regression test for #892: avoid multi-kilobyte inline normal walls in product cells."""
    from dataclasses import replace

    obs = [
        # Duplicate identical runs on same library (different source_identity) collapse in display
        MatchedNormalObservation(
            assay="DNA", library_id="wgs-1", pipeline="DRAGEN", source_identity="pair-a",
            tissue="blood", depth=52, alt_count=0,
        ),
        MatchedNormalObservation(
            assay="DNA", library_id="wgs-1", pipeline="DRAGEN", source_identity="pair-b",
            tissue="blood", depth=52, alt_count=0,
        ),
        MatchedNormalObservation(
            assay="DNA", library_id="wes-1", pipeline="Sarek", source_identity="pair-a",
            tissue="blood", depth=143, alt_count=0,
        ),
        *(
            MatchedNormalObservation(
                assay="RNA", library_id=f"blood-pool-{i}", pipeline="catalogue",
                source_identity=f"pool-{i}", tissue="blood", depth=2000 + i * 100,
                alt_count=2 if i in (1, 3) else 0,
            )
            for i in range(5)
        ),
    ]
    ag = VaccineAntigen(
        kind="mutation",
        amino_acids="MTSAKKMLLWNWKKWEKMRNK",
        targetable_mask=TargetableMask((AminoAcidInterval(13, 14),)),
        tumor_specificity=TumorSpecificityAttestation(
            status="admitted",
            evidence_kind="somatic_variant",
            evidence_source="somatic_vcf",
            patient_specific=True,
            rationale_code="somatic",
        ),
        gene_name="MBD4",
        gene_id="ENSG00000129071",
        source_identifier="mbd4-ncorf",
        normal_observations=tuple(obs),
    )
    ep = _ep("KKMLLWNWKK", ag.amino_acids, "pred-mbd4", offset=4)
    vp = VaccinePeptide(antigen=ag, epitopes=[ep])
    ranked = [(ag, [vp])]
    products = {
        "peptide": assemble_peptide_constructs(
            ranked,
            PeptideConstructConfig(min_antigen_length_aa=15, max_antigen_length_aa=25),
            target_alleles=["HLA-A*02:01"],
        )
    }
    creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked,
        patient_info=PatientInfo(
            patient_id="P-892",
            vcf_paths=[],
            bam_path=None,
            mhc_alleles=["HLA-A*02:01"],
            num_somatic_variants=1,
            num_coding_effect_variants=1,
            num_variants_with_rna_support=1,
            num_variants_with_vaccine_peptides=1,
        ),
        final_review="",
        reviewers="",
        args_for_report={"manufacturability": False, "wt_epitopes": False},
        input_json_file=None,
        products=products,
    )
    tdata = creator.compute_template_data()
    html_out = tmp_path / "report_892.html"
    txt_out = tmp_path / "report_892.txt"
    make_html_report(tdata, str(html_out))
    make_ascii_report(tdata, str(txt_out))
    html_text = html_out.read_text()
    txt_text = txt_out.read_text()
    assert "Matched normal: DNA 0 alt across 2/2 libraries (52-143x); RNA alt reads in 2/5 libraries (max 2/2100). Review required." in html_text
    assert "Matched normal: DNA 0 alt across 2/2 libraries (52-143x); RNA alt reads in 2/5 libraries (max 2/2100). Review required." in txt_text
    assert "Per-library normal audit (7 libraries)" in html_text
    assert "DRAGEN 0/52, DRAGEN 0/52" not in html_text
    assert "observed_alt_support" not in html_text
    assert "[adequate_coverage]" not in html_text
    # Ensure replace import is used cleanly
    assert replace(ag, gene_name="MBD4").gene_name == "MBD4"
