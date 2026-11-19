"""Tests for within-library Isovar RNA read-set unions on shared vaccine sequences (#825)."""

import json

import pytest
from isovar import RNAFeature, RNASamplingGeometry, RNASamplingProfile, measure_rna_feature
from isovar.rna_path_linkage import OriginalPathMatch
from mhctools import Prediction

from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.cli.arg_parser import external_input_arg_parser
from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
from vaxrank.combined_score_dsl import evaluate_combined_score
from vaxrank.core_logic import SOURCE_AGNOSTIC_RANKING_RULES
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset, dataset_ranking_result
from vaxrank.external_input import ExternalConstructOptions
from vaxrank.external_report import ExternalReport
from vaxrank.patient_info import PatientInfo
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.report import make_ascii_report, make_html_report, TemplateDataCreator
from vaxrank.rna_support import RNAFeatureEvidence, RNASupportArchive, UNION_REF_FIELD
from vaxrank.shared_sequence_rna import (
    attach_shared_sequence_rna_support,
    compute_shared_sequence_rna_unions,
)
from vaxrank.vaccine_antigen import (
    AminoAcidInterval,
    TargetableMask,
    TumorSpecificityAttestation,
    VaccineAntigen,
)
from vaxrank.vaccine_config import VaccineConfig
from vaxrank.vaccine_peptide import VaccinePeptide


def _profile(lib_reads=2_000_000, lib_fragments=1_000_000, paired=True, template_len=120):
    segments = ((0, 45), (55, 100)) if paired else ((0, 100),)
    return RNASamplingProfile(
        "original-library",
        (RNASamplingGeometry(segments, template_len),),
        lib_reads,
        lib_fragments,
        "whole_library",
    )


def _measurement(
    feature_id,
    matches,
    *,
    sample="T0",
    library="lib-T0",
    interval=(100, 130),
    fragment_mode="union",
    profile=None,
):
    feature = RNAFeature(
        feature_id,
        "path-" + feature_id,
        (interval,),
        "coding_interval",
        fragment_mode=fragment_mode,
    )
    return measure_rna_feature(
        feature,
        matches,
        evidence_scope=(sample, library),
        profile=_profile() if profile is None else profile,
    )


def _antigen(gene_name, gene_id, transcript_id, sequence, rna_support=()):
    return VaccineAntigen(
        kind="CTA",
        amino_acids=sequence,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status="admitted",
            evidence_kind="synthetic_fixture",
            evidence_source="test",
            patient_specific=True,
        ),
        gene_name=gene_name,
        gene_id=gene_id,
        transcript_ids=(transcript_id,),
        protein_ids=("ENSP_" + gene_name,),
        species="homo_sapiens",
        source_identifier=gene_id + ":" + transcript_id,
        rna_support=tuple(rna_support),
        source_protein_completeness="complete",
    )


def test_overlapping_disjoint_and_repeated_alignment_unions():
    seq = "SIINFEKLAKK"
    # Overlapping sets: Gene A {q1, q2, q3, q4}, Gene B {q3, q4, q5} -> union = 5 (not 7, not 4)
    m_a = _measurement(
        "feat-a",
        [OriginalPathMatch(("RG", f"q{i}", 0), f"p-a-{i}", 80, 160) for i in (1, 2, 3, 4)],
    )
    m_b = _measurement(
        "feat-b",
        [OriginalPathMatch(("RG", f"q{i}", 0), f"p-b-{i}", 80, 160) for i in (3, 4, 5)],
    )
    ag_a = _antigen(
        "CTAG1A", "ENSG00000184033", "ENST_A", seq,
        (RNAFeatureEvidence(m_a, "window", (0, len(seq))),),
    )
    ag_b = _antigen(
        "CTAG1B", "ENSG00000184034", "ENST_B", seq,
        (RNAFeatureEvidence(m_b, "window", (0, len(seq))),),
    )

    (union_rec,) = compute_shared_sequence_rna_unions([ag_a, ag_b])
    vals = union_rec["values"]
    assert vals["reads"] == len(set(union_rec["read_ids"])) == 5
    assert vals["fragments"] == len(set(union_rec["fragment_ids"])) == 5
    assert vals["sum_of_source_fragments"] == 7
    assert vals["fragments_per_million"] == pytest.approx(5.0)
    assert vals["reads_per_million"] == pytest.approx(2.5)
    assert vals["effective_length_fragments"] is None
    assert union_rec["unavailable"]["fragments_per_effective_kb_per_million"] == (
        "shared_sequence_effective_length_undefined"
    )
    assert "shared_sequence_ambiguous_origin_attribution" in union_rec["attribution_limitations"]
    assert "multi_source_shared_read_or_fragment_identities" in union_rec["attribution_limitations"]
    members = {m["gene_name"]: m for m in union_rec["members"]}
    assert members["CTAG1A"]["fragments_shared_with_other_sources"] == 2
    assert members["CTAG1A"]["fragments_exclusive_to_source"] == 2
    assert members["CTAG1B"]["fragments_shared_with_other_sources"] == 2
    assert members["CTAG1B"]["fragments_exclusive_to_source"] == 1

    # Disjoint sets: {d1, d2} and {d3, d4, d5} -> union = 5, shared = 0
    m_d1 = _measurement(
        "dis-1", [OriginalPathMatch(("RG", f"d{i}", 0), f"pd-{i}", 80, 160) for i in (1, 2)]
    )
    m_d2 = _measurement(
        "dis-2", [OriginalPathMatch(("RG", f"d{i}", 0), f"pd-{i}", 80, 160) for i in (3, 4, 5)]
    )
    (disjoint_rec,) = compute_shared_sequence_rna_unions([
        _antigen(
            "GENE1", "ENSG00000100001", "ENST_1", seq,
            (RNAFeatureEvidence(m_d1, "window", (0, len(seq))),),
        ),
        _antigen(
            "GENE2", "ENSG00000100002", "ENST_2", seq,
            (RNAFeatureEvidence(m_d2, "window", (0, len(seq))),),
        ),
    ])
    assert disjoint_rec["values"]["fragments"] == 5
    assert all(m["fragments_shared_with_other_sources"] == 0 for m in disjoint_rec["members"])
    assert (
        "multi_source_shared_read_or_fragment_identities"
        not in disjoint_rec["attribution_limitations"]
    )

    # Repeated multi-mapping alignments of exact same 3 reads across two paralogs -> union = 3 (not 6)
    m_r1 = _measurement(
        "rep-1", [OriginalPathMatch(("RG", f"r{i}", 0), f"pr1-{i}", 80, 160) for i in (1, 2, 3)]
    )
    m_r2 = _measurement(
        "rep-2", [OriginalPathMatch(("RG", f"r{i}", 0), f"pr2-{i}", 80, 160) for i in (1, 2, 3)]
    )
    (rep_rec,) = compute_shared_sequence_rna_unions([
        _antigen(
            "PARA1", "ENSG00000100011", "ENST_11", seq,
            (RNAFeatureEvidence(m_r1, "window", (0, len(seq))),),
        ),
        _antigen(
            "PARA2", "ENSG00000100012", "ENST_12", seq,
            (RNAFeatureEvidence(m_r2, "window", (0, len(seq))),),
        ),
    ])
    assert rep_rec["values"]["reads"] == 3
    assert rep_rec["values"]["fragments"] == 3
    assert rep_rec["values"]["sum_of_source_fragments"] == 6


def test_joint_paired_mate_fragments_preserved_and_promoted_when_spanned():
    seq = "SIINFEKLAKK"
    paired_prof = _profile(paired=True)
    # Fragment 'joint_only' spans [100, 130) ONLY via Mate 1 [95, 118) U Mate 2 [115, 140)
    m_joint = _measurement(
        "joint-feat",
        [
            OriginalPathMatch(("RG", "joint_only", 64), "m1", 95, 118),
            OriginalPathMatch(("RG", "joint_only", 128), "m2", 115, 140),
        ],
        interval=(100, 130),
        profile=paired_prof,
    )
    assert m_joint["values"]["reads"] == 0
    assert m_joint["values"]["fragments"] == 1
    assert len(m_joint["joint_only_fragment_ids"]) == 1

    # Fragment 'single_span' spans [100, 130) via Mate 1 [90, 140) alone
    m_single = _measurement(
        "single-feat",
        [OriginalPathMatch(("RG", "single_span", 64), "s1", 90, 140)],
        interval=(100, 130),
        profile=paired_prof,
    )
    assert m_single["values"]["reads"] == 1
    assert m_single["values"]["fragments"] == 1

    (union_rec,) = compute_shared_sequence_rna_unions([
        _antigen(
            "G1", "ENSG00000100021", "ENST_21", seq,
            (RNAFeatureEvidence(m_joint, "window", (0, len(seq))),),
        ),
        _antigen(
            "G2", "ENSG00000100022", "ENST_22", seq,
            (RNAFeatureEvidence(m_single, "window", (0, len(seq))),),
        ),
    ])
    # Union has 1 individually spanning read, 2 fragments (1 joint-only + 1 single-read)
    assert union_rec["values"]["reads"] == 1
    assert union_rec["values"]["fragments"] == 2
    assert union_rec["joint_only_fragment_ids"] == m_joint["joint_only_fragment_ids"]

    # If another isoform's feature is spanned individually by Mate 1 of 'joint_only',
    # then 'joint_only' is promoted out of joint_only_fragment_ids in the union.
    m_spanned = _measurement(
        "spanned-feat",
        [OriginalPathMatch(("RG", "joint_only", 64), "m1-span", 95, 135)],
        interval=(100, 130),
        profile=paired_prof,
    )
    (promoted_rec,) = compute_shared_sequence_rna_unions([
        _antigen(
            "G1", "ENSG00000100021", "ENST_21", seq,
            (RNAFeatureEvidence(m_joint, "window", (0, len(seq))),),
        ),
        _antigen(
            "G3", "ENSG00000100023", "ENST_23", seq,
            (RNAFeatureEvidence(m_spanned, "window", (0, len(seq))),),
        ),
    ])
    assert promoted_rec["values"]["reads"] == 1
    assert promoted_rec["values"]["fragments"] == 1
    assert promoted_rec["joint_only_fragment_ids"] == []


def test_cross_library_separation_and_conflicting_denominator_rejection():
    seq = "SIINFEKLAKK"
    m_t0 = _measurement(
        "f-t0", [OriginalPathMatch(("RG", "q1", 0), "p1", 80, 160)], sample="T0", library="lib-T0"
    )
    m_t2 = _measurement(
        "f-t2", [OriginalPathMatch(("RG", "q1", 0), "p2", 80, 160)], sample="T2", library="lib-T2"
    )
    ag_t0 = _antigen(
        "G1", "ENSG00000100031", "ENST_31", seq,
        (RNAFeatureEvidence(m_t0, "window", (0, len(seq))),),
    )
    ag_t2 = _antigen(
        "G2", "ENSG00000100032", "ENST_32", seq,
        (RNAFeatureEvidence(m_t2, "window", (0, len(seq))),),
    )

    unions = compute_shared_sequence_rna_unions([ag_t0, ag_t2])
    assert [u["evidence_scope"] for u in unions] == [["T0", "lib-T0"], ["T2", "lib-T2"]]
    assert [u["values"]["fragments"] for u in unions] == [1, 1]

    vp = VaccinePeptide(
        antigen=ag_t0,
        epitopes=[],
        combined_score_expr="rna_shared_fragments",
        shared_sequence_rna_support=unions,
    )
    with pytest.raises(ValueError, match="Unavailable RNA score bindings"):
        _ = vp.combined_score

    # Conflicting whole-library denominators in the same evidence_scope must raise
    m_conflict = _measurement(
        "f-conflict",
        [OriginalPathMatch(("RG", "q2", 0), "p3", 80, 160)],
        sample="T0",
        library="lib-T0",
        profile=_profile(lib_reads=2_000_000, lib_fragments=1_500_000),
    )
    ag_conflict = _antigen(
        "G3", "ENSG00000100033", "ENST_33", seq,
        (RNAFeatureEvidence(m_conflict, "window", (0, len(seq))),),
    )
    with pytest.raises(ValueError, match="Conflicting library_fragments denominators"):
        compute_shared_sequence_rna_unions([ag_t0, ag_conflict])

    # Conflicting fragment_mode in the same evidence_scope must raise
    m_any = _measurement(
        "f-any",
        [OriginalPathMatch(("RG", "q3", 0), "p4", 80, 160)],
        sample="T0",
        library="lib-T0",
        fragment_mode="any_read",
    )
    ag_any = _antigen(
        "G4", "ENSG00000100034", "ENST_34", seq,
        (RNAFeatureEvidence(m_any, "window", (0, len(seq))),),
    )
    with pytest.raises(ValueError, match="conflicting fragment_mode"):
        compute_shared_sequence_rna_unions([ag_t0, ag_any])


def test_opt_in_shared_sequence_rna_policy_ranking_and_frozen_default(tmp_path):
    # Full proteins (20 aa) sliced into 11-aa windows via max_windows_per_antigen_source=3!
    shared_window = "SIINFEKLAKK"
    solo_window = "GILGFVFTLKK"
    full_p1 = "AAAA" + shared_window + "GGGGG"
    full_p2 = "TTTT" + shared_window + "CCCCC"
    full_solo = "MMMM" + solo_window + "YYYYY"

    # Paralog 1 (4 fragments) and Paralog 2 (3 fragments) share 1 fragment -> union = 6 (naive sum = 7), solo = 5
    m_p1 = _measurement(
        "p1",
        [
            OriginalPathMatch(("RG", "shared_q", 0), "x0", 80, 160),
            OriginalPathMatch(("RG", "p1_1", 0), "x1", 80, 160),
            OriginalPathMatch(("RG", "p1_2", 0), "x2", 80, 160),
            OriginalPathMatch(("RG", "p1_3", 0), "x3", 80, 160),
        ],
    )
    m_p2 = _measurement(
        "p2",
        [
            OriginalPathMatch(("RG", "shared_q", 0), "y0", 80, 160),
            OriginalPathMatch(("RG", "p2_1", 0), "y1", 80, 160),
            OriginalPathMatch(("RG", "p2_2", 0), "y2", 80, 160),
        ],
    )
    m_solo = _measurement(
        "solo", [OriginalPathMatch(("RG", f"s_{i}", 0), f"z{i}", 80, 160) for i in range(5)]
    )

    ag_p1 = _antigen(
        "PARA1", "ENSG00000100041", "ENST_41", full_p1,
        (RNAFeatureEvidence(m_p1, "window", (4, 15)),),
    )
    ag_p2 = _antigen(
        "PARA2", "ENSG00000100042", "ENST_42", full_p2,
        (RNAFeatureEvidence(m_p2, "window", (4, 15)),),
    )
    ag_solo = _antigen(
        "SOLO", "ENSG00000100043", "ENST_43", full_solo,
        (RNAFeatureEvidence(m_solo, "window", (4, 15)),),
    )

    def _ep(ag, pid):
        pep = ag.amino_acids[4:15]
        return CandidateEpitope(
            sequence=pep,
            source_sequence=ag.amino_acids,
            offset=4,
            prediction_id=pid,
            occurs_in_reference=False,
            overlaps_targetable=True,
            predictions=[
                Prediction(
                    peptide=pep,
                    allele="HLA-A*02:01",
                    kind="pMHC_affinity",
                    value=25.0,
                    score=0.9,
                    percentile_rank=0.2,
                    predictor_name="mock",
                    predictor_version="1.0",
                )
            ],
            per_allele_scores={"HLA-A*02:01": 1.0},
        )

    ds = EpitopeDataset.from_predictions(
        [_ep(ag_p1, "pred-p1"), _ep(ag_p2, "pred-p2"), _ep(ag_solo, "pred-solo")],
        antigens={"pred-p1": ag_p1, "pred-p2": ag_p2, "pred-solo": ag_solo},
        config=EpitopeConfig(score_expr="1.0"),
    )
    saved_path = tmp_path / "candidates.tsv"
    ds.save(saved_path)
    loaded = EpitopeDataset.load(saved_path)
    report = ExternalReport(
        "epitopes", str(saved_path),
        report_df=loaded.report_frame(),
        epitopes=loaded.epitopes,
        dataset=loaded,
    )

    # 1) Default source policy: SOLO (5 fragments) outranks PARA1 (4) and PARA2 (3)
    default_vc = VaccineConfig(
        min_peptide_length=11,
        preferred_peptide_length=11,
        max_peptide_length=11,
        combined_score_expr="sqrt(rna_fragments) * target_epitope_score",
        ranking_rules=SOURCE_AGNOSTIC_RANKING_RULES,
    )
    assert default_vc.shared_sequence_rna_policy == "source"
    default_ranked = dataset_ranking_result(
        report, list(loaded.epitopes), options=ExternalConstructOptions(vaccine_config=default_vc)
    ).ranked
    assert [src.gene_name for src, _ in default_ranked] == ["SOLO", "PARA1", "PARA2"]

    # Using rna_shared_fragments under default 'source' policy raises explicitly
    bad_vc = VaccineConfig(
        min_peptide_length=11,
        preferred_peptide_length=11,
        max_peptide_length=11,
        combined_score_expr="sqrt(rna_shared_fragments) * target_epitope_score",
        ranking_rules=SOURCE_AGNOSTIC_RANKING_RULES,
    )
    with pytest.raises(ValueError, match="shared_sequence_rna_policy"):
        dataset_ranking_result(
            report, list(loaded.epitopes), options=ExternalConstructOptions(vaccine_config=bad_vc)
        )

    # 2) Explicit CLI + VaccineConfig opt-in: PARA1 & PARA2 (union = 6 fragments, not 7!) outrank SOLO (5)
    parsed = external_input_arg_parser().parse_args([
        "--shared-sequence-rna-policy", "within_library_union",
        "--vaccine-peptide-length", "11",
    ])
    union_vc = vaccine_config_from_args(
        parsed,
        merged_config={
            "vaccine_peptides": {
                "combined_score_expr": "sqrt(rna_shared_fragments) * target_epitope_score",
                "ranking_rules": list(SOURCE_AGNOSTIC_RANKING_RULES),
            }
        },
    )
    assert union_vc.shared_sequence_rna_policy == "within_library_union"
    union_ranked = dataset_ranking_result(
        report, list(loaded.epitopes), options=ExternalConstructOptions(vaccine_config=union_vc)
    ).ranked
    assert [src.gene_name for src, _ in union_ranked] == ["PARA1", "PARA2", "SOLO"]
    assert union_ranked[0][1][0].combined_score == pytest.approx(6.0 ** 0.5)
    assert union_ranked[2][1][0].combined_score == pytest.approx(5.0 ** 0.5)

    # Effective-length group bindings remain explicitly unavailable even under union policy
    attached = attach_shared_sequence_rna_support(
        [vp for _, vps in union_ranked for vp in vps], policy="within_library_union"
    )
    with pytest.raises(ValueError, match="shared_sequence_effective_length_undefined"):
        evaluate_combined_score("rna_shared_fragments_per_effective_kb_per_million", attached[0])

    # RNASupportArchive deduplicates repeated shared_sequence_rna_union records in-file
    archive = RNASupportArchive()
    encoded = archive.encode_tree([
        union_ranked[0][1][0].shared_sequence_rna_support,
        union_ranked[1][1][0].shared_sequence_rna_support,
    ])
    assert len(archive.unions) == 1
    assert UNION_REF_FIELD in encoded[0][0]
    payload = {}
    archive.write_to_payload(payload)
    restored = RNASupportArchive.from_payload(json.loads(json.dumps(payload)))
    assert list(restored.decode_tree(encoded)[0]) == list(union_ranked[0][1][0].shared_sequence_rna_support)

    # Reports surface shared-sequence RNA union counts, naive sum, and attribution limitations
    peptides = assemble_peptide_constructs(
        union_ranked,
        PeptideConstructConfig(min_antigen_length_aa=11, max_antigen_length_aa=11),
        target_alleles=["HLA-A*02:01"],
    )
    creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=union_ranked,
        patient_info=PatientInfo(
            patient_id="P-SHARED-RNA",
            vcf_paths=[],
            bam_path=None,
            mhc_alleles=["HLA-A*02:01"],
            num_somatic_variants=3,
            num_coding_effect_variants=3,
            num_variants_with_rna_support=3,
            num_variants_with_vaccine_peptides=3,
        ),
        final_review="",
        reviewers=[],
        args_for_report={"manufacturability": True, "wt_epitopes": False},
        input_json_file=None,
        products={"peptide": peptides},
        candidate_details=True,
    )
    template_data = creator.compute_template_data()
    html_path = tmp_path / "report.html"
    ascii_path = tmp_path / "report.txt"
    make_html_report(template_data, str(html_path))
    make_ascii_report(template_data, str(ascii_path))
    for rendered in (html_path.read_text(), ascii_path.read_text()):
        assert "Within-library RNA union: T0/lib-T0: 6 reads, 6 fragments (naive sum 7 fragments" in rendered
        assert "Shared-sequence RNA union 1" in rendered
        assert "shared_sequence_ambiguous_origin_attribution" in rendered
