# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for AI-ready tumor_board_dossier.json, self-contained interactive
HTML review, and deterministic --selection-overrides replay
(openvax/vaxrank#922)."""

import json
from dataclasses import dataclass
from types import SimpleNamespace

import pytest
from jinja2 import Environment, PackageLoader
from mhctools.pred import Prediction

from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.cli.arg_parser import parse_vaxrank_args
from vaxrank.config import load_vaxrank_config
from vaxrank.mhc_allele_status import (
    APMStatus,
    MHCAlleleStatus,
    MHCAlleleStatusBundle,
    apply_mhc_allele_status_to_ranked,
)
from vaxrank.tumor_board import (
    DOSSIER_SCHEMA,
    OVERRIDES_SCHEMA,
    TumorBoardDecision,
    TumorBoardOverrides,
    TumorBoardPolicy,
    TumorBoardValidationError,
    build_tumor_board_dossier,
    load_tumor_board_overrides,
    validate_and_apply_overrides,
    write_tumor_board_dossier,
)
from vaxrank.vaccine_peptide import VaccinePeptide


@dataclass(frozen=True)
class _FakeSource:
    short_description: str
    gene_names: tuple[str, ...]
    contig: str = "17"
    original_start: int = 7675088
    original_ref: str = "C"
    original_alt: str = "T"


def _make_candidate(
    name,
    seq,
    allele_scores,
    *,
    n_alt_reads=20,
    occurs_in_reference=False,
    proxy_tpm=None,
):
    preds = tuple(
        Prediction(
            kind="pMHC_affinity",
            peptide=seq[:9],
            allele=allele,
            score=max(score, 0.5),
            percentile_rank=0.4,
        )
        for allele, score in allele_scores.items()
    )
    ep = CandidateEpitope(
        sequence=seq[:9],
        source_sequence=seq,
        offset=0,
        predictions=preds,
        overlaps_mutation=not occurs_in_reference,
        overlaps_targetable=not occurs_in_reference,
        occurs_in_reference=occurs_in_reference,
        per_allele_scores=dict(allele_scores),
    )
    fragment = SimpleNamespace(
        amino_acids=seq,
        gene_name=name.split("_")[0],
        mutant_amino_acid_start_offset=2,
        mutant_amino_acid_end_offset=3,
        n_alt_reads=n_alt_reads,
        n_ref_reads=10,
        n_overlapping_reads=n_alt_reads + 10,
        n_alt_reads_supporting_protein_sequence=n_alt_reads,
        supporting_reference_transcripts=(),
        proxy_expression_tpm=proxy_tpm,
        proxy_expression_cancer_code="LUAD" if proxy_tpm is not None else "",
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=[ep],
        combined_score_expr="sqrt(max(n_alt_reads, 1)) * target_epitope_score",
    )
    source = _FakeSource(
        short_description=name,
        gene_names=(name.split("_")[0],),
    )
    return (source, [vp])


def test_tumor_board_policy_validation_and_yaml(tmp_path):
    pol = TumorBoardPolicy()
    assert pol.forbid_override_exact_unsafe_self is True
    assert pol.max_normal_alt_vaf == 0.02
    assert pol.require_dossier_sha256_match is True
    assert pol.require_citation_keys is True

    with pytest.raises(ValueError, match="max_normal_alt_vaf"):
        TumorBoardPolicy(max_normal_alt_vaf=1.5)
    with pytest.raises(ValueError, match="min_rna_alt_reads_review_threshold"):
        TumorBoardPolicy(min_rna_alt_reads_review_threshold=-1)
    with pytest.raises(ValueError, match="Unknown tumor_board policy"):
        TumorBoardPolicy.resolve({"bogus": 1})

    cfg_path = tmp_path / "tb_policy.yaml"
    cfg_path.write_text(
        "tumor_board:\n"
        "  max_normal_alt_vaf: 0.01\n"
        "  max_gravy_review_threshold: 1.8\n"
        "  min_rna_alt_reads_review_threshold: 5\n"
        "  require_citation_keys: true\n"
    )
    loaded = load_vaxrank_config(config_path=str(cfg_path))
    resolved = TumorBoardPolicy.resolve(loaded["tumor_board"])
    assert resolved.max_normal_alt_vaf == 0.01
    assert resolved.max_gravy_review_threshold == 1.8
    assert resolved.min_rna_alt_reads_review_threshold == 5


def test_build_dossier_6_layers_hard_gates_and_review_flags(tmp_path):
    # Candidate 1: clean, strong RNA, retained allele
    c1 = _make_candidate("TP53_p.R175H", "EVVRHCPHHERCSDSDGLAPPQHLI", {"HLA-B*07:02": 4.0}, n_alt_reads=30)
    # Candidate 2: normal alt leakage hard gate
    c2 = _make_candidate("ARTIFACT_p.A1V", "ACDEFGHIKLMNPQRSTVWYACDEF", {"HLA-B*07:02": 3.0}, n_alt_reads=15)
    # Candidate 3: all restricting MHC alleles lost (cross-species / HLA-A*02:01 LOH)
    c3 = _make_candidate("LOHONLY_p.L2P", "WYRSTVQPMLKIHGFEDCAWQRSTV", {"HLA-A*02:01": 5.0}, n_alt_reads=20)
    # Candidate 4: DNA-only proxy + hydrophobic -> Review flags (no hard gate)
    c4 = _make_candidate("PROXY_p.V5I", "IIIIIIILLLLLLVVVVVVAAAAAA", {"HLA-B*07:02": 2.0}, n_alt_reads=0, proxy_tpm=12.5)
    # Candidate 5: exact self-proteome only -> hard gate
    c5 = _make_candidate("SELF_p.G1G", "MNPQRSTVWYACDEFGHIKLMNPQR", {"HLA-B*07:02": 2.0}, n_alt_reads=10, occurs_in_reference=True)

    ranked = [c1, c2, c3, c4, c5]
    mhc_bundle = MHCAlleleStatusBundle(
        alleles=(
            MHCAlleleStatus(allele="HLA-A*02:01", dna_loh=True, dna_copy_number=0.0),
            MHCAlleleStatus(allele="HLA-B*07:02", status="retained"),
        ),
        apm=APMStatus(b2m_status="intact"),
    )
    ranked, mhc_audit = apply_mhc_allele_status_to_ranked(ranked, mhc_bundle)

    dossier = build_tumor_board_dossier(
        ranked,
        patient_id="PT-001",
        mhc_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        mhc_allele_status=mhc_audit,
        dna_vaf_by_variant={c2[0]: {"tumor_vaf": 0.25, "normal_alt_vaf": 0.08}},
    )
    assert dossier["schema"] == DOSSIER_SCHEMA
    assert len(dossier["dossier_sha256"]) == 64
    assert dossier["candidate_count"] == 5
    assert dossier["hard_gate_failed_count"] == 3

    by_id = {c["candidate_id"]: c for c in dossier["candidates"]}
    assert by_id["TP53_p.R175H"]["default_recommendation"] == "Accept"
    assert by_id["TP53_p.R175H"]["hard_gate_failures"] == []

    assert by_id["ARTIFACT_p.A1V"]["default_recommendation"] == "Reject"
    assert [g["gate"] for g in by_id["ARTIFACT_p.A1V"]["hard_gate_failures"]] == ["normal_alt_leakage"]

    assert by_id["LOHONLY_p.L2P"]["default_recommendation"] == "Reject"
    assert [g["gate"] for g in by_id["LOHONLY_p.L2P"]["hard_gate_failures"]] == [
        "all_restricting_mhc_lost_or_silenced"
    ]

    assert by_id["SELF_p.G1G"]["default_recommendation"] == "Reject"
    assert [g["gate"] for g in by_id["SELF_p.G1G"]["hard_gate_failures"]] == ["exact_unsafe_self_only"]

    proxy_flags = [f["flag"] for f in by_id["PROXY_p.V5I"]["review_flags"]]
    assert "dna_only_proxy_expression" in proxy_flags
    assert "high_hydrophobicity" in proxy_flags
    assert by_id["PROXY_p.V5I"]["default_recommendation"] == "Review"

    out_path = write_tumor_board_dossier(dossier, tmp_path / "tumor_board_dossier.json")
    reloaded = json.loads(open(out_path).read())
    assert reloaded["dossier_sha256"] == dossier["dossier_sha256"]


def test_validate_and_apply_overrides_human_ai_and_safety_enforcement(tmp_path):
    c1 = _make_candidate("TP53_p.R175H", "EVVRHCPHHERCSDSDGLAPPQHLI", {"HLA-B*07:02": 4.0}, n_alt_reads=30)
    c2 = _make_candidate("ARTIFACT_p.A1V", "ACDEFGHIKLMNPQRSTVWYACDEF", {"HLA-B*07:02": 3.0}, n_alt_reads=15)
    c3 = _make_candidate("KRAS_p.G12D", "MTEYKLVVVGADGVGKSALTIQLIQ", {"HLA-B*07:02": 2.5}, n_alt_reads=25)
    ranked = [c1, c2, c3]

    dossier = build_tumor_board_dossier(
        ranked,
        patient_id="PT-002",
        mhc_alleles=["HLA-B*07:02"],
        dna_vaf_by_variant={c2[0]: {"tumor_vaf": 0.30, "normal_alt_vaf": 0.12}},
    )
    sha = dossier["dossier_sha256"]

    # 1. Valid AI reviewer bundle rejects ARTIFACT and accepts TP53 & KRAS
    overrides_file = tmp_path / "ai_overrides.json"
    overrides_file.write_text(json.dumps({
        "schema": OVERRIDES_SCHEMA,
        "dossier_sha256": sha,
        "reviewer_id": "clinical-ai-v1",
        "reviewer_kind": "ai",
        "decisions": [
            {
                "candidate_id": "TP53_p.R175H",
                "decision": "Accept",
                "confidence": 0.96,
                "cited_evidence_keys": ["transcriptomic_integrity.n_rna_alt_reads"],
                "rationale": "Strong clonal expression and retained MHC restriction.",
            },
            {
                "candidate_id": "ARTIFACT_p.A1V",
                "decision": "Reject",
                "confidence": 0.99,
                "cited_evidence_keys": ["genomic_telemetry.normal_alt_vaf"],
                "rationale": "Normal alt VAF 12% exceeds 2% hard safety gate.",
            },
            {
                "candidate_id": "KRAS_p.G12D",
                "decision": "Review",
                "confidence": 0.78,
                "cited_evidence_keys": ["immunogenicity_and_geometry.target_epitope_count"],
                "rationale": "Retained for clinician tie-break.",
            },
        ],
    }))
    loaded = load_tumor_board_overrides(overrides_file)
    filtered, audit = validate_and_apply_overrides(ranked, dossier, loaded)
    assert [s.short_description for s, _ in filtered] == ["TP53_p.R175H", "KRAS_p.G12D"]
    assert audit["reviewer_id"] == "clinical-ai-v1"
    assert audit["reviewer_kind"] == "ai"
    assert audit["rejected_ids"] == ["ARTIFACT_p.A1V"]
    assert audit["accepted_ids"] == ["TP53_p.R175H"]
    assert audit["review_ids"] == ["KRAS_p.G12D"]

    # 2. Hard safety gate CANNOT be bypassed to Accept
    bypass_attempt = TumorBoardOverrides(
        dossier_sha256=sha,
        reviewer_id="rogue-reviewer",
        reviewer_kind="ai",
        decisions=(
            TumorBoardDecision(
                candidate_id="ARTIFACT_p.A1V",
                decision="Accept",
                confidence=0.9,
                cited_evidence_keys=("genomic_telemetry.normal_alt_vaf",),
                rationale="Trying to accept normal leakage candidate.",
            ),
        ),
    )
    with pytest.raises(TumorBoardValidationError, match="hard gate"):
        validate_and_apply_overrides(ranked, dossier, bypass_attempt)

    # 3. Stale dossier_sha256 rejected
    stale = TumorBoardOverrides(
        dossier_sha256="0" * 64,
        reviewer_id="human-1",
        decisions=(
            TumorBoardDecision(
                candidate_id="TP53_p.R175H",
                decision="Accept",
                cited_evidence_keys=("transcriptomic_integrity.n_rna_alt_reads",),
                rationale="ok",
            ),
        ),
    )
    with pytest.raises(TumorBoardValidationError, match="dossier_sha256 mismatch"):
        validate_and_apply_overrides(ranked, dossier, stale)

    # 4. Invalid citation key rejected
    bad_citation = TumorBoardOverrides(
        dossier_sha256=sha,
        reviewer_id="ai-1",
        reviewer_kind="ai",
        decisions=(
            TumorBoardDecision(
                candidate_id="TP53_p.R175H",
                decision="Accept",
                cited_evidence_keys=("nonexistent.hallucinated_key",),
                rationale="ok",
            ),
        ),
    )
    with pytest.raises(TumorBoardValidationError, match="unknown evidence key"):
        validate_and_apply_overrides(ranked, dossier, bad_citation)


def test_interactive_html_and_audit_templates_render():
    c1 = _make_candidate("TP53_p.R175H", "EVVRHCPHHERCSDSDGLAPPQHLI", {"HLA-B*07:02": 4.0}, n_alt_reads=30)
    c2 = _make_candidate("ARTIFACT_p.A1V", "ACDEFGHIKLMNPQRSTVWYACDEF", {"HLA-B*07:02": 3.0}, n_alt_reads=15)
    ranked = [c1, c2]
    dossier = build_tumor_board_dossier(
        ranked,
        patient_id="PT-003",
        mhc_alleles=["HLA-B*07:02"],
        dna_vaf_by_variant={c2[0]: {"normal_alt_vaf": 0.15}},
    )
    overrides = TumorBoardOverrides(
        dossier_sha256=dossier["dossier_sha256"],
        reviewer_id="dr-smith",
        reviewer_kind="human",
        decisions=(
            TumorBoardDecision(
                candidate_id="ARTIFACT_p.A1V",
                decision="Reject",
                confidence=1.0,
                cited_evidence_keys=("genomic_telemetry.normal_alt_vaf",),
                rationale="High normal VAF.",
            ),
        ),
    )
    _, audit = validate_and_apply_overrides(ranked, dossier, overrides)

    env = Environment(loader=PackageLoader("vaxrank", "templates"))
    interactive = env.get_template("interactive_review.html").render(
        tumor_board_dossier=dossier,
    )
    assert "vaxrank-tumor-board-dossier" in interactive
    assert "Download tumor_board_overrides.json" in interactive
    assert "disabled title=\"Blocked by deterministic hard gate\"" in interactive
    assert dossier["dossier_sha256"] in interactive

    audit_html = env.get_template("tumor_board_audit.html").render(
        tumor_board_overrides_audit=audit,
    )
    assert "dr-smith" in audit_html
    assert "ARTIFACT_p.A1V" in audit_html

    audit_txt = env.get_template("tumor_board_audit.txt").render(
        tumor_board_overrides_audit=audit,
    )
    assert "dr-smith (human)" in audit_txt
    assert "ARTIFACT_p.A1V: Reject" in audit_txt


def test_cli_flags_for_tumor_board_dossier_and_overrides():
    args_direct = parse_vaxrank_args([
        "--vcf", "a.vcf",
        "--bam", "a.bam",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--output-tumor-board-dossier", "dossier.json",
        "--interactive-review-html",
        "--selection-overrides", "overrides.json",
    ])
    assert args_direct.output_tumor_board_dossier == "dossier.json"
    assert args_direct.interactive_review_html is True
    assert args_direct.selection_overrides == "overrides.json"

    args_ext = parse_vaxrank_args([
        "--input-epitopes", "candidate_predictions.tsv",
        "--selection-overrides", "overrides.json",
        "--output-dir", "out/",
    ])
    assert args_ext.selection_overrides == "overrides.json"


def test_hard_gate_blocks_review_and_omission_and_disables_review_in_html():
    c1 = _make_candidate("TP53_p.R175H", "EVVRHCPHHERCSDSDGLAPPQHLI", {"HLA-B*07:02": 4.0}, n_alt_reads=30)
    c2 = _make_candidate("ARTIFACT_p.A1V", "ACDEFGHIKLMNPQRSTVWYACDEF", {"HLA-B*07:02": 3.0}, n_alt_reads=15)
    ranked = [c1, c2]
    dossier = build_tumor_board_dossier(
        ranked,
        patient_id="PT-HARD",
        mhc_alleles=["HLA-B*07:02"],
        dna_vaf_by_variant={c2[0]: {"normal_alt_vaf": 0.15}},
    )
    sha = dossier["dossier_sha256"]

    # Attempting to override a hard-gate failure to Review MUST fail
    bad_review = TumorBoardOverrides(
        dossier_sha256=sha,
        reviewer_id="rev-1",
        reviewer_kind="human",
        decisions=(
            TumorBoardDecision(
                candidate_id="ARTIFACT_p.A1V",
                decision="Review",
                cited_evidence_keys=("genomic_telemetry.normal_alt_vaf",),
                rationale="Try to sneak past hard gate via Review.",
            ),
        ),
    )
    with pytest.raises(TumorBoardValidationError, match="Only 'Reject' is permitted"):
        validate_and_apply_overrides(ranked, dossier, bad_review)

    # Omitting the hard-gate-failed candidate MUST STILL exclude it!
    omitted_overrides = TumorBoardOverrides(
        dossier_sha256=sha,
        reviewer_id="rev-1",
        reviewer_kind="human",
        decisions=(
            TumorBoardDecision(
                candidate_id="TP53_p.R175H",
                decision="Accept",
                cited_evidence_keys=("immunogenicity_and_geometry.best_percentile_rank",),
                rationale="Clean candidate.",
            ),
        ),
    )
    filtered, audit = validate_and_apply_overrides(ranked, dossier, omitted_overrides)
    assert [s.short_description for s, _ in filtered] == ["TP53_p.R175H"]
    assert "ARTIFACT_p.A1V" in audit["rejected_ids"]
    assert audit["auto_rejected_hard_gate_ids"] == ["ARTIFACT_p.A1V"]

    # Interactive HTML must disable BOTH Accept and Review for hard-gate rows
    env = Environment(loader=PackageLoader("vaxrank", "templates"))
    html = env.get_template("interactive_review.html").render(tumor_board_dossier=dossier)
    assert html.count('disabled title="Blocked by deterministic hard gate"') >= 2
    assert "Review (blocked)" in html


def test_all_restricting_mhc_lost_ignores_non_binder_predictions():
    from mhctools.pred import Prediction
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.vaccine_peptide import VaccinePeptide

    aa = "WYRSTVQPMLKIHGFEDCAWQRSTV"
    frag = SimpleNamespace(
        amino_acids=aa,
        gene_name="LOHONLY",
        mutant_amino_acid_start_offset=4,
        mutant_amino_acid_end_offset=5,
        n_alt_reads=20,
        n_ref_reads=5,
        n_overlapping_reads=25,
        n_alt_reads_supporting_protein_sequence=20,
        n_other_reads=0,
        n_mutant_amino_acids=1,
        variant=SimpleNamespace(short_description="LOHONLY_p.L2P"),
    )
    preds = (
        Prediction(
            peptide=aa[:9], allele="HLA-A*02:01", kind="pMHC_affinity",
            value=45.0, percentile_rank=0.4, score=4.0, source_sequence_name="LOHONLY",
        ),
        # 40,000 nM non-binder on retained allele must NOT suppress MHC-LOH hard gate!
        Prediction(
            peptide=aa[:9], allele="HLA-B*07:02", kind="pMHC_affinity",
            value=40000.0, percentile_rank=65.0, score=0.002, source_sequence_name="LOHONLY",
        ),
    )
    ep = CandidateEpitope(
        sequence=aa[:9],
        source_sequence=aa,
        offset=0,
        overlaps_mutation=True,
        overlaps_targetable=True,
        occurs_in_reference=False,
        per_allele_scores={"HLA-A*02:01": 4.0},
        predictions=preds,
    )
    vp = VaccinePeptide(mutant_protein_fragment=frag, epitopes=[ep])
    ranked = [(frag.variant, [vp])]
    bundle = MHCAlleleStatusBundle(
        alleles=(
            MHCAlleleStatus(allele="HLA-A*02:01", dna_loh=True, dna_copy_number=0.0),
            MHCAlleleStatus(allele="HLA-B*07:02", status="retained"),
        )
    )
    ranked, mhc_audit = apply_mhc_allele_status_to_ranked(ranked, bundle)
    dossier = build_tumor_board_dossier(
        ranked,
        patient_id="PT-NONBINDER",
        mhc_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        mhc_allele_status=mhc_audit,
    )
    cand = dossier["candidates"][0]
    assert cand["tumor_evolution_and_mhc"]["all_restricting_mhc_lost_or_silenced"] is True
    assert any(g["gate"] == "all_restricting_mhc_lost_or_silenced" for g in cand["hard_gate_failures"])


def test_zero_vaf_preserved_and_dossier_sha256_invariant_to_products():
    from vaxrank.tumor_board import overrides_from_dict

    c1 = _make_candidate("TP53_p.R175H", "EVVRHCPHHERCSDSDGLAPPQHLI", {"HLA-B*07:02": 4.0}, n_alt_reads=30)
    ranked = [c1]
    pre_assembly = build_tumor_board_dossier(
        ranked,
        patient_id="PT-ZERO",
        mhc_alleles=["HLA-B*07:02"],
        dna_vaf_by_variant={c1[0]: {"tumor_vaf": 0.0, "normal_alt_vaf": 0.0}},
    )
    assert pre_assembly["candidates"][0]["genomic_telemetry"]["tumor_vaf"] == 0.0
    assert pre_assembly["candidates"][0]["genomic_telemetry"]["normal_alt_vaf"] == 0.0

    # Simulate post-assembly product & pool attachment
    fake_construct = SimpleNamespace(
        name="P1",
        sequence="EVVRHCPHHERCSDSDGLAPPQHLI",
        components={"pool": 1},
    )
    post_assembly = build_tumor_board_dossier(
        ranked,
        patient_id="PT-ZERO",
        mhc_alleles=["HLA-B*07:02"],
        dna_vaf_by_variant={c1[0]: {"tumor_vaf": 0.0, "normal_alt_vaf": 0.0}},
        products={"peptide": [fake_construct]},
    )
    assert post_assembly["candidates"][0]["construct_and_manufacturability"]["emitted_products"] == ["P1"]
    assert post_assembly["candidates"][0]["construct_and_manufacturability"]["peptide_pool"] == 1
    assert pre_assembly["dossier_sha256"] == post_assembly["dossier_sha256"]

    # Round-trip overrides_from_dict
    ov = TumorBoardOverrides(
        dossier_sha256=post_assembly["dossier_sha256"],
        reviewer_id="rev-rt",
        reviewer_kind="human",
        decisions=(
            TumorBoardDecision(
                candidate_id="TP53_p.R175H",
                decision="Accept",
                cited_evidence_keys=("genomic_telemetry.normal_alt_vaf",),
                rationale="Zero normal VAF.",
            ),
        ),
    )
    restored = overrides_from_dict(ov.to_dict())
    assert restored == ov


def test_entry_point_flag_registration_for_mhc_status_and_tumor_board():
    from vaxrank.cli.entry_point import _ARG_GROUPS, _INPUT_FILE_FLAGS, _OUTPUT_FILE_FLAGS

    assert "mhc_allele_status" in _INPUT_FILE_FLAGS
    assert "selection_overrides" in _INPUT_FILE_FLAGS
    assert "output_tumor_board_dossier" in _OUTPUT_FILE_FLAGS
    all_grouped = {k for _, keys in _ARG_GROUPS for k in keys}
    for flag in (
        "mhc_allele_status",
        "mhc_loh_alleles",
        "selection_overrides",
        "output_tumor_board_dossier",
        "interactive_review_html",
    ):
        assert flag in all_grouped
