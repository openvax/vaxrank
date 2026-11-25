"""Regression tests for core pipeline edge cases & cross-species MHC naming (#934)."""

from types import SimpleNamespace
import math

import pandas as pd
import pytest
from mhctools.pred import Prediction

from vaxrank.allele_evidence import AllelePolicy, apply_attribution, attribute_frame
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.coverage import _min
from vaxrank.epitope_io import normalize_hla_allele, normalize_mhc_allele
from vaxrank.gene_pathway_check import GenePathwayCheck
from vaxrank.patient_background import _project_pair
from vaxrank.report import _format_report_allele


def test_coverage_min_ignores_nan_and_none():
    assert _min(float("nan"), 0.5) == 0.5
    assert _min(0.5, float("nan")) == 0.5
    assert _min(float("nan"), float("nan")) is None
    assert _min(None, float("nan")) is None
    assert _min(0.2, 0.8) == 0.2


def test_allele_evidence_treats_nan_allele_as_allele_free():
    df = pd.DataFrame([
        {
            "source_sequence_name": "S1",
            "peptide": "SIINFEKL",
            "peptide_offset": 0,
            "allele": "HLA-A*02:01",
            "kind": "pMHC_affinity",
            "prediction_method_name": "netmhcpan",
            "value": 50.0,
            "percentile_rank": 0.4,
            "score": 0.9,
        },
        {
            "source_sequence_name": "S1",
            "peptide": "SIINFEKL",
            "peptide_offset": 0,
            "allele": float("nan"),
            "kind": "proteasome_cleavage",
            "prediction_method_name": "pepsickle",
            "value": 0.8,
            "percentile_rank": float("nan"),
            "score": 0.8,
        },
    ])
    policy = AllelePolicy.parse("top:1")
    attributions = attribute_frame(
        df, policy, group_columns=["source_sequence_name", "peptide", "peptide_offset"]
    )
    assert ("S1", "SIINFEKL", 0) in attributions
    assert attributions[("S1", "SIINFEKL", 0)][0].allele == "HLA-A*02:01"
    rewritten = apply_attribution(
        df, attributions, policy,
        group_columns=["source_sequence_name", "peptide", "peptide_offset"],
    )
    cleavage_rows = rewritten[rewritten["kind"] == "proteasome_cleavage"]
    assert list(cleavage_rows["allele"]) == ["HLA-A*02:01"]


def test_peptide_best_skips_none_and_nan_scores():
    ep = CandidateEpitope(
        sequence="SIINFEKL",
        source_sequence="AASIINFEKLAA",
        offset=2,
        predictions=(
            Prediction(
                peptide="SIINFEKL", allele="HLA-A*02:01", kind="pMHC_affinity",
                value=500.0, score=None, percentile_rank=1.0,
            ),
            Prediction(
                peptide="SIINFEKL", allele="HLA-B*07:02", kind="pMHC_affinity",
                value=100.0, score=0.85, percentile_rank=0.3,
            ),
        ),
    )
    best = ep.best("pMHC_affinity")
    assert best is not None
    assert best.allele == "HLA-B*07:02"
    assert math.isclose(best.score, 0.85)

    ep_all_none = CandidateEpitope(
        sequence="SIINFEKL",
        source_sequence="AASIINFEKLAA",
        offset=2,
        predictions=(
            Prediction(
                peptide="SIINFEKL", allele="HLA-A*02:01", kind="pMHC_affinity",
                value=500.0, score=None, percentile_rank=1.0,
            ),
        ),
    )
    assert ep_all_none.best("pMHC_affinity") is None


def test_patient_background_internal_indel_uses_half_open_right_boundary():
    pair = {
        "pair_id": "p1",
        "translation_id": "t1",
        "transcript_id": "ENST00001",
        "reason": None,
        "mutant": {"coding_interval": [0, 30]},
        "baseline": {"coding_interval": [0, 27], "protein_sequence": "ACDEFGHIK"},
        # 3-nt insertion right at observed_interval [15, 18] (after peptide [0, 5) -> nt [0, 15))
        "removed_edits": [{"observed_interval": [15, 18], "reference_bases": ""}],
    }
    res = _project_pair(pair, 0, 5)
    assert res["state"] == "available"
    assert res["sequence"] == "ACDEF"
    assert res["position_aligned"] is True


def test_gene_pathway_check_handles_none_effect():
    checker = GenePathwayCheck()

    class _EmptyEffects:
        def top_priority_effect(self):
            return None

    fake_variant = SimpleNamespace(
        effects=lambda: _EmptyEffects(),
        gene_ids=[],
    )
    result = checker.make_variant_dict(fake_variant)
    assert list(result.values()) == [False, False, False, False]


def test_cross_species_report_allele_formatting_and_normalize_mhc_allele():
    assert _format_report_allele("HLA-A*02:01") == "A*02:01"
    assert _format_report_allele("DLA-88*034:01") == "DLA-88*034:01"
    assert _format_report_allele("H2-Kb") == "H2-K*b"
    assert _format_report_allele("SLA-1*04:01") == "SLA-1*04:01"
    assert normalize_mhc_allele("HLA-A02:01") == normalize_hla_allele("HLA-A02:01") == "HLA-A*02:01"
    assert normalize_mhc_allele("H2-Kb") == "H2-K*b"
    assert normalize_mhc_allele("DLA-88*034:01") == "DLA-88*034:01"


@pytest.mark.parametrize("bad_doc_token", ["G4S3", "vaxrank-mutation-figure"])
def test_docs_do_not_reintroduce_fixed_contradictions(bad_doc_token):
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    assert bad_doc_token not in (root / "docs" / "vaccine-designs.md").read_text()
    assert bad_doc_token not in (root / "README.md").read_text()
