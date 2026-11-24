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

"""Tests for cross-species MHCAlleleStatus contract and policy masking
(openvax/vaxrank#921)."""

import json
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
    MHCAlleleStatusPolicy,
    apply_mhc_allele_status_to_ranked,
    bundle_from_loh_alleles,
    load_mhc_allele_status,
    merge_mhc_allele_status_bundles,
    normalize_mhc_allele_name,
    resolve_effective_status,
)
from vaxrank.vaccine_peptide import VaccinePeptide


def _make_ranked_item(name, seq, allele_scores, *, expression_reads=25):
    preds = tuple(
        Prediction(
            kind="pMHC_affinity",
            peptide=seq[:9],
            allele=allele,
            score=score,
            percentile_rank=0.4 if score > 0 else 50.0,
        )
        for allele, score in allele_scores.items()
    )
    ep = CandidateEpitope(
        sequence=seq[:9],
        source_sequence=seq,
        offset=0,
        predictions=preds,
        overlaps_mutation=True,
        overlaps_targetable=True,
        per_allele_scores=dict(allele_scores),
    )
    fragment = SimpleNamespace(
        amino_acids=seq,
        gene_name=name.split("_")[0],
        mutant_amino_acid_start_offset=2,
        mutant_amino_acid_end_offset=3,
        n_alt_reads=expression_reads,
        n_ref_reads=10,
        n_overlapping_reads=expression_reads + 10,
        n_alt_reads_supporting_protein_sequence=expression_reads,
        supporting_reference_transcripts=(),
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=[ep],
        combined_score_expr="sqrt(n_alt_reads) * target_epitope_score",
    )
    source = SimpleNamespace(short_description=name, gene_names=[name.split("_")[0]])
    return (source, [vp])


def test_policy_validation_and_yaml_loading(tmp_path):
    pol = MHCAlleleStatusPolicy()
    assert pol.mask_lost_alleles is True
    assert pol.mask_silenced_alleles is True
    assert pol.mask_unresolvable_alleles is False
    assert pol.expression_weight_mode == "binary"
    assert pol.min_scrna_resolvable_umis == 20

    with pytest.raises(ValueError, match="expression_weight_mode"):
        MHCAlleleStatusPolicy(expression_weight_mode="invalid")
    with pytest.raises(ValueError, match="min_scrna_resolvable_umis"):
        MHCAlleleStatusPolicy(min_scrna_resolvable_umis=0)
    with pytest.raises(ValueError, match="Unknown mhc_allele_status"):
        MHCAlleleStatusPolicy.resolve({"bogus": 1})

    cfg_path = tmp_path / "mhc_policy.yaml"
    cfg_path.write_text(
        "mhc_allele_status:\n"
        "  mask_lost_alleles: true\n"
        "  mask_silenced_alleles: true\n"
        "  mask_unresolvable_alleles: false\n"
        "  expression_weight_mode: continuous_ratio\n"
        "  min_scrna_resolvable_umis: 30\n"
        "  silenced_tumor_to_immune_ratio: 0.15\n"
    )
    loaded = load_vaxrank_config(config_path=str(cfg_path))
    resolved = MHCAlleleStatusPolicy.resolve(loaded["mhc_allele_status"])
    assert resolved.expression_weight_mode == "continuous_ratio"
    assert resolved.min_scrna_resolvable_umis == 30
    assert resolved.silenced_tumor_to_immune_ratio == 0.15


def test_cross_species_mhcgnomes_normalization():
    assert normalize_mhc_allele_name("A*02:01") == "HLA-A*02:01"
    assert normalize_mhc_allele_name("HLA-A02:01") == "HLA-A*02:01"
    assert normalize_mhc_allele_name("DLA-88*034:01") == "DLA-88*034:01"
    assert normalize_mhc_allele_name("H2-Kb") == "H2-K*b"
    assert normalize_mhc_allele_name("H-2-Kb") == "H2-K*b"


def test_multi_assay_resolution_dna_bulk_and_scrna():
    pol = MHCAlleleStatusPolicy()

    # 1. DNA copy-number LOH -> lost
    rec_loh = MHCAlleleStatus(
        allele="HLA-A*02:01", dna_copy_number=0.1, dna_loh=True,
        dna_loh_pvalue=1e-6, assay_source="LOHHLA",
    )
    dec_loh = resolve_effective_status(rec_loh, pol)
    assert dec_loh.effective_status == "lost"
    assert dec_loh.badge == "MHC-LOH"
    assert dec_loh.masked is True
    assert dec_loh.ranking_weight == 0.0

    # 2. 10x 5' scRNA tumor-cluster silencing -> silenced
    rec_silenced = MHCAlleleStatus(
        allele="HLA-B*07:02",
        scrna_tumor_cluster_ntpm=8.0,
        scrna_immune_control_ntpm=200.0,
        scrna_resolvable_umis=95,
        scrna_chemistry="10x_5p",
        assay_source="scHLAcount",
    )
    dec_sil = resolve_effective_status(rec_silenced, pol)
    assert dec_sil.effective_status == "silenced"
    assert dec_sil.badge == "MHC-Silenced"
    assert dec_sil.masked is True
    assert dec_sil.ranking_weight == 0.0

    # 3. 10x 3' scRNA low UMIs -> unresolvable (NEVER lost/silenced!)
    rec_3p = MHCAlleleStatus(
        allele="HLA-C*07:01",
        scrna_tumor_cluster_ntpm=0.0,
        scrna_immune_control_ntpm=150.0,
        scrna_resolvable_umis=5,
        scrna_chemistry="10x_3p",
    )
    dec_3p = resolve_effective_status(rec_3p, pol)
    assert dec_3p.effective_status == "unresolvable"
    assert dec_3p.badge == "MHC-Unresolvable"
    assert dec_3p.masked is False
    assert dec_3p.ranking_weight == 1.0
    assert "scrna_insufficient_resolvable_umis" in dec_3p.caveats
    assert "scrna_3p_exon2_3_coverage_gap" in dec_3p.caveats

    # 4. Bulk RNA ASE with low tumor purity -> purity-confounded caveat
    rec_low_purity_low_ase = MHCAlleleStatus(
        allele="HLA-A*24:02",
        bulk_rna_tpm=45.0,
        bulk_rna_allele_ratio=0.05,
        tumor_purity=0.25,
    )
    dec_lp = resolve_effective_status(rec_low_purity_low_ase, pol)
    assert dec_lp.effective_status == "unresolvable"
    assert dec_lp.masked is False
    assert "bulk_rna_purity_confounded" in dec_lp.caveats

    rec_high_purity_low_ase = MHCAlleleStatus(
        allele="HLA-A*24:02",
        bulk_rna_tpm=10.0,
        bulk_rna_allele_ratio=0.04,
        tumor_purity=0.80,
    )
    dec_hp = resolve_effective_status(rec_high_purity_low_ase, pol)
    assert dec_hp.effective_status == "silenced"
    assert dec_hp.masked is True

    # 5. Continuous ratio weighting mode
    pol_cont = MHCAlleleStatusPolicy(expression_weight_mode="continuous_ratio")
    rec_partial = MHCAlleleStatus(
        allele="HLA-B*44:02",
        scrna_tumor_cluster_ntpm=60.0,
        scrna_immune_control_ntpm=200.0,
        scrna_resolvable_umis=80,
        scrna_chemistry="10x_5p",
    )
    dec_part = resolve_effective_status(rec_partial, pol_cont)
    assert dec_part.effective_status == "retained"
    assert dec_part.masked is False
    assert dec_part.ranking_weight == pytest.approx(0.3)


def test_loaders_json_tsv_and_merge(tmp_path):
    json_path = tmp_path / "mhc_status.json"
    json_path.write_text(json.dumps({
        "alleles": [
            {"allele": "HLA-A*02:01", "dna_loh": True, "dna_copy_number": 0.0, "assay_source": "LOHHLA"},
            {"allele": "HLA-B*07:02", "status": "retained"},
        ],
        "apm": {
            "b2m_status": "truncating_mutation_loh",
            "tap1_status": "intact",
        },
    }))
    bundle_json = load_mhc_allele_status(json_path)
    assert len(bundle_json.alleles) == 2
    assert bundle_json.apm is not None
    assert bundle_json.apm.disrupted_genes() == ["B2M"]

    tsv_path = tmp_path / "dog_mhc.tsv"
    tsv_path.write_text(
        "allele\tstatus\tdna_copy_number\tscrna_resolvable_umis\tscrna_chemistry\n"
        "DLA-88*034:01\tlost\t0.0\t\t\n"
        "DLA-88*508:01\t\t1.9\t4\t10x_3p\n"
    )
    bundle_tsv = load_mhc_allele_status(tsv_path)
    assert [a.allele for a in bundle_tsv.alleles] == ["DLA-88*034:01", "DLA-88*508:01"]

    cli_override = bundle_from_loh_alleles("HLA-B*07:02")
    merged = merge_mhc_allele_status_bundles(bundle_json, cli_override)
    assert merged.by_allele()["HLA-B*07:02"].status == "lost"
    assert merged.apm.disrupted_genes() == ["B2M"]


def test_apply_mhc_allele_status_reranks_across_human_dog_and_mouse():
    # Canine example: GENE1 binds lost DLA-88*034:01 strongly (score 5.0);
    # GENE2 binds retained DLA-88*508:01 (score 2.0);
    # GENE3 binds unresolvable DLA-88*012:01 (score 1.5, keeps full credit).
    ranked = [
        _make_ranked_item("GENE1_p.A1V", "ACDEFGHIKLMNPQRSTVWY", {"DLA-88*034:01": 5.0}),
        _make_ranked_item("GENE2_p.B2L", "WYRSTVQPMLKIHGFEDCAW", {"DLA-88*508:01": 2.0}),
        _make_ranked_item("GENE3_p.C3K", "MNPQRSTVWYACDEFGHIKL", {"DLA-88*012:01": 1.5}),
    ]
    bundle = MHCAlleleStatusBundle(
        alleles=(
            MHCAlleleStatus(allele="DLA-88*034:01", dna_loh=True, dna_copy_number=0.0),
            MHCAlleleStatus(allele="DLA-88*508:01", status="retained"),
            MHCAlleleStatus(
                allele="DLA-88*012:01",
                scrna_tumor_cluster_ntpm=0.0,
                scrna_immune_control_ntpm=100.0,
                scrna_resolvable_umis=3,
                scrna_chemistry="10x_3p",
            ),
        ),
        apm=APMStatus(b2m_status="intact"),
    )

    updated, audit = apply_mhc_allele_status_to_ranked(ranked, bundle)
    assert [s.short_description for s, _ in updated] == [
        "GENE2_p.B2L", "GENE3_p.C3K", "GENE1_p.A1V"
    ]
    # GENE1 score zeroed, yet raw Prediction preserved!
    gene1_vp = updated[-1][1][0]
    assert gene1_vp.target_epitope_score == 0.0
    assert gene1_vp.combined_score == 0.0
    raw_preds = gene1_vp.epitopes[0].predictions_flat()
    assert len(raw_preds) == 1
    assert raw_preds[0].score == 5.0

    # Unresolvable allele still credited
    gene3_vp = updated[1][1][0]
    assert gene3_vp.target_epitope_score == pytest.approx(1.5)
    assert audit["masked_alleles"] == ["DLA-88*034:01"]
    assert audit["unresolvable_alleles"] == ["DLA-88*012:01"]
    assert audit["affected_epitope_count"] == 1
    assert audit["total_credit_removed"] == pytest.approx(5.0)

    # Mouse allele matching across spelling variants (H-2-Kb vs H2-Kb)
    mouse_ranked = [
        _make_ranked_item("M1_p.R1Q", "ACDEFGHIKLMNPQRSTVWY", {"H-2-Kb": 4.0}),
        _make_ranked_item("M2_p.K2E", "WYRSTVQPMLKIHGFEDCAW", {"H-2-Db": 1.0}),
    ]
    mouse_bundle = bundle_from_loh_alleles("H2-Kb")
    mouse_updated, mouse_audit = apply_mhc_allele_status_to_ranked(mouse_ranked, mouse_bundle)
    assert [s.short_description for s, _ in mouse_updated] == ["M2_p.K2E", "M1_p.R1Q"]
    assert mouse_audit["masked_alleles"] == ["H2-K*b"]


def test_mhc_allele_status_report_templates_render():
    ranked = [
        _make_ranked_item("TP53_p.R175H", "EVVRHCPHHERCSDSDGLAPPQHLI", {"HLA-A*02:01": 3.0}),
    ]
    bundle = MHCAlleleStatusBundle(
        alleles=(
            MHCAlleleStatus(allele="HLA-A*02:01", dna_loh=True, dna_copy_number=0.0, assay_source="LOHHLA"),
            MHCAlleleStatus(
                allele="HLA-B*07:02",
                scrna_tumor_cluster_ntpm=5.0,
                scrna_immune_control_ntpm=180.0,
                scrna_resolvable_umis=75,
                scrna_chemistry="10x_5p",
            ),
        ),
        apm=APMStatus(b2m_status="lof_mutated", tap1_status="intact"),
    )
    _, audit = apply_mhc_allele_status_to_ranked(ranked, bundle)

    env = Environment(loader=PackageLoader("vaxrank", "templates"))
    html_out = env.get_template("mhc_allele_status.html").render(mhc_allele_status=audit)
    assert "MHC-LOH" in html_out
    assert "MHC-Silenced" in html_out
    assert "APM disruption: B2M" in html_out

    txt_out = env.get_template("mhc_allele_status.txt").render(mhc_allele_status=audit)
    assert "HLA-A*02:01: MHC-LOH" in txt_out
    assert "HLA-B*07:02: MHC-Silenced" in txt_out
    assert "DISRUPTED: B2M" in txt_out


def test_cli_flags_for_mhc_allele_status():
    args_direct = parse_vaxrank_args([
        "--vcf", "a.vcf",
        "--bam", "a.bam",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01,HLA-B*07:02",
        "--mhc-allele-status", "status.json",
        "--mhc-loh-alleles", "HLA-A*02:01",
    ])
    assert args_direct.mhc_allele_status == "status.json"
    assert args_direct.mhc_loh_alleles == "HLA-A*02:01"

    args_ext = parse_vaxrank_args([
        "--input-pvacseq", "p.tsv",
        "--mhc-loh-alleles", "DLA-88*034:01",
    ])
    assert args_ext.mhc_loh_alleles == "DLA-88*034:01"
