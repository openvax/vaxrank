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
    assert "Loss of APM: B2M" in html_out
    assert "DNA Status" in html_out
    assert "RNA Status" in html_out

    txt_out = env.get_template("mhc_allele_status.txt").render(mhc_allele_status=audit)
    assert "HLA-A*02:01: MHC-LOH" in txt_out
    assert "HLA-B*07:02: MHC-Silenced" in txt_out
    assert "LOSS OF APM: B2M" in txt_out


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


def test_apm_negation_words_not_falsely_disrupted():
    from vaxrank.mhc_allele_status import APMGeneStatus

    clean_apm = APMStatus(
        b2m_status="no_loh",
        tap1_status="not_deleted",
        tap2_status="unmutated",
        nlrc5_status="non-truncating",
        jak1_status="wildtype",
        jak2_status="intact",
        genes=(
            APMGeneStatus(gene="CALR", dna_status="not_mutated", rna_status="normal"),
            APMGeneStatus(gene="ERAP1", dna_status="intact", rna_status="retained"),
        ),
    )
    assert clean_apm.disrupted_genes() == []

    disrupted_apm = APMStatus(
        b2m_status="frameshift_lof",
        tap1_status="intact",
        genes=(
            APMGeneStatus(
                gene="NLRC5",
                dna_status="intact",
                rna_status="silenced",
                rna_ntpm=0.4,
                immune_control_ntpm=45.0,
            ),
        ),
    )
    assert disrupted_apm.disrupted_genes() == ["B2M", "NLRC5"]


def test_low_umi_scrna_does_not_short_circuit_high_purity_bulk_ase():
    rec = MHCAlleleStatus(
        allele="HLA-A*02:01",
        scrna_tumor_cluster_ntpm=0.0,
        scrna_immune_control_ntpm=100.0,
        scrna_resolvable_umis=2,
        scrna_chemistry="10x_3p",
        bulk_rna_allele_ratio=0.02,
        tumor_purity=0.80,
    )
    dec = resolve_effective_status(rec)
    assert dec.effective_status == "silenced"
    assert dec.rna_status == "silenced"
    assert dec.dna_status == "unknown"
    assert dec.loss_mechanism == "rna_silenced"
    assert dec.masked is True
    assert "scrna_insufficient_resolvable_umis" in dec.caveats


def test_dna_somatic_lof_vs_missense_and_separate_dna_rna_axes():
    lof_rec = MHCAlleleStatus(
        allele="DLA-88*034:01",
        dna_somatic_variants=({"consequence": "frameshift", "hgvs_p": "p.R102fs"},),
        bulk_rna_allele_ratio=0.45,
        tumor_purity=0.70,
    )
    lof_dec = resolve_effective_status(lof_rec)
    assert lof_dec.dna_status == "truncating_lof"
    assert lof_dec.rna_status == "expressed"
    assert lof_dec.effective_status == "lost"
    assert lof_dec.loss_mechanism == "dna_truncating_lof"
    assert lof_dec.masked is True

    missense_rec = MHCAlleleStatus(
        allele="H2-Kb",
        dna_somatic_variants=({"consequence": "missense", "hgvs_p": "p.V55A"},),
        bulk_rna_allele_ratio=0.48,
        tumor_purity=0.75,
    )
    mis_dec = resolve_effective_status(missense_rec)
    assert mis_dec.dna_status == "missense_or_inframe"
    assert mis_dec.rna_status == "expressed"
    assert mis_dec.effective_status == "retained"
    assert "dna_somatic_missense_variant" in mis_dec.caveats
    assert mis_dec.masked is False


def test_masking_preserves_self_epitopes_and_strips_ms_increment():
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.vaccine_peptide import VaccinePeptide

    aa = "EVVRHCPHHERCSDSDGLAPPQHLI"
    frag = SimpleNamespace(
        amino_acids=aa,
        gene_name="TP53",
        mutant_amino_acid_start_offset=4,
        mutant_amino_acid_end_offset=5,
        n_alt_reads=20,
        n_ref_reads=10,
        n_overlapping_reads=30,
        n_alt_reads_supporting_protein_sequence=20,
        n_other_reads=0,
        n_mutant_amino_acids=1,
        variant=SimpleNamespace(short_description="TP53_p.R175H"),
    )
    target_preds = (
        Prediction(
            peptide=aa[:9],
            allele="HLA-A*02:01",
            kind="pMHC_affinity",
            value=50.0,
            percentile_rank=0.4,
            score=1.0,
            source_sequence_name="TP53",
        ),
    )
    target_ep = CandidateEpitope(
        sequence=aa[:9],
        source_sequence=aa,
        offset=0,
        overlaps_mutation=True,
        overlaps_targetable=True,
        occurs_in_reference=False,
        per_allele_scores={"HLA-A*02:01": 1.0},
        ms_scoring={"HLA-A*02:01": {"baseline": 1.0, "adjusted": 2.5}},
        predictions=target_preds,
    )
    self_preds = (
        Prediction(
            peptide=aa[10:19],
            allele="HLA-A*02:01",
            kind="pMHC_affinity",
            value=40.0,
            percentile_rank=0.3,
            score=3.0,
            source_sequence_name="TP53",
        ),
    )
    self_ep = CandidateEpitope(
        sequence=aa[10:19],
        source_sequence=aa,
        offset=10,
        overlaps_mutation=False,
        overlaps_targetable=False,
        occurs_in_reference=True,
        per_allele_scores={"HLA-A*02:01": 3.0},
        predictions=self_preds,
    )
    vp = VaccinePeptide(mutant_protein_fragment=frag, epitopes=[target_ep, self_ep])
    assert vp.target_epitope_score == pytest.approx(2.5)
    assert vp.self_epitope_score == pytest.approx(3.0)

    bundle = bundle_from_loh_alleles("HLA-A*02:01")
    updated, audit = apply_mhc_allele_status_to_ranked([(frag.variant, [vp])], bundle)
    new_vp = updated[0][1][0]
    assert new_vp.target_epitope_score == pytest.approx(0.0)
    assert new_vp.target_epitopes[0].ms_increment == pytest.approx(0.0)
    assert "HLA-A*02:01" not in new_vp.target_epitopes[0].per_allele_scores
    # Normal-tissue self-epitope score MUST remain completely unmasked
    assert new_vp.self_epitope_score == pytest.approx(3.0)
    assert audit["loss_of_mhc"]["has_loss"] is True
    assert audit["loss_of_mhc"]["dna_lost_alleles"] == ["HLA-A*02:01"]


def test_resolve_target_alleles_excludes_masked_mhc_alleles():
    from vaxrank.cli.entry_point import resolve_target_alleles

    args = SimpleNamespace(
        vcf="a.vcf",
        _declared_mhc_alleles=["HLA-A*02:01", "HLA-B*07:02", "DLA-88*034:01"],
        _masked_mhc_alleles=["A0201", "DLA-88*034:01"],
    )
    assert resolve_target_alleles(args) == ["HLA-B*07:02"]


def test_bundle_from_dict_roundtrip_preserves_variants_and_apm_genes():
    from vaxrank.mhc_allele_status import APMGeneStatus, bundle_from_dict

    bundle = MHCAlleleStatusBundle(
        alleles=(
            MHCAlleleStatus(
                allele="HLA-A*02:01",
                dna_loh=True,
                dna_somatic_variants=({"gene": "HLA-A", "consequence": "nonsense"},),
            ),
        ),
        apm=APMStatus(
            b2m_status="lof_mutated",
            genes=(
                APMGeneStatus(
                    gene="B2M",
                    dna_status="truncating_lof",
                    rna_status="silenced",
                    dna_variants=({"consequence": "frameshift"},),
                    rna_ntpm=1.2,
                ),
            ),
        ),
        source_path="in_memory.json",
    )
    restored = bundle_from_dict(bundle.to_dict())
    assert len(restored.alleles) == 1
    assert restored.alleles[0].dna_somatic_variants == ({"gene": "HLA-A", "consequence": "nonsense"},)
    assert restored.apm is not None
    assert len(restored.apm.genes) == 1
    assert restored.apm.genes[0].gene == "B2M"
    assert restored.apm.genes[0].dna_status == "truncating_lof"
    assert restored.apm.genes[0].rna_ntpm == pytest.approx(1.2)
