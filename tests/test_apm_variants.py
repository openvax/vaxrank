"""Tests for cross-species DNA somatic MHC & APM variant/LoF detector (#935)."""

from types import SimpleNamespace

from jinja2 import Environment, PackageLoader
import pytest

from vaxrank.apm_variants import (
    APMVariantPolicy,
    classify_gene_role,
    classify_variant_effect_impact,
    detect_dna_mhc_and_apm_variants,
)
from vaxrank.cli.arg_parser import parse_vaxrank_args
from vaxrank.config import load_vaxrank_config
from vaxrank.mhc_allele_status import (
    APMGeneStatus,
    APMStatus,
    MHCAlleleStatusBundle,
    apply_mhc_allele_status_to_ranked,
    merge_mhc_allele_status_bundles,
    resolve_effective_status,
)
from vaxrank.tumor_board import build_tumor_board_dossier


def _fake_effect(cls_name: str, gene_name: str, short_desc: str):
    effect_cls = type(cls_name, (), {})
    eff = effect_cls()
    eff.gene_name = gene_name
    eff.short_description = short_desc
    return eff


def _fake_variant(desc: str, gene_name: str, effect_cls: str, effect_desc: str, *, species: str = "homo_sapiens"):
    eff = _fake_effect(effect_cls, gene_name, effect_desc)

    class _Effects:
        def top_priority_effect(self):
            return eff

    return SimpleNamespace(
        short_description=desc,
        gene_names=(gene_name,),
        genome=SimpleNamespace(species=species),
        effects=lambda: _Effects(),
    )


def test_cross_species_gene_role_classification():
    # Human, mouse, canine, swine, and dynamic mhcgnomes fallback
    assert classify_gene_role("B2M") == ("B2M", "apm", None)
    assert classify_gene_role("B2m") == ("B2M", "apm", None)
    assert classify_gene_role("Tapbp") == ("TAPBP", "apm", None)
    assert classify_gene_role("Psmb8") == ("PSMB8", "apm", None)
    assert classify_gene_role("Jak1") == ("JAK1", "ifng", None)
    assert classify_gene_role("Ifngr2") == ("IFNGR2", "ifng", None)
    assert classify_gene_role("Stat1") == ("STAT1", "ifng", None)
    assert classify_gene_role("HLA-A") == ("HLA-A", "mhc_class_i", ("HLA", "A"))
    assert classify_gene_role("DLA-88") == ("DLA-88", "mhc_class_i", ("DLA", "88"))
    assert classify_gene_role("H2-K1") == ("H2-K", "mhc_class_i", ("H2", "K"))
    assert classify_gene_role("H2-D1") == ("H2-D", "mhc_class_i", ("H2", "D"))
    assert classify_gene_role("SLA-1") == ("SLA-1", "mhc_class_i", ("SLA", "1"))
    assert classify_gene_role("TP53") == (None, None, None)

    # Extra trial-configured gene
    pol = APMVariantPolicy(extra_apm_genes=("CUSTOMAPM",))
    assert classify_gene_role("customapm", pol) == ("CUSTOMAPM", "apm", None)


def test_variant_effect_impact_classification():
    assert classify_variant_effect_impact(_fake_effect("FrameShift", "B2M", "p.S14fs"))[0] == "truncating_lof"
    assert classify_variant_effect_impact(_fake_effect("PrematureStop", "B2M", "p.Q22*"))[0] == "truncating_lof"
    assert classify_variant_effect_impact(_fake_effect("StartLoss", "B2M", "p.M1?"))[0] == "truncating_lof"
    assert classify_variant_effect_impact(_fake_effect("SpliceAcceptor", "NLRC5", "splice_acceptor"))[0] == "truncating_lof"
    assert classify_variant_effect_impact(_fake_effect("Substitution", "TAP1", "p.R412Q"))[0] == "missense_or_inframe"
    assert classify_variant_effect_impact(_fake_effect("Deletion", "ERAP1", "p.K100del"))[0] == "missense_or_inframe"
    assert classify_variant_effect_impact(_fake_effect("Silent", "JAK1", "p.L10L"))[0] == "silent_or_noncoding"
    assert classify_variant_effect_impact(_fake_effect("Intronic", "JAK2", "intronic"))[0] == "silent_or_noncoding"


def test_human_apm_detection_separates_lof_missense_silent_and_rna_confirmation():
    v_b2m = _fake_variant("15:44711580:CA>C", "B2M", "FrameShift", "p.L15fs")
    v_tap1 = _fake_variant("6:32845000:G>A", "TAP1", "Substitution", "p.R412Q")
    v_jak1 = _fake_variant("1:64850000:C>T", "JAK1", "Silent", "p.A50A")
    v_tp53 = _fake_variant("17:7675088:C>T", "TP53", "Substitution", "p.R175H")

    iso_b2m = SimpleNamespace(variant=v_b2m, num_alt_supporting_reads=0, num_ref_supporting_reads=12)
    iso_tap1 = SimpleNamespace(variant=v_tap1, num_alt_supporting_reads=14, num_ref_supporting_reads=8)

    det = detect_dna_mhc_and_apm_variants(
        [v_b2m, v_tap1, v_jak1, v_tp53],
        isovar_results=[iso_b2m, iso_tap1],
        dna_vaf_by_variant={v_b2m.short_description: {"tumor_vaf": 0.42}},
    )
    assert det.has_findings is True
    assert det.has_disruptive_findings is True
    assert len(det.variants) == 2
    assert len(det.silent_variants) == 1  # JAK1 silent tracked separately, never disrupts APM

    by_gene = {v.gene: v for v in det.variants}
    assert by_gene["B2M"].impact_class == "truncating_lof"
    assert by_gene["B2M"].dna_vaf == pytest.approx(0.42)
    assert by_gene["B2M"].rna_alt_reads == 0
    assert by_gene["B2M"].rna_confirmed is False  # e.g. NMD — DNA LoF still recorded!
    assert by_gene["TAP1"].impact_class == "missense_or_inframe"
    assert by_gene["TAP1"].rna_confirmed is True

    bundle = det.to_mhc_allele_status_bundle()
    assert bundle.apm is not None
    # B2M truncating_lof disrupts APM; TAP1 missense does NOT falsely flag as disrupted; JAK1 silent excluded
    assert bundle.apm.disrupted_genes() == ["B2M"]
    genes_by_name = {g.gene: g for g in bundle.apm.genes}
    assert genes_by_name["B2M"].dna_status == "truncating_lof"
    assert genes_by_name["B2M"].rna_status == "unknown"  # strictly separated from DNA!
    assert genes_by_name["TAP1"].dna_status == "missense_or_inframe"
    assert genes_by_name["TAP1"].is_disrupted is False


def test_canine_and_murine_classical_mhc_locus_homozygous_vs_heterozygous_phasing():
    # Canine homozygous DLA-88*034:01 with nonsense somatic variant -> auto-masks allele!
    v_dla = _fake_variant(
        "12:123456:C>T", "DLA-88", "PrematureStop", "p.Q112*",
        species="canis_lupus_familiaris",
    )
    det_hom = detect_dna_mhc_and_apm_variants(
        [v_dla],
        patient_mhc_alleles=["DLA-88*034:01", "DLA-12*001:01"],
    )
    assert det_hom.variants[0].phasing_status == "homozygous_mhc_locus"
    bundle_hom = det_hom.to_mhc_allele_status_bundle()
    assert len(bundle_hom.alleles) == 1
    assert bundle_hom.alleles[0].allele == "DLA-88*034:01"
    dec = resolve_effective_status(bundle_hom.alleles[0])
    assert dec.dna_status == "truncating_lof"
    assert dec.effective_status == "lost"
    assert dec.loss_mechanism == "dna_truncating_lof"
    assert dec.masked is True

    # Canine heterozygous DLA-88 -> do NOT guess which allele was hit!
    det_het = detect_dna_mhc_and_apm_variants(
        [v_dla],
        patient_mhc_alleles=["DLA-88*034:01", "DLA-88*508:01"],
    )
    assert det_het.variants[0].phasing_status == "unphased_heterozygous_mhc_locus"
    bundle_het = det_het.to_mhc_allele_status_bundle()
    assert bundle_het.alleles == ()  # does not blindly mask either heterozygous allele
    assert bundle_het.apm is not None
    assert "DLA-88" in bundle_het.apm.disrupted_genes()
    assert "unphased across heterozygous alleles" in bundle_het.apm.genes[0].notes

    # Mouse H2-K1 frameshift on homozygous H2-Kb
    v_mouse = _fake_variant(
        "17:33996000:GA>G", "H2-K1", "FrameShift", "p.R45fs",
        species="mus_musculus",
    )
    det_mouse = detect_dna_mhc_and_apm_variants(
        [v_mouse],
        patient_mhc_alleles=["H2-Kb", "H2-Db"],
    )
    bundle_mouse = det_mouse.to_mhc_allele_status_bundle()
    assert len(bundle_mouse.alleles) == 1
    assert bundle_mouse.alleles[0].allele == "H2-K*b"
    assert resolve_effective_status(bundle_mouse.alleles[0]).masked is True


def test_merge_preserves_both_dna_detected_and_user_rna_apm_evidence():
    v_b2m = _fake_variant("15:44711580:CA>C", "B2M", "FrameShift", "p.L15fs")
    dna_bundle = detect_dna_mhc_and_apm_variants([v_b2m]).to_mhc_allele_status_bundle()

    rna_bundle = MHCAlleleStatusBundle(
        apm=APMStatus(
            genes=(
                APMGeneStatus(
                    gene="B2M",
                    dna_status="unknown",
                    rna_status="silenced",
                    rna_ntpm=0.5,
                    immune_control_ntpm=85.0,
                    assay_source="scrna_tumor_cluster",
                ),
            ),
        ),
        source_path="user_rna.json",
    )
    merged = merge_mhc_allele_status_bundles(dna_bundle, rna_bundle)
    assert merged is not None and merged.apm is not None
    b2m = merged.apm.genes[0]
    assert b2m.dna_status == "truncating_lof"
    assert b2m.rna_status == "silenced"
    assert b2m.rna_ntpm == pytest.approx(0.5)
    assert len(b2m.dna_variants) == 1


def test_end_to_end_apm_detection_renders_in_reports_and_tumor_board_dossier():
    from .test_tumor_board import _make_candidate

    v_b2m = _fake_variant("15:44711580:CA>C", "B2M", "FrameShift", "p.L15fs")
    bundle = detect_dna_mhc_and_apm_variants([v_b2m]).to_mhc_allele_status_bundle()

    c1 = _make_candidate("TP53_p.R175H", "EVVRHCPHHERCSDSDGLAPPQHLI", {"HLA-A*02:01": 4.0}, n_alt_reads=25)
    ranked, audit = apply_mhc_allele_status_to_ranked([c1], bundle)
    assert audit["loss_of_apm"]["has_loss"] is True
    assert audit["loss_of_apm"]["disrupted_genes"] == ["B2M"]

    env = Environment(loader=PackageLoader("vaxrank", "templates"))
    html = env.get_template("mhc_allele_status.html").render(mhc_allele_status=audit)
    txt = env.get_template("mhc_allele_status.txt").render(mhc_allele_status=audit)
    assert "Loss of APM: B2M" in html
    assert "truncating_lof" in html
    assert "LOSS OF APM: B2M" in txt

    dossier = build_tumor_board_dossier(
        ranked,
        patient_id="PT-APM",
        mhc_alleles=["HLA-A*02:01"],
        mhc_allele_status=audit,
    )
    cand = dossier["candidates"][0]
    assert cand["tumor_evolution_and_mhc"]["apm_disrupted_genes"] == ["B2M"]
    assert any(f["flag"] == "apm_pathway_disruption" for f in cand["review_flags"])


def test_apm_variants_yaml_policy_and_cli_opt_out():
    cfg = load_vaxrank_config(config_path="builtin:default")
    pol = APMVariantPolicy.resolve(cfg["apm_variants"])
    assert pol.enabled is True
    assert pol.auto_mask_mhc_gene_lof_when_homozygous_locus is True

    with pytest.raises(ValueError, match="Unknown apm_variants keys"):
        APMVariantPolicy.resolve({"bad_key": True})

    args_default = parse_vaxrank_args([
        "--vcf", "a.vcf", "--bam", "a.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
    ])
    assert args_default.detect_apm_variants is True

    args_off = parse_vaxrank_args([
        "--vcf", "a.vcf", "--bam", "a.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--no-detect-apm-variants",
    ])
    assert args_off.detect_apm_variants is False
