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

"""Regression tests for DNA-only CLI mode (#907), compact dna_effect_selection (#908),
and OncoRef whole-gene cohort proxy expression (including 5' upstream gene selection
for complex / structural variants).
"""

import math
from types import SimpleNamespace

import pytest
import varcode
from mhctools import RandomBindingPredictor
from vaxrank.epitope_config import EpitopeConfig

from vaxrank.cli import entry_point
from vaxrank.cli.arg_parser import parse_vaxrank_args
from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
from vaxrank.core_logic import vaccine_peptides_for_variant
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.proxy_expression import (
    OncoRefProxyExpressionLookup,
    resolve_proxy_expression_gene,
    validate_human_host_for_proxy_expression,
)
from vaxrank.report import TemplateDataCreator
from vaxrank.vaccine_config import (
    PROXY_EXPRESSION_COMBINED_SCORE_EXPR,
    SOURCE_AGNOSTIC_COMBINED_SCORE_EXPR,
    VaccineConfig,
)


def _walk_keys(obj):
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield key, value
            yield from _walk_keys(value)
    elif isinstance(obj, list):
        for item in obj:
            yield from _walk_keys(item)


# ---------------------------------------------------------------------------
# #908: Compact dna_effect_selection (no per-base genomic_positions bloat)
# ---------------------------------------------------------------------------

def test_dna_effect_selection_compacts_genomic_positions(human_genome_grch37):
    variant = varcode.Variant("7", 140453136, "A", "T", human_genome_grch37)  # BRAF V600E
    fragment = MutantProteinFragment.from_variant_dna(variant, protein_sequence_length=35)
    assert fragment is not None
    selection = fragment.dna_effect_selection
    assert selection is not None

    keys = {k for k, _ in _walk_keys(selection)}
    assert "genomic_positions" not in keys
    assert "genomic_positions_count" in keys
    assert "genomic_positions_sha256" in keys
    assert "exon_runs" in keys

    # Exact deterministic round-trip matching still works after both JSON codecs
    for restored in (
        MutantProteinFragment.from_json(fragment.to_json()),
        from_native_json(to_native_json(fragment), MutantProteinFragment),
    ):
        assert restored.dna_effect_selection == selection
        eff = restored.predicted_effect()
        assert eff is not None
        assert eff.short_description == "p.V600E"


# ---------------------------------------------------------------------------
# #907: DNA-only CLI without --bam + clean default combined_score_expr
# ---------------------------------------------------------------------------

def test_dna_only_without_bam_requires_opt_in(tmp_path, capsys):
    vcf = tmp_path / "one.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "7\t140453136\t.\tA\tT\t.\tPASS\t.\n"
    )
    with pytest.raises(SystemExit) as exc:
        entry_point.main([
            "--vcf", str(vcf),
            "--genome", "GRCh37",
            "--ensembl-release", "75",
            "--mhc-predictor", "random",
            "--mhc-alleles", "HLA-A*02:01",
            "--output-csv", str(tmp_path / "out.csv"),
            "--log-path", str(tmp_path / "vaxrank.log"),
        ])
    assert exc.value.code == 1
    err = capsys.readouterr().err
    assert "No --bam file supplied" in err
    assert "--allow-dna-only-fallback" in err


def test_proxy_expression_rejects_bam_combination(tmp_path, capsys):
    vcf = tmp_path / "one.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    bam = tmp_path / "fake.bam"
    bam.write_bytes(b"")
    with pytest.raises(SystemExit) as exc:
        entry_point.main([
            "--vcf", str(vcf),
            "--bam", str(bam),
            "--proxy-expression-cancer-type", "SKCM",
            "--mhc-predictor", "random",
            "--mhc-alleles", "HLA-A*02:01",
            "--output-csv", str(tmp_path / "out.csv"),
            "--log-path", str(tmp_path / "vaxrank.log"),
        ])
    assert exc.value.code == 1
    assert "cannot be combined with --bam" in capsys.readouterr().err


def test_dna_only_cli_default_combined_score_expr_resolution():
    # Pure DNA-only without proxy -> target_epitope_score
    args_dna = parse_vaxrank_args([
        "--vcf", "somatic.vcf",
        "--allow-dna-only-fallback",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert args_dna.bam is None
    assert vaccine_config_from_args(args_dna).combined_score_expr == SOURCE_AGNOSTIC_COMBINED_SCORE_EXPR

    # DNA-only with OncoRef proxy -> sqrt(proxy_expression_tpm) * target_epitope_score
    args_proxy = parse_vaxrank_args([
        "--vcf", "somatic.vcf",
        "--proxy-expression-cancer-type", "melanoma",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert vaccine_config_from_args(args_proxy).combined_score_expr == PROXY_EXPRESSION_COMBINED_SCORE_EXPR

    # Explicit user override still wins
    args_override = parse_vaxrank_args([
        "--vcf", "somatic.vcf",
        "--proxy-expression-cancer-type", "melanoma",
        "--set", "vaccine_peptides.combined_score_expr=log1p(proxy_expression_tpm) * target_epitope_score",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert (
        vaccine_config_from_args(args_override).combined_score_expr
        == "log1p(proxy_expression_tpm) * target_epitope_score"
    )


def test_dna_only_cli_runs_without_bam_and_scores_nonzero(tmp_path):
    vcf = tmp_path / "braf.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "7\t140453136\t.\tA\tT\t.\tPASS\t.\n"
    )
    out_csv = tmp_path / "ranked.csv"
    out_html = tmp_path / "report.html"
    entry_point.main([
        "--vcf", str(vcf),
        "--genome", "GRCh37",
        "--ensembl-release", "75",
        "--allow-dna-only-fallback",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*08:01",
        "--mhc-peptide-lengths", "8,9",
        "--set", "epitopes.affinity_cutoff=50000.0",
        "--set", "epitopes.min_score=0.0",
        "--no-processing-aware-annotation",
        "--output-csv", str(out_csv),
        "--output-html-report", str(out_html),
        "--log-path", str(tmp_path / "vaxrank.log"),
    ])
    import pandas as pd

    df = pd.read_csv(out_csv)
    assert not df.empty
    assert (df["combined_score"] > 0).any()
    assert math.isclose(df.iloc[0]["combined_score"], df.iloc[0]["target_epitope_score"], rel_tol=1e-6)
    html = out_html.read_text()
    assert "DNA-only" in html


# ---------------------------------------------------------------------------
# OncoRef whole-gene cohort proxy expression & 5' upstream complex variants
# ---------------------------------------------------------------------------

def test_oncoref_proxy_expression_lookup_and_guards():
    lookup = OncoRefProxyExpressionLookup("melanoma")
    assert lookup.cancer_code == "SKCM"
    sf3b1 = lookup.lookup(gene_id="ENSG00000115524.12", gene_name="SF3B1")
    gria1 = lookup.lookup(gene_name="GRIA1")
    assert sf3b1 is not None and sf3b1.tpm > 50.0
    assert gria1 is not None and gria1.tpm < 0.5

    with pytest.raises(ValueError, match="Unknown cancer type"):
        OncoRefProxyExpressionLookup("NOT_A_REAL_CANCER_TYPE")

    with pytest.raises(ValueError, match="no whole-gene reference expression"):
        OncoRefProxyExpressionLookup("BCC")

    validate_human_host_for_proxy_expression(9606, "SKCM")
    with pytest.raises(ValueError, match="human OncoRef"):
        validate_human_host_for_proxy_expression(10090, "SKCM")


def test_proxy_expression_weights_dna_only_and_preserves_provenance(human_genome_grch37):
    lookup = OncoRefProxyExpressionLookup("melanoma")
    variant = varcode.Variant("7", 140453136, "A", "T", human_genome_grch37)  # BRAF V600E
    isovar_result = SimpleNamespace(
        variant=variant,
        passes_all_filters=False,
        predicted_effect=variant.effects().top_priority_effect(),
    )
    predictor = RandomBindingPredictor(
        alleles=["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"],
        default_peptide_lengths=[8, 9],
    )
    ec = EpitopeConfig(binding_affinity_cutoff=50000.0, min_epitope_score=0.0)
    vc = VaccineConfig(combined_score_expr=PROXY_EXPRESSION_COMBINED_SCORE_EXPR)

    peptides = vaccine_peptides_for_variant(
        isovar_result=isovar_result,
        mhc_predictor=predictor,
        epitope_config=ec,
        vaccine_config=vc,
        allow_dna_only_fallback=True,
        proxy_expression_lookup=lookup,
    )
    assert peptides
    vp = peptides[0]
    frag = vp.mutant_protein_fragment
    # Strictly preserves zero patient RNA reads while attaching cohort proxy expression
    assert frag.n_alt_reads == 0
    assert frag.rna_evidence_method == ""
    assert frag.proxy_expression_cancer_code == "SKCM"
    assert frag.proxy_expression_tpm is not None and frag.proxy_expression_tpm > 1.0
    assert math.isclose(
        vp.combined_score,
        math.sqrt(frag.proxy_expression_tpm) * vp.target_epitope_score,
        rel_tol=1e-6,
    )

    # Report renders explicit cohort proxy badge (never patient RNA)
    from vaxrank.patient_info import PatientInfo
    patient = PatientInfo(
        patient_id="P1",
        vcf_paths=["braf.vcf"],
        bam_path="",
        mhc_alleles=["HLA-A*02:01"],
        num_somatic_variants=1,
        num_coding_effect_variants=1,
        num_variants_with_rna_support=0,
        num_variants_with_vaccine_peptides=1,
    )
    tdata = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=[(variant, peptides)],
        patient_info=patient,
        final_review="",
        reviewers="",
        args_for_report={
            "vcf": ["braf.vcf"],
            "bam": None,
            "allow_dna_only_fallback": True,
            "proxy_expression_cancer_type": "melanoma",
            "_proxy_expression_summary": {
                "cancer_code": "SKCM",
                "source_label": lookup.source_label,
            },
            "vaccine_type": ["peptide"],
            "manufacturability": False,
            "wt_epitopes": True,
        },
        input_json_file=None,
        candidate_details=True,
    ).compute_template_data()
    ranking_row = tdata["ranked_candidates"][0]
    assert "Proxy SKCM:" in ranking_row["rna_summary"]
    assert "OncoRef" in ranking_row["rna_summary"]
    variant_meta = tdata["variants"][0]["variant_data"]
    assert "Cohort proxy expression (not patient RNA)" in variant_meta
    assert "SKCM:" in variant_meta["Cohort proxy expression (not patient RNA)"]


def test_complex_fusion_uses_five_prime_upstream_gene_expression(human_genome_grch37, monkeypatch):
    """For a gene fusion anchored on the 3' mate, cohort proxy expression must
    resolve the 5' (upstream-of-breakpoint) promoter gene rather than the 3' partner."""
    from varcode import MutantTranscript, StructuralVariant
    from varcode.effects import GeneFusion

    # Anchor the SV on TP53 (3' partner), with CFTR as the 5' upstream promoter gene
    five_prime_tx = human_genome_grch37.transcript_by_id("ENST00000003084")  # CFTR
    three_prime_tx = human_genome_grch37.transcript_by_id("ENST00000269305")  # TP53
    variant = StructuralVariant(
        "17", 7577121, "BND", mate_contig="7", mate_start=117171168,
        genome=human_genome_grch37,
    )
    reference = three_prime_tx.protein_sequence
    fusion_protein = reference[:40] + "WQWQWQWQWQ" + reference[50:]
    model = MutantTranscript(
        reference_transcript=three_prime_tx,
        mutant_protein_sequence=fusion_protein,
        annotator_name="synthetic-fusion-test",
    )
    fusion_effect = GeneFusion(
        variant,
        transcript=three_prime_tx,
        partner_transcript=five_prime_tx,
        mutant_transcript=model,
        five_prime_transcript=five_prime_tx,
        three_prime_transcript=three_prime_tx,
    )
    monkeypatch.setattr(StructuralVariant, "effects", lambda self: [fusion_effect])

    fragment = MutantProteinFragment.from_variant_dna(variant, protein_sequence_length=35)
    assert fragment is not None
    # Fragment's primary gene_name comes from the anchored 3' transcript (TP53),
    # while proxy expression MUST resolve the 5' upstream-of-breakpoint gene (CFTR).
    assert fragment.gene_name == "TP53"
    gene_id, gene_name, is_upstream_complex = resolve_proxy_expression_gene(fragment)
    assert is_upstream_complex is True
    assert gene_name == "CFTR"
    assert gene_id == five_prime_tx.gene_id

    # Verify core_logic attaches the 5' upstream gene's proxy TPM and provenance
    fake_lookup = SimpleNamespace(
        cancer_code="SKCM",
        source_label="OncoRef · Synthetic",
        lookup=lambda *, gene_id=None, gene_name=None: (
            SimpleNamespace(tpm=42.5, q1_tpm=20.0, q3_tpm=60.0, cancer_code="SKCM")
            if gene_name == "CFTR"
            else SimpleNamespace(tpm=0.1, q1_tpm=0.0, q3_tpm=0.2, cancer_code="SKCM")
        ),
    )
    isovar_result = SimpleNamespace(
        variant=variant,
        passes_all_filters=False,
        predicted_effect=fusion_effect,
    )
    predictor = RandomBindingPredictor(
        alleles=["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"],
        default_peptide_lengths=[8, 9],
    )
    ec = EpitopeConfig(binding_affinity_cutoff=50000.0, min_epitope_score=0.0)
    peptides = vaccine_peptides_for_variant(
        isovar_result=isovar_result,
        mhc_predictor=predictor,
        epitope_config=ec,
        vaccine_config=VaccineConfig(combined_score_expr=PROXY_EXPRESSION_COMBINED_SCORE_EXPR),
        allow_dna_only_fallback=True,
        proxy_expression_lookup=fake_lookup,
    )
    assert peptides
    frag = peptides[0].mutant_protein_fragment
    assert frag.proxy_expression_gene_name == "CFTR"
    assert math.isclose(frag.proxy_expression_tpm, 42.5)
    assert "5' upstream gene CFTR" in frag.proxy_expression_source
    # Confirm JSON round-trip preserves proxy_expression_gene_name
    restored = from_native_json(to_native_json(frag), MutantProteinFragment)
    assert restored.proxy_expression_gene_name == "CFTR"
    assert math.isclose(restored.proxy_expression_tpm, 42.5)


def test_multi_segment_complex_variant_uses_first_upstream_segment(human_genome_grch37):
    """When five_prime_transcript is absent on a multi-segment MutantTranscript,
    resolve_proxy_expression_gene picks the first (5'-most) reference segment."""
    from varcode.mutant_transcript import MutantTranscript, ReferenceSegment

    upstream_tx = human_genome_grch37.transcript_by_id("ENST00000003084")  # CFTR
    downstream_tx = human_genome_grch37.transcript_by_id("ENST00000269305")  # TP53
    model = MutantTranscript(
        reference_transcript=downstream_tx,
        reference_segments=(
            ReferenceSegment(source=upstream_tx, start=0, end=300, strand="+", label="5prime"),
            ReferenceSegment(source=downstream_tx, start=200, end=600, strand="+", label="3prime"),
        ),
    )
    fake_effect = SimpleNamespace(transcript=downstream_tx, mutant_transcript=model)
    fake_fragment = SimpleNamespace(
        gene_name="TP53",
        supporting_reference_transcripts=[downstream_tx, upstream_tx],
        dna_effect_selection={"transcript_id": downstream_tx.id, "partner_transcript_id": upstream_tx.id},
        variant=SimpleNamespace(genome=human_genome_grch37),
        predicted_effect=lambda: fake_effect,
    )
    gene_id, gene_name, is_upstream_complex = resolve_proxy_expression_gene(fake_fragment)
    assert is_upstream_complex is True
    assert gene_name == "CFTR"
    assert gene_id == upstream_tx.gene_id
