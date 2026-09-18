import json
from pathlib import Path
import struct

from vaxrank.complex_variant_visualization import (
    _balanced_chunks,
    generate_complex_variant_results,
)


def test_summary_chunks_are_balanced_without_exceeding_page_capacity():
    chunks = _balanced_chunks(list(range(22)))
    assert [len(chunk) for chunk in chunks] == [6, 6, 5, 5]


def example_record(outcome="selected"):
    selected = outcome == "selected"
    return {
        "combined_score": 4.25,
        "decision_reason": "The result follows the evidence and selection gates.",
        "dna_only": "DNA does not establish the expressed mutant protein.",
        "id": "GENE1-complex-event",
        "interpretation": "Assembly establishes what Vaxrank may predict and rank.",
        "limitation": "Research use only.",
        "outcome": outcome,
        "protein_sequence": "MABCDEFGHIJK",
        "rna_evidence": {
            "alt_fragments": 7,
            "other_fragments": 0,
            "passes_vaxrank_filters": selected,
            "platform": "Long-read RNA",
            "ref_fragments": 12,
            "source": "T1 fixture",
            "top_protein_fragments": 6,
        },
        "sample": "T1",
        "selected_long_peptide": "ABCDEFGHIJK" if selected else None,
        "source_labels": ["fixture", "predictor"],
        "source_urls": ["https://example.org/evidence"],
        "summary_rna": "Assembly establishes a mutant protein.",
        "target_epitope_score": 2.5,
        "title": "GENE1 complex event",
        "top_epitopes": ([{
            "allele": "HLA-A*02:01",
            "epitope_score": 0.9,
            "ic50_nm": 32.1,
            "percentile_rank": 0.2,
            "sequence": "ABCDEFGHI",
        }] if selected else []),
        "variant": "chr1:g.100_120del",
        "variant_class": "Large deletion",
    }


def example_payload(records):
    return {
        "metadata": {
            "analysis_date": "2026-09-16",
            "hla_alleles": ["HLA-A*02:01"],
            "predictor_version": "4.2c",
            "software_versions": {
                "isovar": "1.17.0",
                "vaxrank": "3.18.2",
            },
        },
        "records": records,
    }


def add_platform_comparison(record):
    record["platform_comparison"] = {
        "conclusion": "Both platforms agree on the mutant tail.",
        "dna_model": "DNA predicts a frameshift but not its expressed context.",
        "vaccine_impact": "Long reads provide a 25-aa construct window.",
        "overview": {
            "long": "Longer translated context",
            "short": "Same tail, shorter context",
            "vaccine": "Same epitopes; long-read construct only",
        },
        "rows": [
            {
                "label": "Oxford Nanopore RNA",
                "evidence": "Nine fragments support the translation.",
                "transcript_nt": "AAA[TT]CCC",
                "protein": "ABCDEFG[HIJKL]*",
                "peptide_pool": "20 class-I windows; five 25-aa windows",
            },
            {
                "label": "Illumina RNA",
                "evidence": "Two fragments support the same translation.",
                "transcript_nt": "A[TT]CCC",
                "protein": "EFG[HIJKL]*",
                "peptide_pool": "20 class-I windows; no 25-aa window",
            },
        ],
    }
    return record


def test_generate_compact_pdf_and_high_resolution_pages(tmp_path):
    records = [
        example_record("selected"),
        example_record("no_target_binder"),
        example_record("held_out"),
        example_record("no_translation"),
    ]
    for index, record in enumerate(records):
        record["id"] += "-%s" % index
        if record["outcome"] == "no_translation":
            record["protein_sequence"] = None
    source = tmp_path / "results.json"
    source.write_text(json.dumps(example_payload(records)))
    combined = tmp_path / "combined.pdf"
    run = generate_complex_variant_results(
        source,
        tmp_path / "runs",
        timestamp="2026-09-16T120000Z",
        combined_output=combined,
    )
    assert combined.read_bytes().startswith(b"%PDF")
    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["page_count"] == 6
    assert manifest["png_dimensions"] == {"width": 3600, "height": 2280}
    pngs = sorted(run.glob("*/result.png"))
    assert len(pngs) == 6
    assert struct.unpack(">II", pngs[0].read_bytes()[16:24]) == (3600, 2280)
    records = sorted(run.glob("*/record.json"))
    assert len(records) == 4


def test_generate_paginates_large_decision_matrix(tmp_path):
    records = [example_record() for _ in range(8)]
    for index, record in enumerate(records):
        record["id"] += "-%s" % index
    source = tmp_path / "results.json"
    source.write_text(json.dumps(example_payload(records)))

    run = generate_complex_variant_results(
        source,
        tmp_path / "runs",
        timestamp="2026-09-16T120001Z",
        png_scale=1,
    )

    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["page_count"] == 11
    assert [page["name"] for page in manifest["pages"][:2]] == [
        "summary", "summary-2"]
    assert len(list(run.glob("*/record.json"))) == 8


def test_generate_adds_one_platform_overview_and_replaces_result_layout(tmp_path):
    record = add_platform_comparison(example_record())
    source = tmp_path / "results.json"
    source.write_text(json.dumps(example_payload([record])))

    run = generate_complex_variant_results(
        source,
        tmp_path / "runs",
        timestamp="2026-09-16T120002Z",
        png_scale=1,
    )

    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["page_count"] == 4
    assert [page["name"] for page in manifest["pages"]] == [
        "summary", "platform-overview", "GENE1-complex-event", "provenance"]


def test_committed_complex_inputs_cover_requested_cases_and_provenance():
    source = (
        Path(__file__).parents[1]
        / "examples" / "osteosarc_complex_results" / "source")
    antigens = json.loads((source / "assembled_antigens.json").read_text())
    rankings = json.loads((source / "assembled_rankings.json").read_text())
    sv_audits = json.loads((source / "sv_audits.json").read_text())
    audit = json.loads((source / "additional_candidate_audit.json").read_text())
    by_id = {record["id"]: record for record in antigens["antigens"]}
    ranked_by_id = {record["id"]: record for record in rankings["records"]}

    expected = {
        "GTF3C5-p-Glu503_Glu506del",
        "TECPR1-p-Thr259fs",
        "GLIS3-p-Ser775fs",
        "RNF213-p-Ile1070delLeu",
        "KTN1-p-Leu340fs",
        "DYNC1H1-p-Val314Ile",
        "DYNC1H1-p-Gln3267His",
        "NAV2-p-Ala1809Val",
        "MAP2-compound-haplotype",
        "NTF3-compound-haplotype",
        "CD109-long-read-phased",
        "ZNF436-long-read-context",
    }
    assert set(by_id) == expected
    assert set(ranked_by_id) == expected
    assert all(record["evidence_gate"] == "pass" for record in by_id.values())
    assert all(record["source_metadata"] for record in by_id.values())
    assert by_id["CD109-long-read-phased"]["targetable_intervals"] == [
        [2, 3], [71, 72]]
    assert ranked_by_id["CD109-long-read-phased"][
        "selected_long_peptide"] is None
    assert "DRAH" in by_id["MAP2-compound-haplotype"]["amino_acids"]
    assert "DVSENY" in by_id["NTF3-compound-haplotype"]["amino_acids"]
    assert {
        record["id"] for record in sv_audits["records"]
    } == {
        "FOXO3-STRADA-CCDC47-unresolved",
        "PARD3B-CDKN2B-unresolved",
        "AMPH-internal-deletion-DNA-only",
        "chr21-unnamed-long-read-rich-unresolved",
        "DLG5-start-region-deletion",
        "AFF3-intronic-deletion",
        "KEAP1-intronic-noncoding-deletion",
        "GABBR1-SLC29A1-unresolved",
        "OTUD7A-FMN1-unresolved",
    }
    inventory = {
        record["id"]: record for record in sv_audits["inventory_only"]}
    assert inventory["MYO15B-chr17-75589953"][
        "published_vaccine_peptide"] == "AGRRAQAPTRVLGLAPP"
    assert "DLG5--DLG5" not in inventory
    assert by_id["GTF3C5-p-Glu503_Glu506del"]["source_metadata"][
        "direct_alt_fragments_t2_ont"] == "128"
    audit_indels = {record["gene"]: record for record in audit["indels"]}
    assert audit_indels["KTN1"]["fragment_counts"]["2bc1fc291308debb"][
        "alt"] == 10
    audit_svs = {
        record["id"]: record for record in audit["structural_variants"]}
    assert audit_svs["OTUD7A--FMN1"]["tagged_ont_paths"][
        "2bc1fc291308debb"]["cell_umi_keys"] == 6
    comparison = json.loads(
        (source / "platform_comparison_audit.json").read_text())
    comparisons = {record["id"]: record for record in comparison["records"]}
    assert comparisons["KTN1-p-Leu340fs"]["rows"][0][
        "protein_sequence"] == "QDALKKSSKGELTTLIHQLQEKDKFYSLL"
    assert len(comparisons["KTN1-p-Leu340fs"]["mhc_windows_8_11"]) == 20
    assert comparisons["KEAP1-intronic-noncoding-deletion"]["rows"][1][
        "protein"].endswith("processed transcript has no CDS")
    orf_audit = json.loads((source / "orf_platform_audit.json").read_text())
    full_by_platform = {
        record["platform"]: record
        for record in orf_audit["full_matrix"]["platform_summary"]}
    assert full_by_platform["ILMN"][
        "variants_with_validated_protein_window"] == 39
    assert full_by_platform["ONT"][
        "variants_with_validated_protein_window"] == 30
    assert full_by_platform["PacBio"][
        "variants_with_validated_protein_window"] == 3
    paired = {
        (record["platform"], record["mode"]): record
        for record in orf_audit["paired_corpus"]["platform_mode_summary"]}
    assert paired[("ILMN", "assembly_on")][
        "validated_local_protein_windows"] == 27
    assert paired[("ILMN", "assembly_off")][
        "validated_local_protein_windows"] == 27
    assert orf_audit["varcode"][
        "variants_with_exact_isolated_effect_protein"] == 44
