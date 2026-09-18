#!/usr/bin/env python3
"""Combine base panels with reproducible assembled-antigen rankings."""

import json
from hashlib import sha256
from pathlib import Path


HERE = Path(__file__).parent
SOURCE = HERE / "source"


def build_results():
    component_names = (
        "results.json",
        "assembled_antigens.json",
        "assembled_rankings.json",
        "additional_candidate_audit.json",
        "orf_platform_audit.json",
        "platform_comparison_audit.json",
        "panel_metadata.json",
        "sv_audits.json",
    )
    component_bytes = {
        name: (SOURCE / name).read_bytes() for name in component_names}
    base = json.loads(component_bytes["results.json"])
    inputs = json.loads(component_bytes["assembled_antigens.json"])
    rankings = json.loads(component_bytes["assembled_rankings.json"])
    panels = json.loads(component_bytes["panel_metadata.json"])
    platform_comparisons = json.loads(
        component_bytes["platform_comparison_audit.json"])
    orf_platform_audit = json.loads(
        component_bytes["orf_platform_audit.json"])
    sv_audits = json.loads(component_bytes["sv_audits.json"])
    input_by_id = {record["id"]: record for record in inputs["antigens"]}
    ranking_by_id = {record["id"]: record for record in rankings["records"]}
    expanded = []
    for panel in panels["records"]:
        antigen = input_by_id[panel["id"]]
        ranking = ranking_by_id[panel["id"]]
        record = dict(panel)
        record["protein_sequence"] = antigen["amino_acids"]
        record["assembled_transcript_ids"] = antigen["transcript_ids"]
        record["targetable_intervals"] = antigen["targetable_intervals"]
        record["osteosarc_protein_context"] = antigen.get(
            "osteosarc_protein_context")
        record["osteosarc_vaccine_sequence"] = antigen.get(
            "osteosarc_vaccine_sequence")
        record["osteosarc_sequence_status"] = antigen.get(
            "osteosarc_sequence_status")
        record["selected_long_peptide"] = ranking["selected_long_peptide"]
        record["combined_score"] = ranking["rank_score"]
        record["target_epitope_score"] = ranking["target_epitope_score"]
        record["top_epitopes"] = ranking["top_epitopes"]
        record["score_basis"] = rankings["score_basis"]
        record["outcome"] = (
            "selected" if ranking["selected_long_peptide"]
            else "no_target_binder")
        expanded.append(record)
    base["metadata"]["ranking_note"] = (
        "Expanded explicit-antigen panels rank windows by source-agnostic "
        "target epitope score; RNA evidence is retained separately.")
    base["metadata"]["self_reference_screen"] = (
        "Not evaluated for explicit assembled-antigen rankings because no "
        "full human reference proteome was supplied; exploratory binding "
        "ranks only.")
    base["metadata"]["source_component_sha256"] = {
        name: sha256(data).hexdigest()
        for name, data in component_bytes.items()
    }
    base["records"] = [
        base["records"][0],
        *expanded,
        *sv_audits["records"],
        *base["records"][1:],
    ]
    base["orf_platform_audit"] = orf_platform_audit
    comparison_by_id = {
        record["id"]: record
        for record in platform_comparisons["records"]
    }
    for record in base["records"]:
        comparison = comparison_by_id.get(record["id"])
        if comparison is not None:
            record["platform_comparison"] = comparison
    return base


def main():
    output = SOURCE / "all_results.json"
    output.write_text(json.dumps(build_results(), indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    main()
