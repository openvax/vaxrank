#!/usr/bin/env python3
"""Rank the committed evidence-gated assembled osteosarcoma antigens."""

import argparse
import json
from pathlib import Path
import sys

from mhctools import NetMHCpan

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from vaxrank import (  # noqa: E402
    AminoAcidInterval,
    EpitopeConfig,
    TumorSpecificityAttestation,
    VaccineAntigen,
    vaccine_peptides_for_antigen,
)
from vaxrank.vaccine_antigen import ATTESTATION_ADMITTED  # noqa: E402


DEFAULT_ALLELES = (
    "HLA-A*01:01",
    "HLA-B*08:01",
    "HLA-B*27:05",
    "HLA-C*01:02",
    "HLA-C*07:01",
)


def _top_epitopes(vaccine_peptide):
    epitopes = sorted(
        vaccine_peptide.target_epitopes,
        key=lambda epitope: max(epitope.per_allele_scores.values()),
        reverse=True,
    )
    rows = []
    for epitope in epitopes[:8]:
        prediction = min(
            epitope.predictions_for("pMHC_affinity"),
            key=lambda value: value.value,
        )
        rows.append({
            "sequence": epitope.sequence,
            "allele": str(prediction.allele),
            "ic50_nm": prediction.value,
            "percentile_rank": prediction.percentile_rank,
            "epitope_score": max(epitope.per_allele_scores.values()),
        })
    return rows


def rank_antigens(source, alleles=DEFAULT_ALLELES):
    """Return auditable ranking records for evidence-admitted inputs."""
    payload = json.loads(Path(source).read_text())
    predictor = NetMHCpan(
        list(alleles), default_peptide_lengths=[8, 9, 10, 11])
    attestation = TumorSpecificityAttestation(
        status=ATTESTATION_ADMITTED,
        evidence_kind="sample_specific_rna_translation",
        evidence_source="committed osteosarcoma evidence record",
        patient_specific=True,
        rationale_code="rna_evidence_gate_passed",
    )
    results = []
    for record in payload["antigens"]:
        if record["evidence_gate"] != "pass":
            continue
        metadata = tuple(record["source_metadata"].items())
        antigen = VaccineAntigen.from_assembled_mutation_sequence(
            amino_acids=record["amino_acids"],
            targetable_intervals=tuple(
                AminoAcidInterval(*interval)
                for interval in record["targetable_intervals"]),
            tumor_specificity=attestation,
            gene_name=record["gene_name"],
            transcript_ids=tuple(record["transcript_ids"]),
            source_identifier=record["id"],
            source_metadata=metadata,
        )
        peptides = vaccine_peptides_for_antigen(
            antigen=antigen,
            mhc_predictor=predictor,
            vaccine_peptide_length=25,
            max_vaccine_peptides=1,
            epitope_config=EpitopeConfig(),
        )
        result = {
            "id": record["id"],
            "evidence_gate": record["evidence_gate"],
            "selected_long_peptide": None,
            "rank_score": None,
            "target_epitope_score": None,
            "top_epitopes": [],
        }
        if peptides:
            peptide = peptides[0]
            result.update({
                "selected_long_peptide": peptide.amino_acids,
                "rank_score": peptide.combined_score,
                "target_epitope_score": peptide.target_epitope_score,
                "top_epitopes": _top_epitopes(peptide),
            })
        results.append(result)
    return {
        "schema_version": 1,
        "score_basis": "source-agnostic target epitope score",
        "self_reference_screen": (
            "not evaluated; no full human reference proteome supplied"),
        "hla_alleles": list(alleles),
        "records": results,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        default=str(Path(__file__).parent / "source" / "assembled_antigens.json"),
    )
    parser.add_argument("--output")
    args = parser.parse_args(argv)
    result = rank_antigens(args.source)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
