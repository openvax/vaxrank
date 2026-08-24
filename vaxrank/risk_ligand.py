# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Patient-HLA ligand index over oncoref-derived tissue-risk proteins."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from typing import Any, Optional

from serializable import DataclassSerializable
from topiary import TopiaryPredictor

from .amino_acids import (
    has_only_standard_amino_acids,
    validate_amino_acid_sequence,
)
from .identifiers import normalize_ensembl_gene_id, normalize_mhc_allele
from .prediction_input import finite_prediction_value, prediction_integer
from .tissue_risk import (
    TissueRiskAssessment,
    TissueRiskEvidence,
)


PROTEIN_RESOLUTION_ALL_ISOFORMS = "all_protein_coding_isoforms"
PROTEIN_RESOLUTION_LONGEST_ISOFORM = "longest_protein_coding_isoform"
PROTEIN_RESOLUTION_POLICIES = frozenset({
    PROTEIN_RESOLUTION_ALL_ISOFORMS,
    PROTEIN_RESOLUTION_LONGEST_ISOFORM,
})
RISK_COMBINATION_UNION = "union"
RISK_COMBINATION_INTERSECTION = "intersection"
RISK_COMBINATION_POLICIES = frozenset({
    RISK_COMBINATION_UNION,
    RISK_COMBINATION_INTERSECTION,
})


class RiskLigandError(RuntimeError):
    """Raised when a risk-ligand index cannot be built unambiguously."""


@dataclass(frozen=True)
class ResolvedRiskProtein(DataclassSerializable):
    """One distinct risk-protein sequence with all transcript provenance."""

    gene_id: str
    gene_name: str
    amino_acids: str
    transcript_ids: tuple[str, ...]
    protein_ids: tuple[str, ...]
    tissue_evidence: tuple[TissueRiskEvidence, ...]

    def __post_init__(self):
        gene_id = normalize_ensembl_gene_id(self.gene_id)
        if not gene_id:
            raise ValueError(
                "Resolved risk protein requires a canonical amino-acid sequence"
            )
        try:
            validate_amino_acid_sequence(
                self.amino_acids, "Resolved risk protein"
            )
        except ValueError as error:
            raise ValueError(
                "Resolved risk protein requires a canonical amino-acid sequence"
            ) from error
        if not self.transcript_ids or not self.tissue_evidence:
            raise ValueError(
                "Resolved risk protein requires transcript and tissue provenance"
            )
        object.__setattr__(self, "gene_id", gene_id)
        object.__setattr__(
            self, "transcript_ids", tuple(sorted(set(self.transcript_ids)))
        )
        object.__setattr__(
            self,
            "protein_ids",
            tuple(sorted(set(protein_id for protein_id in self.protein_ids
                             if protein_id))),
        )
        object.__setattr__(
            self,
            "tissue_evidence",
            tuple(sorted(set(self.tissue_evidence), key=lambda item: (
                item.tissue_group,
                item.hpa_tissue,
                item.cell_type,
                item.level,
                item.reliability,
            ))),
        )


@dataclass(frozen=True)
class ProteinResolutionCoverage(DataclassSerializable):
    requested_gene_count: int
    resolved_gene_count: int
    unresolved_gene_ids: tuple[str, ...]
    protein_coding_transcript_count: int
    invalid_sequence_transcript_count: int
    distinct_sequence_count: int

    def __post_init__(self):
        object.__setattr__(
            self, "unresolved_gene_ids", tuple(sorted(set(self.unresolved_gene_ids)))
        )
        if self.requested_gene_count < self.resolved_gene_count:
            raise ValueError("Resolved gene count exceeds requested genes")
        if (
            self.resolved_gene_count
            + len(self.unresolved_gene_ids)
            != self.requested_gene_count
        ):
            raise ValueError("Protein-resolution coverage does not reconcile")


@dataclass(frozen=True)
class RiskProteinSequenceSet(DataclassSerializable):
    resolution_policy: str
    genome_release: str
    species: str
    tissue_risk_policy_sha256: str
    hpa_source_sha256: str
    cta_unfiltered_gene_ids_sha256: str
    proteins: tuple[ResolvedRiskProtein, ...]
    coverage: ProteinResolutionCoverage

    def __post_init__(self):
        if self.resolution_policy not in PROTEIN_RESOLUTION_POLICIES:
            raise ValueError("Unknown protein resolution policy")
        for digest in (
            self.tissue_risk_policy_sha256,
            self.hpa_source_sha256,
            self.cta_unfiltered_gene_ids_sha256,
        ):
            if len(digest) != 64:
                raise ValueError("Risk protein provenance requires SHA-256 digests")
        proteins = tuple(sorted(self.proteins, key=lambda protein: (
            protein.gene_id,
            protein.transcript_ids,
            protein.amino_acids,
        )))
        identities = [
            (protein.gene_id, protein.amino_acids) for protein in proteins
        ]
        if len(identities) != len(set(identities)):
            raise ValueError("Risk protein set has duplicate gene/sequence entries")
        if len(proteins) != self.coverage.distinct_sequence_count:
            raise ValueError("Risk protein sequence count disagrees with coverage")
        object.__setattr__(self, "proteins", proteins)

    def named_sequences(self) -> dict[str, str]:
        return {
            f"risk_{index:06d}_{protein.gene_id}": protein.amino_acids
            for index, protein in enumerate(self.proteins)
        }

    def sources_by_name(self) -> dict[str, ResolvedRiskProtein]:
        return {
            name: protein
            for name, protein in zip(self.named_sequences(), self.proteins)
        }


@dataclass(frozen=True)
class RiskPercentileThreshold(DataclassSerializable):
    kind: str
    max_percentile_rank: float = 1.0

    def __post_init__(self):
        if not self.kind:
            raise ValueError("Risk prediction kind is required")
        if not 0.0 <= self.max_percentile_rank <= 100.0:
            raise ValueError("Risk percentile threshold must be between 0 and 100")


@dataclass(frozen=True)
class RiskLigandSelectionPolicy(DataclassSerializable):
    thresholds: tuple[RiskPercentileThreshold, ...] = (
        RiskPercentileThreshold("pMHC_affinity", 1.0),
        RiskPercentileThreshold("pMHC_presentation", 1.0),
    )
    combination: str = RISK_COMBINATION_UNION

    def __post_init__(self):
        thresholds = tuple(sorted(self.thresholds, key=lambda item: item.kind))
        kinds = [threshold.kind for threshold in thresholds]
        if not thresholds or len(kinds) != len(set(kinds)):
            raise ValueError("Risk selection requires unique prediction kinds")
        if self.combination not in RISK_COMBINATION_POLICIES:
            raise ValueError("Unknown risk model-combination policy")
        object.__setattr__(self, "thresholds", thresholds)

    @property
    def threshold_by_kind(self) -> dict[str, float]:
        return {
            threshold.kind: threshold.max_percentile_rank
            for threshold in self.thresholds
        }


@dataclass(frozen=True)
class RiskLigandPrediction(DataclassSerializable):
    kind: str
    predictor_name: str
    predictor_version: str
    allele: str
    score: Optional[float]
    value: Optional[float]
    percentile_rank: Optional[float]

    def __post_init__(self):
        if not self.kind or not self.predictor_name or not self.predictor_version:
            raise ValueError(
                "Risk prediction requires kind, predictor name, and version"
            )
        for label in ("score", "value", "percentile_rank"):
            value = getattr(self, label)
            if value is not None and not math.isfinite(value):
                raise ValueError(f"Risk prediction {label} must be finite or None")
        if (
            self.percentile_rank is not None
            and not 0.0 <= self.percentile_rank <= 100.0
        ):
            raise ValueError("Risk prediction percentile rank must be 0..100")
        object.__setattr__(self, "allele", normalize_mhc_allele(self.allele))

    @classmethod
    def from_prediction_row(cls, row: Any) -> "RiskLigandPrediction":
        """Build risk-ligand evidence from one topiary prediction row."""
        value = finite_prediction_value(row.get("affinity"))
        if value is None:
            value = finite_prediction_value(row.get("value"))
        return cls(
            kind=str(row.get("kind") or ""),
            predictor_name=str(row.get("prediction_method_name") or ""),
            predictor_version=str(row.get("predictor_version") or ""),
            allele=str(row.get("allele") or ""),
            score=finite_prediction_value(row.get("score")),
            value=value,
            percentile_rank=finite_prediction_value(row.get("percentile_rank")),
        )

    @property
    def identity(self) -> tuple[str, str, str, str]:
        return (
            self.kind,
            self.predictor_name,
            self.predictor_version,
            self.allele,
        )


@dataclass(frozen=True)
class RiskLigandSource(DataclassSerializable):
    source_sequence_name: str
    protein_sequence_sha256: str
    gene_id: str
    gene_name: str
    transcript_ids: tuple[str, ...]
    protein_ids: tuple[str, ...]
    protein_offsets: tuple[int, ...]
    tissue_evidence: tuple[TissueRiskEvidence, ...]

    def __post_init__(self):
        if not self.source_sequence_name or len(self.protein_sequence_sha256) != 64:
            raise ValueError("Risk ligand source requires protein sequence identity")
        object.__setattr__(
            self, "gene_id", normalize_ensembl_gene_id(self.gene_id)
        )
        object.__setattr__(
            self, "transcript_ids", tuple(sorted(set(self.transcript_ids)))
        )
        object.__setattr__(self, "protein_ids", tuple(sorted(set(self.protein_ids))))
        object.__setattr__(
            self, "protein_offsets", tuple(sorted(set(self.protein_offsets)))
        )
        object.__setattr__(
            self,
            "tissue_evidence",
            tuple(sorted(set(self.tissue_evidence), key=lambda item: (
                item.tissue_group,
                item.hpa_tissue,
                item.cell_type,
                item.level,
                item.reliability,
            ))),
        )


@dataclass(frozen=True)
class RiskLigand(DataclassSerializable):
    peptide: str
    allele: str
    predictions: tuple[RiskLigandPrediction, ...]
    inclusion_kinds: tuple[str, ...]
    sources: tuple[RiskLigandSource, ...]

    def __post_init__(self):
        if not self.peptide or not self.allele or not self.sources:
            raise ValueError("Risk ligand requires peptide, allele, and sources")
        allele = normalize_mhc_allele(self.allele)
        predictions = tuple(sorted(self.predictions, key=lambda item: item.identity))
        identities = [prediction.identity for prediction in predictions]
        if len(identities) != len(set(identities)):
            raise ValueError("Risk ligand has duplicate prediction identities")
        if any(prediction.allele != allele for prediction in predictions):
            raise ValueError("Risk ligand and prediction alleles disagree")
        object.__setattr__(self, "allele", allele)
        object.__setattr__(self, "predictions", predictions)
        object.__setattr__(
            self, "inclusion_kinds", tuple(sorted(set(self.inclusion_kinds)))
        )
        object.__setattr__(
            self,
            "sources",
            tuple(sorted(self.sources, key=lambda item: (
                item.gene_id,
                item.source_sequence_name,
                item.transcript_ids,
                item.protein_offsets,
            ))),
        )

    @property
    def peptide_length(self) -> int:
        return len(self.peptide)


@dataclass(frozen=True)
class RiskCoverageGap(DataclassSerializable):
    allele: str
    peptide_length: int
    kind: str
    expected_window_count: int
    observed_window_count: int

    def __post_init__(self):
        if (
            isinstance(self.expected_window_count, bool)
            or isinstance(self.observed_window_count, bool)
            or self.observed_window_count < 0
            or self.expected_window_count <= self.observed_window_count
        ):
            raise ValueError("Risk coverage gaps require missing prediction windows")
        object.__setattr__(self, "allele", normalize_mhc_allele(self.allele))

    @property
    def missing_window_count(self) -> int:
        return self.expected_window_count - self.observed_window_count


@dataclass(frozen=True)
class RiskLigandCoverage(DataclassSerializable):
    prediction_row_count: int
    configured_prediction_row_count: int
    passing_prediction_row_count: int
    retained_ligand_count: int
    missing_combinations: tuple[RiskCoverageGap, ...] = field(default_factory=tuple)

    def __post_init__(self):
        counts = (
            self.prediction_row_count,
            self.configured_prediction_row_count,
            self.passing_prediction_row_count,
            self.retained_ligand_count,
        )
        if any(isinstance(count, bool) or count < 0 for count in counts):
            raise ValueError("Risk ligand coverage counts cannot be negative")
        if not (
            self.prediction_row_count >= self.configured_prediction_row_count
            >= self.passing_prediction_row_count
        ):
            raise ValueError("Risk ligand prediction coverage does not reconcile")
        object.__setattr__(
            self,
            "missing_combinations",
            tuple(sorted(set(self.missing_combinations), key=lambda item: (
                item.allele, item.peptide_length, item.kind
            ))),
        )


@dataclass(frozen=True)
class RiskLigandIndexProvenance(DataclassSerializable):
    genome_release: str
    species: str
    protein_resolution_policy: str
    tissue_risk_policy_sha256: str
    hpa_source_sha256: str
    cta_unfiltered_gene_ids_sha256: str
    predictor_identities: tuple[str, ...]
    cache_identity_sha256: str

    def __post_init__(self):
        object.__setattr__(
            self, "predictor_identities", tuple(sorted(set(self.predictor_identities)))
        )
        for digest in (
            self.tissue_risk_policy_sha256,
            self.hpa_source_sha256,
            self.cta_unfiltered_gene_ids_sha256,
            self.cache_identity_sha256,
        ):
            if len(digest) != 64:
                raise ValueError("Risk ligand provenance requires SHA-256 digests")


@dataclass(frozen=True)
class PatientHLARiskLigandIndex(DataclassSerializable):
    alleles: tuple[str, ...]
    peptide_lengths: tuple[int, ...]
    selection_policy: RiskLigandSelectionPolicy
    ligands: tuple[RiskLigand, ...]
    protein_resolution_coverage: ProteinResolutionCoverage
    coverage: RiskLigandCoverage
    provenance: RiskLigandIndexProvenance

    def __post_init__(self):
        alleles = tuple(sorted({
            normalize_mhc_allele(allele) for allele in self.alleles
        }))
        lengths = tuple(sorted(set(self.peptide_lengths)))
        if not alleles or not lengths or any(length <= 0 for length in lengths):
            raise ValueError("Risk index requires alleles and peptide lengths")
        ligands = tuple(sorted(self.ligands, key=lambda item: (
            item.allele, len(item.peptide), item.peptide
        )))
        identities = [
            (ligand.allele, ligand.peptide) for ligand in ligands
        ]
        if len(identities) != len(set(identities)):
            raise ValueError("Risk index has duplicate allele/peptide entries")
        if len(ligands) != self.coverage.retained_ligand_count:
            raise ValueError("Risk index ligand count disagrees with coverage")
        for ligand in ligands:
            if ligand.allele not in alleles or len(ligand.peptide) not in lengths:
                raise ValueError("Risk ligand is outside the configured index")
        object.__setattr__(self, "alleles", alleles)
        object.__setattr__(self, "peptide_lengths", lengths)
        object.__setattr__(self, "ligands", ligands)

    def for_allele_and_length(self, allele: str, peptide_length: int):
        normalized = normalize_mhc_allele(allele)
        return tuple(
            ligand for ligand in self.ligands
            if ligand.allele == normalized and len(ligand.peptide) == peptide_length
        )

    def to_report_dict(self) -> dict:
        return json.loads(json.dumps(asdict(self)))


def resolve_tissue_risk_protein_sequences(
    tissue_risk: TissueRiskAssessment,
    genome,
    *,
    resolution_policy: str = PROTEIN_RESOLUTION_ALL_ISOFORMS,
) -> RiskProteinSequenceSet:
    """Resolve risk genes to Ensembl proteins without inventing an HPA isoform.

    HPA IHC is gene-level. The safety-conservative default therefore scans every
    distinct protein-coding isoform; callers may choose the deterministic longest
    isoform policy as an explicit performance tradeoff.
    """
    if resolution_policy not in PROTEIN_RESOLUTION_POLICIES:
        raise ValueError("Unknown protein resolution policy")
    proteins = []
    unresolved = []
    transcript_count = 0
    invalid_count = 0
    resolved_gene_count = 0
    for risk_protein in tissue_risk.proteins:
        try:
            gene = genome.gene_by_id(risk_protein.gene_id)
            transcripts = tuple(gene.transcripts)
        except Exception:
            transcripts = ()
        by_sequence: dict[str, dict[str, set[str]]] = {}
        for transcript in transcripts:
            if not bool(getattr(transcript, "is_protein_coding", False)):
                continue
            transcript_count += 1
            try:
                sequence = (transcript.protein_sequence or "").rstrip("*")
            except Exception:
                sequence = ""
            if not sequence or not has_only_standard_amino_acids(sequence):
                invalid_count += 1
                continue
            slot = by_sequence.setdefault(
                sequence, {"transcript_ids": set(), "protein_ids": set()}
            )
            transcript_id = str(
                getattr(transcript, "transcript_id", "")
                or getattr(transcript, "id", "")
            )
            if transcript_id:
                slot["transcript_ids"].add(transcript_id)
            protein_id = str(getattr(transcript, "protein_id", "") or "")
            if protein_id:
                slot["protein_ids"].add(protein_id)
        if not by_sequence:
            unresolved.append(risk_protein.gene_id)
            continue
        resolved_gene_count += 1
        sequence_items = list(by_sequence.items())
        if resolution_policy == PROTEIN_RESOLUTION_LONGEST_ISOFORM:
            sequence_items = [sorted(
                sequence_items,
                key=lambda item: (
                    -len(item[0]),
                    tuple(sorted(item[1]["transcript_ids"])),
                    item[0],
                ),
            )[0]]
        for sequence, source in sequence_items:
            proteins.append(ResolvedRiskProtein(
                gene_id=risk_protein.gene_id,
                gene_name=risk_protein.gene_name,
                amino_acids=sequence,
                transcript_ids=tuple(source["transcript_ids"]),
                protein_ids=tuple(source["protein_ids"]),
                tissue_evidence=risk_protein.evidence,
            ))

    policy_payload = {
        "policy": asdict(tissue_risk.policy),
        "cta_count": tissue_risk.cta_unfiltered_gene_count,
    }
    policy_sha256 = hashlib.sha256(json.dumps(
        policy_payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    species = str(getattr(getattr(genome, "species", None), "latin_name", "") or "")
    return RiskProteinSequenceSet(
        resolution_policy=resolution_policy,
        genome_release=str(getattr(genome, "release", "") or ""),
        species=species,
        tissue_risk_policy_sha256=policy_sha256,
        hpa_source_sha256=tissue_risk.hpa_provenance.source_sha256,
        cta_unfiltered_gene_ids_sha256=(
            tissue_risk.cta_unfiltered_gene_ids_sha256
        ),
        proteins=tuple(proteins),
        coverage=ProteinResolutionCoverage(
            requested_gene_count=len(tissue_risk.proteins),
            resolved_gene_count=resolved_gene_count,
            unresolved_gene_ids=tuple(unresolved),
            protein_coding_transcript_count=transcript_count,
            invalid_sequence_transcript_count=invalid_count,
            distinct_sequence_count=len(proteins),
        ),
    )


def risk_ligand_index_from_prediction_frame(
    predictions_df,
    *,
    protein_sequences: RiskProteinSequenceSet,
    alleles: tuple[str, ...],
    peptide_lengths: tuple[int, ...],
    selection_policy: Optional[RiskLigandSelectionPolicy] = None,
) -> PatientHLARiskLigandIndex:
    """Build a patient-specific index from unfiltered topiary output."""
    selection_policy = selection_policy or RiskLigandSelectionPolicy()
    try:
        normalized_alleles = tuple(sorted({
            normalize_mhc_allele(allele) for allele in alleles
        }))
    except ValueError as error:
        raise RiskLigandError(str(error)) from error
    try:
        lengths = tuple(sorted({
            prediction_integer(length, "Configured peptide length")
            for length in peptide_lengths
        }))
    except ValueError as error:
        raise RiskLigandError(str(error)) from error
    if not normalized_alleles or not lengths or any(length <= 0 for length in lengths):
        raise ValueError("Risk index requires alleles and positive peptide lengths")
    sources_by_name = protein_sequences.sources_by_name()
    thresholds = selection_policy.threshold_by_kind
    groups: dict[tuple[str, str], dict[str, Any]] = {}
    observed_windows: dict[tuple[str, int, str], set[tuple[str, int]]] = {}
    predictor_identities = set()
    total_rows = 0
    configured_rows = 0
    passing_rows = 0
    if predictions_df is not None:
        for _, row in predictions_df.iterrows():
            total_rows += 1
            source_name = str(row.get("source_sequence_name") or "")
            source = sources_by_name.get(source_name)
            if source is None:
                raise RiskLigandError(
                    "Risk prediction source does not match resolved proteins"
                )
            peptide = str(row.get("peptide") or "")
            try:
                length = prediction_integer(
                    row.get("peptide_length"), "Peptide length"
                )
                offset = prediction_integer(
                    row.get("peptide_offset"), "Peptide offset"
                )
            except ValueError as error:
                raise RiskLigandError(str(error)) from error
            if length != len(peptide) or length not in lengths:
                raise RiskLigandError(
                    "Risk prediction peptide length is inconsistent or unconfigured"
                )
            if (
                offset < 0
                or offset + length > len(source.amino_acids)
                or source.amino_acids[offset:offset + length] != peptide
            ):
                raise RiskLigandError(
                    "Risk prediction peptide does not match its source protein"
                )
            try:
                allele = normalize_mhc_allele(row.get("allele"))
            except ValueError as error:
                raise RiskLigandError(str(error)) from error
            if allele not in normalized_alleles:
                raise RiskLigandError("Risk predictor emitted an unconfigured allele")
            try:
                prediction = RiskLigandPrediction.from_prediction_row(row)
            except ValueError as error:
                raise RiskLigandError(str(error)) from error
            if prediction.kind not in thresholds:
                continue
            configured_rows += 1
            if prediction.percentile_rank is not None:
                observed_windows.setdefault(
                    (allele, length, prediction.kind), set()
                ).add((source_name, offset))
            predictor_identities.add(
                "|".join((
                    prediction.kind,
                    prediction.predictor_name,
                    prediction.predictor_version,
                ))
            )
            passed = (
                prediction.percentile_rank is not None
                and prediction.percentile_rank <= thresholds[prediction.kind]
            )
            if passed:
                passing_rows += 1
            group = groups.setdefault(
                (peptide, allele),
                {
                    "predictions": {},
                    "passing_kinds": set(),
                    "sources": {},
                },
            )
            existing_prediction = group["predictions"].get(prediction.identity)
            if existing_prediction is not None and existing_prediction != prediction:
                raise RiskLigandError(
                    "Conflicting risk predictions for one model/allele/peptide"
                )
            group["predictions"][prediction.identity] = prediction
            if passed:
                group["passing_kinds"].add(prediction.kind)
            source_slot = group["sources"].setdefault(
                source_name, {"source": source, "offsets": set()}
            )
            source_slot["offsets"].add(offset)

    required_kinds = set(thresholds)
    ligands = []
    for (peptide, allele), group in groups.items():
        passing_kinds = group["passing_kinds"]
        included = (
            bool(passing_kinds)
            if selection_policy.combination == RISK_COMBINATION_UNION
            else required_kinds <= passing_kinds
        )
        if not included:
            continue
        sources = tuple(
            RiskLigandSource(
                source_sequence_name=source_name,
                protein_sequence_sha256=hashlib.sha256(
                    slot["source"].amino_acids.encode("ascii")
                ).hexdigest(),
                gene_id=slot["source"].gene_id,
                gene_name=slot["source"].gene_name,
                transcript_ids=slot["source"].transcript_ids,
                protein_ids=slot["source"].protein_ids,
                protein_offsets=tuple(slot["offsets"]),
                tissue_evidence=slot["source"].tissue_evidence,
            )
            for slot in group["sources"].values()
        )
        ligands.append(RiskLigand(
            peptide=peptide,
            allele=allele,
            predictions=tuple(group["predictions"].values()),
            inclusion_kinds=tuple(passing_kinds),
            sources=sources,
        ))

    expected_windows_by_length = {
        length: sum(
            max(0, len(protein.amino_acids) - length + 1)
            for protein in protein_sequences.proteins
        )
        for length in lengths
    }
    missing = tuple(
        RiskCoverageGap(
            allele=allele,
            peptide_length=length,
            kind=kind,
            expected_window_count=expected_windows_by_length[length],
            observed_window_count=len(observed_windows.get((allele, length, kind), ())),
        )
        for allele in normalized_alleles
        for length in lengths
        for kind in sorted(required_kinds)
        if len(observed_windows.get((allele, length, kind), ()))
        < expected_windows_by_length[length]
    )
    cache_payload = {
        "protein_provenance": {
            "resolution_policy": protein_sequences.resolution_policy,
            "genome_release": protein_sequences.genome_release,
            "species": protein_sequences.species,
            "tissue_policy": protein_sequences.tissue_risk_policy_sha256,
            "hpa": protein_sequences.hpa_source_sha256,
            "cta": protein_sequences.cta_unfiltered_gene_ids_sha256,
            "proteins": [
                {
                    "gene_id": protein.gene_id,
                    "transcripts": protein.transcript_ids,
                    "proteins": protein.protein_ids,
                    "sequence_sha256": hashlib.sha256(
                        protein.amino_acids.encode("ascii")
                    ).hexdigest(),
                }
                for protein in protein_sequences.proteins
            ],
        },
        "alleles": normalized_alleles,
        "peptide_lengths": lengths,
        "selection_policy": asdict(selection_policy),
        "predictors": sorted(predictor_identities),
    }
    cache_identity_sha256 = hashlib.sha256(json.dumps(
        cache_payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    provenance = RiskLigandIndexProvenance(
        genome_release=protein_sequences.genome_release,
        species=protein_sequences.species,
        protein_resolution_policy=protein_sequences.resolution_policy,
        tissue_risk_policy_sha256=protein_sequences.tissue_risk_policy_sha256,
        hpa_source_sha256=protein_sequences.hpa_source_sha256,
        cta_unfiltered_gene_ids_sha256=(
            protein_sequences.cta_unfiltered_gene_ids_sha256
        ),
        predictor_identities=tuple(predictor_identities),
        cache_identity_sha256=cache_identity_sha256,
    )
    return PatientHLARiskLigandIndex(
        alleles=normalized_alleles,
        peptide_lengths=lengths,
        selection_policy=selection_policy,
        ligands=tuple(ligands),
        protein_resolution_coverage=protein_sequences.coverage,
        coverage=RiskLigandCoverage(
            prediction_row_count=total_rows,
            configured_prediction_row_count=configured_rows,
            passing_prediction_row_count=passing_rows,
            retained_ligand_count=len(ligands),
            missing_combinations=missing,
        ),
        provenance=provenance,
    )


def build_patient_hla_risk_ligand_index(
    mhc_predictor,
    protein_sequences: RiskProteinSequenceSet,
    *,
    alleles: tuple[str, ...],
    peptide_lengths: tuple[int, ...],
    selection_policy: Optional[RiskLigandSelectionPolicy] = None,
) -> PatientHLARiskLigandIndex:
    """Predict and select patient-HLA risk ligands without raw-value mixing."""
    predictor = (
        mhc_predictor
        if isinstance(mhc_predictor, TopiaryPredictor)
        else TopiaryPredictor(models=[mhc_predictor])
    )
    try:
        predictions_df = predictor.predict_from_named_sequences(
            protein_sequences.named_sequences()
        )
    except Exception as error:
        raise RiskLigandError("MHC tissue-risk prediction failed") from error
    return risk_ligand_index_from_prediction_frame(
        predictions_df,
        protein_sequences=protein_sequences,
        alleles=alleles,
        peptide_lengths=peptide_lengths,
        selection_policy=selection_policy,
    )
