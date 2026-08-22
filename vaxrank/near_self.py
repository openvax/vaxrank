# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Policy-neutral similarity evidence against patient-HLA risk ligands.

BLOSUM scores are half-bit log-odds values from Henikoff & Henikoff (1992),
doi:10.1073/pnas.89.22.10915. The embedded 20-residue BLOSUM62 table is checked
against the stable NCBI representation documented at
https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/sm__blosum62_8c.html.

Hamming/BLOSUM similarity is a screening heuristic. It is not evidence that
two peptide-MHC complexes do or do not share T-cell receptor recognition.
This module records similarity facts; enforcement thresholds belong to the
safety-policy layer.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from numbers import Integral
from typing import Optional

import numpy as np
from serializable import DataclassSerializable

from .amino_acids import STANDARD_AMINO_ACIDS, validate_amino_acid_sequence
from .identifiers import normalize_mhc_allele
from .risk_ligand import (
    PatientHLARiskLigandIndex,
    RiskCoverageGap,
    RiskLigand,
)


NEAR_SELF_REASON_NO_RISK_LIGANDS = "no_same_allele_length_risk_ligands"
NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX = "query_outside_risk_index"
NEAR_SELF_REASON_PREDICTION_COVERAGE = "risk_prediction_coverage_gap"
NEAR_SELF_REASON_PROTEIN_RESOLUTION = "risk_protein_resolution_gap"

BLOSUM62_ALPHABET = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62_ROWS = (
    (4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0),
    (-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3),
    (-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3),
    (-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3),
    (0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1),
    (-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2),
    (-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2),
    (0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3),
    (-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3),
    (-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3),
    (-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1),
    (-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2),
    (-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1),
    (-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1),
    (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2),
    (1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2),
    (0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0),
    (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3),
    (-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1),
    (0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4),
)


class NearSelfError(RuntimeError):
    """Raised when similarity evidence cannot be computed unambiguously."""


@dataclass(frozen=True)
class AminoAcidSubstitutionMatrix(DataclassSerializable):
    """Symmetric, versioned substitution matrix over canonical residues."""

    name: str
    version: str
    alphabet: str
    rows: tuple[tuple[int, ...], ...]
    citation: str
    source_url: str

    def __post_init__(self):
        if not self.name or not self.version or not self.citation or not self.source_url:
            raise ValueError("Substitution matrix requires identity and provenance")
        if set(self.alphabet) != STANDARD_AMINO_ACIDS or len(self.alphabet) != 20:
            raise ValueError("Substitution matrix must cover 20 canonical amino acids")
        if len(self.rows) != len(self.alphabet):
            raise ValueError("Substitution matrix row count does not match alphabet")
        rows = tuple(tuple(row) for row in self.rows)
        if any(len(row) != len(self.alphabet) for row in rows):
            raise ValueError("Substitution matrix must be square")
        if any(
            rows[i][j] != rows[j][i]
            for i in range(len(rows))
            for j in range(len(rows))
        ):
            raise ValueError("Substitution matrix must be symmetric")
        object.__setattr__(self, "rows", rows)

    def score(self, left: str, right: str) -> int:
        try:
            return self.rows[self.alphabet.index(left)][self.alphabet.index(right)]
        except (ValueError, IndexError) as error:
            raise NearSelfError(
                f"Substitution matrix has no score for {left!r}/{right!r}"
            ) from error

    @property
    def sha256(self) -> str:
        payload = json.dumps(
            {"alphabet": self.alphabet, "rows": self.rows},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
        return hashlib.sha256(payload).hexdigest()


BLOSUM62 = AminoAcidSubstitutionMatrix(
    name="BLOSUM62",
    version="Henikoff-Henikoff-1992-half-bit",
    alphabet=BLOSUM62_ALPHABET,
    rows=BLOSUM62_ROWS,
    citation="Henikoff S, Henikoff JG. PNAS 1992;89:10915-10919.",
    source_url="https://doi.org/10.1073/pnas.89.22.10915",
)


@dataclass(frozen=True)
class SubstitutionEvidence(DataclassSerializable):
    """One zero-based differing position and its matrix score."""

    position: int
    target_residue: str
    risk_residue: str
    matrix_score: int

    def __post_init__(self):
        if (
            isinstance(self.position, bool)
            or not isinstance(self.position, Integral)
            or self.position < 0
        ):
            raise ValueError("Substitution position must be a non-negative integer")
        if (
            self.target_residue not in STANDARD_AMINO_ACIDS
            or self.risk_residue not in STANDARD_AMINO_ACIDS
            or self.target_residue == self.risk_residue
        ):
            raise ValueError("Substitution evidence requires differing residues")
        object.__setattr__(self, "position", int(self.position))


@dataclass(frozen=True)
class SimilarityComparatorProvenance(DataclassSerializable):
    metric: str
    matrix_name: str
    matrix_version: str
    matrix_sha256: str
    matrix_citation: str
    matrix_source_url: str

    def __post_init__(self):
        if not self.metric or len(self.matrix_sha256) != 64:
            raise ValueError("Similarity comparator provenance is incomplete")


@dataclass(frozen=True)
class HammingRiskComparator:
    """Equal-length Hamming distance with per-difference matrix evidence."""

    substitution_matrix: AminoAcidSubstitutionMatrix = BLOSUM62
    metric: str = "hamming"

    @property
    def provenance(self) -> SimilarityComparatorProvenance:
        matrix = self.substitution_matrix
        return SimilarityComparatorProvenance(
            metric=self.metric,
            matrix_name=matrix.name,
            matrix_version=matrix.version,
            matrix_sha256=matrix.sha256,
            matrix_citation=matrix.citation,
            matrix_source_url=matrix.source_url,
        )

    def distances(self, target: str, risk_peptides: tuple[str, ...]) -> np.ndarray:
        """Vectorized exact Hamming distances for one equal-length group."""
        if any(len(peptide) != len(target) for peptide in risk_peptides):
            raise NearSelfError("Hamming distance requires equal-length peptides")
        if not risk_peptides:
            return np.array([], dtype=np.int16)
        target_array = np.frombuffer(target.encode("ascii"), dtype=np.uint8)
        risk_array = np.frombuffer(
            "".join(risk_peptides).encode("ascii"), dtype=np.uint8
        ).reshape(len(risk_peptides), len(target))
        return np.count_nonzero(risk_array != target_array, axis=1)

    def substitutions(
        self, target: str, risk_peptide: str
    ) -> tuple[SubstitutionEvidence, ...]:
        if len(target) != len(risk_peptide):
            raise NearSelfError("Hamming distance requires equal-length peptides")
        return tuple(
            SubstitutionEvidence(
                position=position,
                target_residue=target_residue,
                risk_residue=risk_residue,
                matrix_score=self.substitution_matrix.score(
                    target_residue, risk_residue
                ),
            )
            for position, (target_residue, risk_residue) in enumerate(
                zip(target, risk_peptide)
            )
            if target_residue != risk_residue
        )


DEFAULT_NEAR_SELF_COMPARATOR = HammingRiskComparator()


@dataclass(frozen=True)
class NearSelfQuery(DataclassSerializable):
    peptide: str
    allele: str
    target_id: str = ""
    source_name: str = ""
    source_offset: int = 0

    def __post_init__(self):
        try:
            validate_amino_acid_sequence(self.peptide, "Near-self query")
        except ValueError:
            raise ValueError("Near-self query requires a canonical peptide")
        if (
            isinstance(self.source_offset, bool)
            or not isinstance(self.source_offset, Integral)
            or self.source_offset < 0
        ):
            raise ValueError("Near-self source offset must be non-negative")
        object.__setattr__(self, "source_offset", int(self.source_offset))


@dataclass(frozen=True)
class NearestRiskLigandHit(DataclassSerializable):
    risk_ligand: RiskLigand
    hamming_distance: int
    substitutions: tuple[SubstitutionEvidence, ...]
    substitution_score_sum: int

    def __post_init__(self):
        substitutions = tuple(sorted(self.substitutions, key=lambda item: item.position))
        if self.hamming_distance != len(substitutions):
            raise ValueError("Hamming distance disagrees with substitutions")
        if self.substitution_score_sum != sum(
            substitution.matrix_score for substitution in substitutions
        ):
            raise ValueError("Substitution score sum disagrees with evidence")
        object.__setattr__(self, "substitutions", substitutions)


@dataclass(frozen=True)
class NearSelfAssessment(DataclassSerializable):
    query: NearSelfQuery
    normalized_allele: str
    nearest_distance: Optional[int]
    nearest_hits: tuple[NearestRiskLigandHit, ...]
    prediction_coverage_gaps: tuple[RiskCoverageGap, ...]
    reason_codes: tuple[str, ...]
    comparator: SimilarityComparatorProvenance
    risk_index_cache_identity_sha256: str

    def __post_init__(self):
        hits = tuple(sorted(self.nearest_hits, key=lambda hit: (
            hit.risk_ligand.peptide,
            tuple(
                (source.gene_id, source.source_sequence_name)
                for source in hit.risk_ligand.sources
            ),
        )))
        if hits and self.nearest_distance is None:
            raise ValueError("Near-self hits require a nearest distance")
        if not hits and self.nearest_distance is not None:
            raise ValueError("Near-self distance requires at least one hit")
        if any(hit.hamming_distance != self.nearest_distance for hit in hits):
            raise ValueError("Near-self result contains a non-nearest hit")
        if len(self.risk_index_cache_identity_sha256) != 64:
            raise ValueError("Near-self result requires risk-index provenance")
        object.__setattr__(self, "nearest_hits", hits)
        object.__setattr__(
            self,
            "prediction_coverage_gaps",
            tuple(sorted(self.prediction_coverage_gaps, key=lambda gap: (
                gap.allele, gap.peptide_length, gap.kind
            ))),
        )
        object.__setattr__(self, "reason_codes", tuple(sorted(set(self.reason_codes))))

    @property
    def has_complete_coverage(self) -> bool:
        incomplete_reasons = {
            NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX,
            NEAR_SELF_REASON_PREDICTION_COVERAGE,
            NEAR_SELF_REASON_PROTEIN_RESOLUTION,
        }
        return incomplete_reasons.isdisjoint(self.reason_codes)

    def to_report_dict(self) -> dict:
        return json.loads(json.dumps(asdict(self)))


def assess_near_self_queries(
    queries: tuple[NearSelfQuery, ...],
    risk_index: PatientHLARiskLigandIndex,
    *,
    comparator=DEFAULT_NEAR_SELF_COMPARATOR,
) -> tuple[NearSelfAssessment, ...]:
    """Find all equally nearest risk ligands for each target query.

    Comparisons are restricted to the same normalized patient allele and exact
    peptide length. Coverage gaps remain attached to results and do not suppress
    similarity evidence that can still be computed.
    """
    try:
        comparator_provenance = comparator.provenance
    except Exception as error:
        raise NearSelfError("Similarity comparator lacks provenance") from error
    if not isinstance(comparator_provenance, SimilarityComparatorProvenance):
        raise NearSelfError("Similarity comparator provenance has the wrong type")
    risk_by_group: dict[tuple[str, int], tuple[RiskLigand, ...]] = {}
    for allele in risk_index.alleles:
        for length in risk_index.peptide_lengths:
            risk_by_group[(allele, length)] = risk_index.for_allele_and_length(
                allele, length
            )

    protein_coverage = risk_index.protein_resolution_coverage
    protein_resolution_gap = bool(
        protein_coverage.unresolved_gene_ids
        or protein_coverage.invalid_sequence_transcript_count
    )
    results = []
    for query in queries:
        try:
            validate_amino_acid_sequence(query.peptide, "Target peptide")
            allele = normalize_mhc_allele(query.allele)
        except ValueError as error:
            raise NearSelfError(str(error)) from error
        group_key = (allele, len(query.peptide))
        gaps = tuple(
            gap for gap in risk_index.coverage.missing_combinations
            if gap.allele == allele and gap.peptide_length == len(query.peptide)
        )
        reasons = []
        if protein_resolution_gap:
            reasons.append(NEAR_SELF_REASON_PROTEIN_RESOLUTION)
        if gaps:
            reasons.append(NEAR_SELF_REASON_PREDICTION_COVERAGE)
        if (
            allele not in risk_index.alleles
            or len(query.peptide) not in risk_index.peptide_lengths
        ):
            reasons.append(NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX)
            risk_ligands = ()
        else:
            risk_ligands = risk_by_group[group_key]

        hits = []
        nearest_distance = None
        if risk_ligands:
            risk_peptides = tuple(ligand.peptide for ligand in risk_ligands)
            try:
                distances = np.asarray(
                    comparator.distances(query.peptide, risk_peptides)
                )
            except Exception as error:
                raise NearSelfError("Similarity comparison failed") from error
            if (
                distances.shape != (len(risk_ligands),)
                or not np.issubdtype(distances.dtype, np.integer)
                or np.any(distances < 0)
            ):
                raise NearSelfError(
                    "Similarity comparator returned invalid distances"
                )
            nearest_distance = int(np.min(distances))
            for index in np.flatnonzero(distances == nearest_distance):
                risk_ligand = risk_ligands[int(index)]
                substitutions = comparator.substitutions(
                    query.peptide, risk_ligand.peptide
                )
                expected_differences = tuple(
                    (position, target_residue, risk_residue)
                    for position, (target_residue, risk_residue) in enumerate(
                        zip(query.peptide, risk_ligand.peptide)
                    )
                    if target_residue != risk_residue
                )
                observed_differences = tuple(
                    (
                        substitution.position,
                        substitution.target_residue,
                        substitution.risk_residue,
                    )
                    for substitution in substitutions
                )
                if observed_differences != expected_differences:
                    raise NearSelfError(
                        "Similarity substitutions disagree with peptide sequences"
                    )
                hits.append(NearestRiskLigandHit(
                    risk_ligand=risk_ligand,
                    hamming_distance=nearest_distance,
                    substitutions=substitutions,
                    substitution_score_sum=sum(
                        substitution.matrix_score
                        for substitution in substitutions
                    ),
                ))
        else:
            reasons.append(NEAR_SELF_REASON_NO_RISK_LIGANDS)

        results.append(NearSelfAssessment(
            query=query,
            normalized_allele=allele,
            nearest_distance=nearest_distance,
            nearest_hits=tuple(hits),
            prediction_coverage_gaps=gaps,
            reason_codes=tuple(reasons),
            comparator=comparator_provenance,
            risk_index_cache_identity_sha256=(
                risk_index.provenance.cache_identity_sha256
            ),
        ))
    return tuple(results)


def assess_near_self(
    peptide: str,
    allele: str,
    risk_index: PatientHLARiskLigandIndex,
    *,
    comparator=DEFAULT_NEAR_SELF_COMPARATOR,
) -> NearSelfAssessment:
    """Convenience wrapper for one target peptide and patient allele."""
    return assess_near_self_queries(
        (NearSelfQuery(peptide=peptide, allele=allele),),
        risk_index,
        comparator=comparator,
    )[0]
