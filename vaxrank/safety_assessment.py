# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Source-agnostic safety observations for complete antigen windows.

This module deliberately does not rank, filter, repair, or exclude vaccine
content.  It records every prediction row emitted for a concrete amino-acid
window so later safety policy can act on complete evidence rather than the
subset that passed an immunogenicity filter.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from numbers import Integral
from typing import Any, Optional

from serializable import DataclassSerializable
from topiary import TopiaryPredictor

from .reference_proteome import ReferenceProteome
from .vaccine_antigen import SelfReferenceMatch, VaccineAntigen


SAFETY_REASON_NO_PREDICTIONS = "no_predictions_emitted"


class SafetyAssessmentError(RuntimeError):
    """Raised when a safety inventory cannot be produced unambiguously."""


def _finite_or_none(value):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise SafetyAssessmentError(
            f"Prediction value {value!r} is not numeric"
        ) from error
    return value if math.isfinite(value) else None


def _required_int(value, label: str) -> int:
    if isinstance(value, bool):
        raise SafetyAssessmentError(f"{label} must be an integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise SafetyAssessmentError(f"{label} must be an integer") from error
    try:
        equivalent = float(value) == result
    except (TypeError, ValueError):
        equivalent = False
    if not equivalent:
        raise SafetyAssessmentError(f"{label} must be an integer")
    return result


@dataclass(frozen=True)
class ConstructBoundary(DataclassSerializable):
    """A junction within the scanned window, using a half-open offset."""

    offset: int
    left_source: str = ""
    right_source: str = ""

    def __post_init__(self):
        if (
            isinstance(self.offset, bool)
            or not isinstance(self.offset, Integral)
            or self.offset <= 0
        ):
            raise ValueError("Construct boundary offset must be a positive integer")
        object.__setattr__(self, "offset", int(self.offset))


@dataclass(frozen=True)
class SafetyPrediction(DataclassSerializable):
    """One model/allele prediction retained exactly as safety evidence."""

    kind: str
    predictor_name: str
    predictor_version: str
    allele: str
    score: Optional[float] = None
    value: Optional[float] = None
    percentile_rank: Optional[float] = None

    def __post_init__(self):
        if not self.kind:
            raise ValueError("Safety prediction kind is required")
        for label in ("score", "value", "percentile_rank"):
            value = getattr(self, label)
            if value is not None and not math.isfinite(value):
                raise ValueError(f"Safety prediction {label} must be finite or None")

    @property
    def identity(self) -> tuple[str, str, str, str]:
        return (
            self.kind,
            self.predictor_name,
            self.predictor_version,
            self.allele,
        )


@dataclass(frozen=True)
class SafetyPredictionCoverage(DataclassSerializable):
    """Observed output coverage for one concrete predictor/version/kind."""

    kind: str
    predictor_name: str
    predictor_version: str
    alleles: tuple[str, ...]
    peptide_lengths: tuple[int, ...]
    prediction_count: int

    def __post_init__(self):
        if not self.kind:
            raise ValueError("Prediction coverage kind is required")
        object.__setattr__(self, "alleles", tuple(sorted(set(self.alleles))))
        object.__setattr__(
            self, "peptide_lengths", tuple(sorted(set(self.peptide_lengths)))
        )
        if not self.peptide_lengths or any(
            isinstance(length, bool) or length <= 0
            for length in self.peptide_lengths
        ):
            raise ValueError("Prediction coverage peptide lengths must be positive")
        if self.prediction_count <= 0:
            raise ValueError("Prediction coverage requires at least one record")


@dataclass(frozen=True)
class EmittedSafetyLigand(DataclassSerializable):
    """All safety facts for one peptide position emitted by a predictor."""

    peptide: str
    window_start_offset: int
    window_end_offset: int
    antigen_start_offset: int
    antigen_end_offset: int
    overlaps_targetable: bool
    self_reference_match: SelfReferenceMatch
    predictions: tuple[SafetyPrediction, ...]
    crosses_antigen_boundary: bool = False
    crosses_window_boundary: bool = False
    crossed_construct_boundaries: tuple[ConstructBoundary, ...] = field(
        default_factory=tuple
    )

    def __post_init__(self):
        if not self.peptide:
            raise ValueError("Safety ligand peptide is required")
        if self.window_start_offset < 0:
            raise ValueError("Safety ligand window offset cannot be negative")
        if self.window_end_offset - self.window_start_offset != len(self.peptide):
            raise ValueError("Safety ligand window coordinates do not match peptide")
        if self.antigen_end_offset - self.antigen_start_offset != len(self.peptide):
            raise ValueError("Safety ligand antigen coordinates do not match peptide")
        if self.self_reference_match.peptide != self.peptide:
            raise ValueError("Safety ligand and self-reference peptide differ")
        if not self.predictions:
            raise ValueError("Safety ligand requires prediction evidence")
        predictions = tuple(sorted(self.predictions, key=lambda prediction: (
            prediction.identity,
            prediction.score is None,
            prediction.score or 0.0,
            prediction.value is None,
            prediction.value or 0.0,
            prediction.percentile_rank is None,
            prediction.percentile_rank or 0.0,
        )))
        identities = [prediction.identity for prediction in predictions]
        if len(identities) != len(set(identities)):
            raise ValueError(
                "Safety ligand has duplicate predictor/kind/allele evidence"
            )
        object.__setattr__(self, "predictions", predictions)
        object.__setattr__(
            self,
            "crossed_construct_boundaries",
            tuple(sorted(
                self.crossed_construct_boundaries,
                key=lambda boundary: (
                    boundary.offset,
                    boundary.left_source,
                    boundary.right_source,
                ),
            )),
        )

    @property
    def occurs_in_self_reference(self) -> bool:
        return self.self_reference_match.occurs

    @property
    def is_non_target(self) -> bool:
        return not self.overlaps_targetable

    @property
    def crosses_construct_boundary(self) -> bool:
        return bool(self.crossed_construct_boundaries)


@dataclass(frozen=True)
class WindowSafetyAssessment(DataclassSerializable):
    """Immutable complete prediction inventory for one antigen window."""

    antigen: VaccineAntigen
    window_start_offset: int
    window_end_offset: int
    ligands: tuple[EmittedSafetyLigand, ...]
    coverage: tuple[SafetyPredictionCoverage, ...]
    genome_release: str = ""
    reason_codes: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self):
        if self.window_start_offset < 0:
            raise ValueError("Window start cannot be negative")
        if not self.window_start_offset < self.window_end_offset:
            raise ValueError("Window must have positive length")
        if self.window_end_offset > len(self.antigen.amino_acids):
            raise ValueError("Window extends beyond vaccine antigen")
        ligands = tuple(sorted(
            self.ligands,
            key=lambda ligand: (
                ligand.window_start_offset,
                ligand.window_end_offset,
                ligand.peptide,
            ),
        ))
        positions = [
            (ligand.peptide, ligand.window_start_offset)
            for ligand in ligands
        ]
        if len(positions) != len(set(positions)):
            raise ValueError("Safety assessment contains duplicate ligand positions")
        window_sequence = self.antigen.amino_acids[
            self.window_start_offset:self.window_end_offset
        ]
        for ligand in ligands:
            expected_antigen_start = (
                self.window_start_offset + ligand.window_start_offset
            )
            if ligand.antigen_start_offset != expected_antigen_start:
                raise ValueError(
                    "Safety ligand window and antigen coordinates disagree"
                )
            observed = window_sequence[
                ligand.window_start_offset:ligand.window_end_offset
            ]
            if observed != ligand.peptide:
                raise ValueError(
                    "Safety ligand does not match the scanned antigen window"
                )
        object.__setattr__(self, "ligands", ligands)
        coverage = tuple(sorted(
            self.coverage,
            key=lambda item: (
                item.kind,
                item.predictor_name,
                item.predictor_version,
            ),
        ))
        coverage_keys = [
            (item.kind, item.predictor_name, item.predictor_version)
            for item in coverage
        ]
        if len(coverage_keys) != len(set(coverage_keys)):
            raise ValueError("Safety assessment has duplicate coverage records")
        expected_coverage: dict[tuple[str, str, str], dict[str, Any]] = {}
        for ligand in ligands:
            for prediction in ligand.predictions:
                key = (
                    prediction.kind,
                    prediction.predictor_name,
                    prediction.predictor_version,
                )
                expected = expected_coverage.setdefault(
                    key,
                    {"alleles": set(), "lengths": set(), "count": 0},
                )
                expected["alleles"].add(prediction.allele)
                expected["lengths"].add(len(ligand.peptide))
                expected["count"] += 1
        for item in coverage:
            expected = expected_coverage.pop(
                (item.kind, item.predictor_name, item.predictor_version),
                None,
            )
            if expected is None or (
                set(item.alleles) != expected["alleles"]
                or set(item.peptide_lengths) != expected["lengths"]
                or item.prediction_count != expected["count"]
            ):
                raise ValueError(
                    "Safety prediction coverage disagrees with ligand evidence"
                )
        if expected_coverage:
            raise ValueError("Safety prediction coverage omits ligand evidence")
        object.__setattr__(self, "coverage", coverage)
        object.__setattr__(self, "reason_codes", tuple(sorted(set(self.reason_codes))))
        if self.ligands and SAFETY_REASON_NO_PREDICTIONS in self.reason_codes:
            raise ValueError("Non-empty assessment cannot claim no predictions")
        if not self.ligands and SAFETY_REASON_NO_PREDICTIONS not in self.reason_codes:
            raise ValueError("Empty assessment must explicitly record no predictions")

    @property
    def window_sequence(self) -> str:
        return self.antigen.amino_acids[
            self.window_start_offset:self.window_end_offset
        ]

    def to_report_dict(self) -> dict:
        """Plain JSON-compatible tree without serializer type metadata."""
        return _json_native(asdict(self))

    def prediction_rows(self) -> list[dict]:
        """Flat, primitive-valued rows for CSV/DataFrame report inputs."""
        rows = []
        for ligand in self.ligands:
            base = {
                "antigen_kind": self.antigen.kind,
                "antigen_source_identifier": self.antigen.source_identifier,
                "genome_release": self.genome_release,
                "window_start_offset": self.window_start_offset,
                "window_end_offset": self.window_end_offset,
                "peptide": ligand.peptide,
                "peptide_window_start_offset": ligand.window_start_offset,
                "peptide_window_end_offset": ligand.window_end_offset,
                "peptide_antigen_start_offset": ligand.antigen_start_offset,
                "peptide_antigen_end_offset": ligand.antigen_end_offset,
                "overlaps_targetable": ligand.overlaps_targetable,
                "occurs_in_self_reference": ligand.occurs_in_self_reference,
                "self_reference_excluded_gene_ids": ";".join(
                    ligand.self_reference_match.excluded_gene_ids
                ),
                "crosses_antigen_boundary": ligand.crosses_antigen_boundary,
                "crosses_window_boundary": ligand.crosses_window_boundary,
                "crosses_construct_boundary": ligand.crosses_construct_boundary,
                "crossed_construct_boundary_offsets": ";".join(
                    str(boundary.offset)
                    for boundary in ligand.crossed_construct_boundaries
                ),
            }
            for prediction in ligand.predictions:
                rows.append(base | {
                    "prediction_kind": prediction.kind,
                    "predictor_name": prediction.predictor_name,
                    "predictor_version": prediction.predictor_version,
                    "allele": prediction.allele,
                    "score": prediction.score,
                    "value": prediction.value,
                    "percentile_rank": prediction.percentile_rank,
                })
        return rows


def _json_native(value):
    if isinstance(value, dict):
        return {key: _json_native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_native(item) for item in value]
    return value


def _safety_prediction_from_row(row: Any) -> SafetyPrediction:
    value = _finite_or_none(row.get("affinity"))
    if value is None:
        value = _finite_or_none(row.get("value"))
    return SafetyPrediction(
        kind=str(row.get("kind") or "pMHC_affinity"),
        predictor_name=str(row.get("prediction_method_name") or ""),
        predictor_version=str(row.get("predictor_version") or ""),
        allele=str(row.get("allele") or ""),
        score=_finite_or_none(row.get("score")),
        value=value,
        percentile_rank=_finite_or_none(row.get("percentile_rank")),
    )


def safety_assessment_from_prediction_frame(
    predictions_df,
    *,
    antigen: VaccineAntigen,
    reference_proteome: ReferenceProteome,
    window_start: int = 0,
    window_end: Optional[int] = None,
    construct_boundaries: tuple[ConstructBoundary, ...] = (),
    genome_release: str = "",
    expected_source_name: Optional[str] = None,
) -> WindowSafetyAssessment:
    """Convert an unfiltered topiary prediction frame into safety evidence."""
    if window_end is None:
        window_end = len(antigen.amino_acids)
    if window_start < 0 or window_end > len(antigen.amino_acids):
        raise ValueError("Safety window lies outside the vaccine antigen")
    if window_start >= window_end:
        raise ValueError("Safety window must have positive length")
    window_sequence = antigen.amino_acids[window_start:window_end]
    boundaries = tuple(construct_boundaries)
    for boundary in boundaries:
        if boundary.offset >= len(window_sequence):
            raise ValueError("Construct boundary lies outside the safety window")

    if predictions_df is None or predictions_df.empty:
        return WindowSafetyAssessment(
            antigen=antigen,
            window_start_offset=window_start,
            window_end_offset=window_end,
            ligands=(),
            coverage=(),
            genome_release=genome_release,
            reason_codes=(SAFETY_REASON_NO_PREDICTIONS,),
        )

    groups: dict[tuple[str, int], dict[str, Any]] = {}
    coverage_records: dict[tuple[str, str, str], dict[str, Any]] = {}
    for _, row in predictions_df.iterrows():
        if expected_source_name is not None:
            observed_source_name = row.get("source_sequence_name")
            if observed_source_name != expected_source_name:
                raise SafetyAssessmentError(
                    "Prediction source does not match the scanned window"
                )
        peptide = str(row.get("peptide") or "")
        offset = _required_int(row.get("peptide_offset"), "Peptide offset")
        peptide_length = _required_int(
            row.get("peptide_length"), "Peptide length"
        )
        if peptide_length != len(peptide):
            raise SafetyAssessmentError(
                "Prediction peptide length does not match its sequence"
            )
        end = offset + peptide_length
        if offset < 0 or end > len(window_sequence):
            raise SafetyAssessmentError(
                "Prediction coordinates lie outside the scanned window"
            )
        if window_sequence[offset:end] != peptide:
            raise SafetyAssessmentError(
                "Prediction peptide does not match the scanned window"
            )

        prediction = _safety_prediction_from_row(row)
        key = (peptide, offset)
        group = groups.setdefault(key, {"predictions": [], "identities": set()})
        if prediction.identity in group["identities"]:
            raise SafetyAssessmentError(
                "Duplicate predictor/kind/allele evidence for one ligand"
            )
        group["identities"].add(prediction.identity)
        group["predictions"].append(prediction)

        coverage_key = (
            prediction.kind,
            prediction.predictor_name,
            prediction.predictor_version,
        )
        coverage = coverage_records.setdefault(
            coverage_key,
            {"alleles": set(), "peptide_lengths": set(), "count": 0},
        )
        coverage["alleles"].add(prediction.allele)
        coverage["peptide_lengths"].add(peptide_length)
        coverage["count"] += 1

    ligands = []
    for (peptide, offset), group in groups.items():
        end = offset + len(peptide)
        antigen_start = window_start + offset
        antigen_end = window_start + end
        crossed_boundaries = tuple(
            boundary
            for boundary in boundaries
            if offset < boundary.offset < end
        )
        occurs = reference_proteome.contains(peptide)
        ligands.append(EmittedSafetyLigand(
            peptide=peptide,
            window_start_offset=offset,
            window_end_offset=end,
            antigen_start_offset=antigen_start,
            antigen_end_offset=antigen_end,
            overlaps_targetable=antigen.interval_is_targetable(
                antigen_start, antigen_end
            ),
            self_reference_match=antigen.self_reference_match(
                peptide,
                occurs,
                genome_release=genome_release,
            ),
            predictions=tuple(group["predictions"]),
            crossed_construct_boundaries=crossed_boundaries,
        ))

    coverage = tuple(
        SafetyPredictionCoverage(
            kind=key[0],
            predictor_name=key[1],
            predictor_version=key[2],
            alleles=tuple(record["alleles"]),
            peptide_lengths=tuple(record["peptide_lengths"]),
            prediction_count=record["count"],
        )
        for key, record in coverage_records.items()
    )
    return WindowSafetyAssessment(
        antigen=antigen,
        window_start_offset=window_start,
        window_end_offset=window_end,
        ligands=tuple(ligands),
        coverage=coverage,
        genome_release=genome_release,
    )


def assess_vaccine_antigen_window(
    mhc_predictor,
    antigen: VaccineAntigen,
    *,
    genome=None,
    reference_proteome: Optional[ReferenceProteome] = None,
    window_start: int = 0,
    window_end: Optional[int] = None,
    construct_boundaries: tuple[ConstructBoundary, ...] = (),
) -> WindowSafetyAssessment:
    """Predict and inventory every ligand emitted for an antigen window.

    Unlike :func:`vaxrank.predict_epitopes`, this API applies no ranking or
    immunogenicity filter. Predictor failures raise ``SafetyAssessmentError``;
    they are never converted into a reassuring empty result.
    """
    if genome is not None and reference_proteome is not None:
        raise ValueError("Pass genome or reference_proteome, not both")
    if reference_proteome is None:
        if genome is None:
            raise ValueError(
                "Safety assessment requires a genome or reference proteome"
            )
        if antigen.self_reference_excluded_gene_ids:
            reference_proteome = ReferenceProteome.from_genome(
                genome,
                exclude_gene_ids=antigen.self_reference_excluded_gene_ids,
            )
        else:
            reference_proteome = ReferenceProteome(genome)
    else:
        expected_exclusions = frozenset(antigen.self_reference_excluded_gene_ids)
        actual_exclusions = frozenset(
            getattr(reference_proteome, "excluded_gene_ids", ())
        )
        if actual_exclusions != expected_exclusions:
            raise ValueError(
                "Reference proteome exclusions do not match antigen policy"
            )

    if window_end is None:
        window_end = len(antigen.amino_acids)
    if window_start < 0 or window_end > len(antigen.amino_acids):
        raise ValueError("Safety window lies outside the vaccine antigen")
    if window_start >= window_end:
        raise ValueError("Safety window must have positive length")
    sequence = antigen.amino_acids[window_start:window_end]
    source_name = antigen.source_identifier or antigen.gene_name or "antigen"
    predictor = (
        mhc_predictor
        if isinstance(mhc_predictor, TopiaryPredictor)
        else TopiaryPredictor(models=[mhc_predictor])
    )
    try:
        predictions_df = predictor.predict_from_named_sequences({
            source_name: sequence,
        })
    except Exception as error:
        raise SafetyAssessmentError(
            f"MHC safety prediction failed for {source_name!r}"
        ) from error

    genome_release = str(getattr(genome, "release", "") or "")
    return safety_assessment_from_prediction_frame(
        predictions_df,
        antigen=antigen,
        reference_proteome=reference_proteome,
        window_start=window_start,
        window_end=window_end,
        construct_boundaries=construct_boundaries,
        genome_release=genome_release,
        expected_source_name=source_name,
    )
