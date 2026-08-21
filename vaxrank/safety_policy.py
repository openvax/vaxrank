# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Public, deterministic policy decisions over immutable safety evidence.

This module decides whether already-computed safety facts are safe, risky, or
coverage-incomplete. It does not shift, trim, exclude, rank, or build vaccine
constructs; those actions consume the decisions produced here.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Optional

from serializable import DataclassSerializable

from .near_self import BLOSUM62, AminoAcidSubstitutionMatrix, NearSelfAssessment
from .safety_assessment import (
    SAFETY_REASON_NO_PREDICTIONS,
    AntigenSafetyAssessment,
    EmittedSafetyLigand,
    SafetyPrediction,
)
from .vaccine_antigen import STANDARD_AMINO_ACIDS


SAFETY_MODE_OFF = "off"
SAFETY_MODE_AUDIT = "audit"
SAFETY_MODE_ENFORCE = "enforce"
SAFETY_MODES = frozenset({SAFETY_MODE_OFF, SAFETY_MODE_AUDIT, SAFETY_MODE_ENFORCE})

SAFETY_MODEL_COMBINATION_ANY = "any"
SAFETY_MODEL_COMBINATION_ALL = "all"
SAFETY_MODEL_COMBINATIONS = frozenset({
    SAFETY_MODEL_COMBINATION_ANY,
    SAFETY_MODEL_COMBINATION_ALL,
})

SAFETY_FINDING_DISABLED = "disabled"
SAFETY_FINDING_SAFE = "safe"
SAFETY_FINDING_UNSAFE = "unsafe"
SAFETY_FINDING_INCOMPLETE = "incomplete"

SAFETY_ACTION_ALLOW = "allow"
SAFETY_ACTION_REPAIR = "repair"
SAFETY_ACTION_ERROR = "error"

SAFETY_REASON_DISABLED = "safety_disabled"
SAFETY_REASON_NON_TARGET = "threshold_non_target_ligand"
SAFETY_REASON_EXACT_SELF = "threshold_exact_self_ligand"
SAFETY_REASON_NEAR_SELF_DISTANCE = "near_self_below_minimum_distance"
SAFETY_REASON_NEAR_SELF_CONSERVATIVE = "near_self_conservative_substitution"
SAFETY_REASON_MISSING_PREDICTION_KIND = "missing_safety_prediction_kind"
SAFETY_REASON_INCOMPLETE_NEAR_SELF = "incomplete_near_self_coverage"


class SafetyPolicyError(RuntimeError):
    """Fail-closed enforcement error retaining the complete decision."""

    def __init__(self, message: str, decision=None):
        super().__init__(message)
        self.decision = decision


@dataclass(frozen=True)
class SafetyPercentileThreshold(DataclassSerializable):
    """Maximum per-model percentile rank that enters safety evaluation."""

    kind: str
    max_percentile_rank: float = 1.0

    def __post_init__(self):
        if not self.kind:
            raise ValueError("Safety percentile threshold requires a prediction kind")
        if not 0.0 <= self.max_percentile_rank <= 100.0:
            raise ValueError("Safety percentile threshold must be between 0 and 100")


@dataclass(frozen=True)
class SafetyEnforcementPolicy(DataclassSerializable):
    """Resolved policy shared by off, audit, and enforce modes."""

    mode: str = SAFETY_MODE_ENFORCE
    prediction_thresholds: tuple[SafetyPercentileThreshold, ...] = (
        SafetyPercentileThreshold("pMHC_affinity", 1.0),
        SafetyPercentileThreshold("pMHC_presentation", 1.0),
    )
    model_combination: str = SAFETY_MODEL_COMBINATION_ANY
    minimum_distance: int = 2
    conservative_minimum_distance: int = 3
    substitution_matrix: AminoAcidSubstitutionMatrix = BLOSUM62
    conservative_substitution_min_score: int = 2
    additional_conservative_pairs: tuple[tuple[str, str], ...] = field(
        default_factory=tuple
    )

    def __post_init__(self):
        if self.mode not in SAFETY_MODES:
            raise ValueError(f"Unknown safety mode {self.mode!r}")
        thresholds = tuple(sorted(
            self.prediction_thresholds,
            key=lambda threshold: threshold.kind,
        ))
        kinds = [threshold.kind for threshold in thresholds]
        if not thresholds or len(kinds) != len(set(kinds)):
            raise ValueError("Safety policy requires unique prediction kinds")
        if self.model_combination not in SAFETY_MODEL_COMBINATIONS:
            raise ValueError("Unknown safety model-combination policy")
        if (
            isinstance(self.minimum_distance, bool)
            or isinstance(self.conservative_minimum_distance, bool)
            or self.minimum_distance < 0
            or self.conservative_minimum_distance < self.minimum_distance
        ):
            raise ValueError("Safety distance thresholds are inconsistent")
        pairs = []
        for pair in self.additional_conservative_pairs:
            if (
                len(pair) != 2
                or pair[0] not in STANDARD_AMINO_ACIDS
                or pair[1] not in STANDARD_AMINO_ACIDS
                or pair[0] == pair[1]
            ):
                raise ValueError(
                    "Additional conservative pairs require two differing residues"
                )
            pairs.append(tuple(sorted(pair)))
        object.__setattr__(self, "prediction_thresholds", thresholds)
        object.__setattr__(
            self,
            "additional_conservative_pairs",
            tuple(sorted(set(pairs))),
        )

    @property
    def threshold_by_kind(self) -> dict[str, float]:
        return {
            threshold.kind: threshold.max_percentile_rank
            for threshold in self.prediction_thresholds
        }

    def threshold_predictions(
        self,
        predictions: tuple[SafetyPrediction, ...],
    ) -> tuple[SafetyPrediction, ...]:
        """Predictions whose configured per-model percentile gate is met."""
        thresholds = self.threshold_by_kind
        return tuple(
            prediction
            for prediction in predictions
            if (
                prediction.kind in thresholds
                and prediction.percentile_rank is not None
                and prediction.percentile_rank <= thresholds[prediction.kind]
            )
        )

    def missing_prediction_kinds(
        self,
        predictions: tuple[SafetyPrediction, ...],
    ) -> tuple[str, ...]:
        """Configured kinds lacking a finite percentile for this ligand/allele."""
        observed = {
            prediction.kind
            for prediction in predictions
            if (
                prediction.kind in self.threshold_by_kind
                and prediction.percentile_rank is not None
            )
        }
        return tuple(sorted(set(self.threshold_by_kind) - observed))

    def prediction_gate_is_met(
        self,
        predictions: tuple[SafetyPrediction, ...],
    ) -> bool:
        """Whether the configured union/intersection prediction gate is met."""
        passing_kinds = {
            prediction.kind
            for prediction in self.threshold_predictions(predictions)
        }
        if self.model_combination == SAFETY_MODEL_COMBINATION_ANY:
            return bool(passing_kinds)
        return passing_kinds == set(self.threshold_by_kind)

    def near_self_risk_reason_codes(
        self,
        assessment: NearSelfAssessment,
    ) -> tuple[str, ...]:
        """Apply the minimum-distance and conservative-substitution gate."""
        if assessment.comparator.matrix_sha256 != self.substitution_matrix.sha256:
            raise SafetyPolicyError(
                "Near-self evidence uses a different substitution matrix"
            )
        distance = assessment.nearest_distance
        if distance is None:
            return ()
        if distance < self.minimum_distance:
            return (SAFETY_REASON_NEAR_SELF_DISTANCE,)
        if distance >= self.conservative_minimum_distance:
            return ()
        additional_pairs = set(self.additional_conservative_pairs)
        has_conservative_substitution = any(
            (
                substitution.matrix_score
                >= self.conservative_substitution_min_score
                or tuple(sorted((
                    substitution.target_residue,
                    substitution.risk_residue,
                ))) in additional_pairs
            )
            for hit in assessment.nearest_hits
            for substitution in hit.substitutions
        )
        return (
            (SAFETY_REASON_NEAR_SELF_CONSERVATIVE,)
            if has_conservative_substitution
            else ()
        )


DEFAULT_SAFETY_ENFORCEMENT_POLICY = SafetyEnforcementPolicy()


@dataclass(frozen=True)
class SafetyLigandDecision(DataclassSerializable):
    """Policy result for one emitted ligand and patient allele."""

    ligand: EmittedSafetyLigand
    allele: str
    threshold_predictions: tuple[SafetyPrediction, ...]
    near_self_assessment: Optional[NearSelfAssessment]
    missing_prediction_kinds: tuple[str, ...]
    risk_reason_codes: tuple[str, ...]
    coverage_reason_codes: tuple[str, ...]

    def __post_init__(self):
        predictions = tuple(sorted(
            self.threshold_predictions,
            key=lambda prediction: prediction.identity,
        ))
        if any(prediction.allele != self.allele for prediction in predictions):
            raise ValueError("Ligand decision prediction allele mismatch")
        if any(prediction not in self.ligand.predictions for prediction in predictions):
            raise ValueError("Ligand decision contains foreign prediction evidence")
        if self.near_self_assessment is not None and (
            self.near_self_assessment.normalized_allele != self.allele
            or self.near_self_assessment.query.peptide != self.ligand.peptide
            or self.near_self_assessment.query.source_offset
            != self.ligand.antigen_start_offset
        ):
            raise ValueError("Ligand decision near-self evidence mismatch")
        object.__setattr__(self, "threshold_predictions", predictions)
        object.__setattr__(
            self,
            "missing_prediction_kinds",
            tuple(sorted(set(self.missing_prediction_kinds))),
        )
        object.__setattr__(
            self,
            "risk_reason_codes",
            tuple(sorted(set(self.risk_reason_codes))),
        )
        object.__setattr__(
            self,
            "coverage_reason_codes",
            tuple(sorted(set(self.coverage_reason_codes))),
        )

    @property
    def has_safety_risk(self) -> bool:
        return bool(self.risk_reason_codes)

    @property
    def has_complete_coverage(self) -> bool:
        return not self.coverage_reason_codes


@dataclass(frozen=True)
class AntigenSafetyPolicyDecision(DataclassSerializable):
    """Policy finding and non-mutating action for one antigen window."""

    assessment: AntigenSafetyAssessment
    policy: SafetyEnforcementPolicy
    ligand_decisions: tuple[SafetyLigandDecision, ...]
    finding: str
    action: str
    reason_codes: tuple[str, ...]

    def __post_init__(self):
        decisions = tuple(sorted(
            self.ligand_decisions,
            key=lambda decision: (
                decision.ligand.antigen_start_offset,
                decision.ligand.peptide,
                decision.allele,
            ),
        ))
        identities = [
            (
                decision.ligand.peptide,
                decision.ligand.antigen_start_offset,
                decision.allele,
            )
            for decision in decisions
        ]
        if len(identities) != len(set(identities)):
            raise ValueError("Safety policy decision has duplicate ligand/allele rows")
        expected_findings = {
            SAFETY_FINDING_DISABLED,
            SAFETY_FINDING_SAFE,
            SAFETY_FINDING_UNSAFE,
            SAFETY_FINDING_INCOMPLETE,
        }
        if self.finding not in expected_findings:
            raise ValueError("Unknown antigen safety finding")
        if self.action not in {
            SAFETY_ACTION_ALLOW,
            SAFETY_ACTION_REPAIR,
            SAFETY_ACTION_ERROR,
        }:
            raise ValueError("Unknown antigen safety action")
        if self.policy.mode == SAFETY_MODE_OFF:
            if (
                self.finding != SAFETY_FINDING_DISABLED
                or self.action != SAFETY_ACTION_ALLOW
                or decisions
                or SAFETY_REASON_DISABLED not in self.reason_codes
            ):
                raise ValueError("Off mode requires one disabled, non-mutating result")
        else:
            if self.finding == SAFETY_FINDING_DISABLED:
                raise ValueError("Audit/enforce modes cannot report safety disabled")
            expected_action = (
                SAFETY_ACTION_ERROR
                if (
                    self.policy.mode == SAFETY_MODE_ENFORCE
                    and self.finding == SAFETY_FINDING_INCOMPLETE
                )
                else SAFETY_ACTION_REPAIR
                if (
                    self.policy.mode == SAFETY_MODE_ENFORCE
                    and self.finding == SAFETY_FINDING_UNSAFE
                )
                else SAFETY_ACTION_ALLOW
            )
            if self.action != expected_action:
                raise ValueError("Safety finding, mode, and action are inconsistent")
        object.__setattr__(self, "ligand_decisions", decisions)
        object.__setattr__(self, "reason_codes", tuple(sorted(set(self.reason_codes))))

    @property
    def would_require_repair(self) -> bool:
        return any(decision.has_safety_risk for decision in self.ligand_decisions)

    @property
    def has_complete_coverage(self) -> bool:
        return self.policy.mode != SAFETY_MODE_OFF and all(
            decision.has_complete_coverage for decision in self.ligand_decisions
        ) and SAFETY_REASON_NO_PREDICTIONS not in self.reason_codes

    def to_report_dict(self) -> dict:
        """JSON-native decision tree for audit and counterfactual reports."""
        return json.loads(json.dumps(asdict(self)))


def evaluate_antigen_safety_policy(
    assessment: AntigenSafetyAssessment,
    policy: SafetyEnforcementPolicy = DEFAULT_SAFETY_ENFORCEMENT_POLICY,
) -> AntigenSafetyPolicyDecision:
    """Evaluate immutable safety facts without repairing or excluding content."""
    if policy.mode == SAFETY_MODE_OFF:
        return AntigenSafetyPolicyDecision(
            assessment=assessment,
            policy=policy,
            ligand_decisions=(),
            finding=SAFETY_FINDING_DISABLED,
            action=SAFETY_ACTION_ALLOW,
            reason_codes=(SAFETY_REASON_DISABLED,),
        )

    near_self_by_identity = {
        (
            item.query.peptide,
            item.query.source_offset,
            item.normalized_allele,
        ): item
        for item in assessment.near_self_assessments
    }
    ligand_decisions = []
    for ligand in assessment.window_assessment.ligands:
        for allele in sorted({
            prediction.allele for prediction in ligand.predictions
        }):
            predictions = tuple(
                prediction
                for prediction in ligand.predictions
                if prediction.allele == allele
            )
            threshold_predictions = policy.threshold_predictions(predictions)
            missing_kinds = policy.missing_prediction_kinds(predictions)
            coverage_reasons = (
                (SAFETY_REASON_MISSING_PREDICTION_KIND,)
                if missing_kinds
                else ()
            )
            near_self_assessment = near_self_by_identity.get((
                ligand.peptide,
                ligand.antigen_start_offset,
                allele,
            ))
            risk_reasons = []
            if policy.prediction_gate_is_met(predictions):
                if ligand.is_non_target:
                    risk_reasons.append(SAFETY_REASON_NON_TARGET)
                if ligand.occurs_in_self_reference:
                    risk_reasons.append(SAFETY_REASON_EXACT_SELF)
                elif ligand.overlaps_targetable and near_self_assessment is not None:
                    risk_reasons.extend(
                        policy.near_self_risk_reason_codes(near_self_assessment)
                    )
                    if not near_self_assessment.has_complete_coverage:
                        coverage_reasons = (
                            *coverage_reasons,
                            SAFETY_REASON_INCOMPLETE_NEAR_SELF,
                        )
            ligand_decisions.append(SafetyLigandDecision(
                ligand=ligand,
                allele=allele,
                threshold_predictions=threshold_predictions,
                near_self_assessment=near_self_assessment,
                missing_prediction_kinds=missing_kinds,
                risk_reason_codes=tuple(risk_reasons),
                coverage_reason_codes=coverage_reasons,
            ))

    reason_codes = {
        code
        for decision in ligand_decisions
        for code in (*decision.risk_reason_codes, *decision.coverage_reason_codes)
    }
    if SAFETY_REASON_NO_PREDICTIONS in assessment.window_assessment.reason_codes:
        reason_codes.add(SAFETY_REASON_NO_PREDICTIONS)
    incomplete = any(
        not decision.has_complete_coverage for decision in ligand_decisions
    ) or SAFETY_REASON_NO_PREDICTIONS in reason_codes
    unsafe = any(decision.has_safety_risk for decision in ligand_decisions)
    finding = (
        SAFETY_FINDING_INCOMPLETE
        if incomplete
        else SAFETY_FINDING_UNSAFE
        if unsafe
        else SAFETY_FINDING_SAFE
    )
    action = (
        SAFETY_ACTION_ERROR
        if policy.mode == SAFETY_MODE_ENFORCE and incomplete
        else SAFETY_ACTION_REPAIR
        if policy.mode == SAFETY_MODE_ENFORCE and unsafe
        else SAFETY_ACTION_ALLOW
    )
    decision = AntigenSafetyPolicyDecision(
        assessment=assessment,
        policy=policy,
        ligand_decisions=tuple(ligand_decisions),
        finding=finding,
        action=action,
        reason_codes=tuple(reason_codes),
    )
    if action == SAFETY_ACTION_ERROR:
        raise SafetyPolicyError(
            "Safety enforcement requires complete prediction and near-self coverage",
            decision=decision,
        )
    return decision
