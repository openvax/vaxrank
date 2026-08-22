# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Direct-oncoref CTA admission and antigen construction."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from typing import Optional

import pandas as pd
from serializable import DataclassSerializable

from .identifiers import normalize_ensembl_gene_id
from .vaccine_antigen import (
    ANTIGEN_KIND_CTA,
    ATTESTATION_ADMITTED,
    ATTESTATION_HELD_OUT,
    ATTESTATION_OVERRIDDEN,
    AminoAcidInterval,
    TargetableMask,
    TumorSpecificityAttestation,
    TumorSpecificityEvidence,
    VaccineAntigen,
)


CTA_REASON_CANONICAL_EXPRESSION_PASS = "cta_canonical_expression_pass"
CTA_REASON_EXPRESSION_BELOW_THRESHOLD = "cta_expression_below_threshold"
CTA_REASON_NONCANONICAL_HELD_OUT = "cta_noncanonical_held_out"
CTA_REASON_EXPLICIT_OVERRIDE = "cta_explicit_evidence_override"


class CTAAdmissionError(RuntimeError):
    """Raised when CTA admission evidence cannot be resolved unambiguously."""


def _finite_nonnegative(value, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{label} must be finite and non-negative")
    return result


def _clean(value) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value)


def _required_bool(value, label: str) -> bool:
    if isinstance(value, bool):
        return value
    if value in (0, 1):
        return bool(value)
    text = _clean(value).strip().casefold()
    if text == "true":
        return True
    if text == "false":
        return False
    raise CTAAdmissionError(f"oncoref CTA {label} is not boolean")


def _fingerprint(value) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class CTAAdmissionPolicy(DataclassSerializable):
    """Run-specific quantitative tumor-expression admission rule."""

    min_tumor_expression: float
    expression_unit: str

    def __post_init__(self):
        object.__setattr__(
            self,
            "min_tumor_expression",
            _finite_nonnegative(
                self.min_tumor_expression, "Minimum tumor expression"
            ),
        )
        if self.min_tumor_expression <= 0:
            raise ValueError("Minimum tumor expression must be greater than zero")
        if not self.expression_unit:
            raise ValueError("CTA expression policy requires an explicit unit")


@dataclass(frozen=True)
class PatientTumorExpressionEvidence(DataclassSerializable):
    gene_id: str
    sample_id: str
    value: float
    unit: str
    evidence_source: str
    evidence_version: str
    assay: str = ""

    def __post_init__(self):
        gene_id = normalize_ensembl_gene_id(self.gene_id)
        if not gene_id or not self.sample_id:
            raise ValueError("Tumor expression requires gene and patient sample IDs")
        if (
            not self.unit
            or not self.evidence_source
            or not self.evidence_version
            or not self.assay
        ):
            raise ValueError(
                "Tumor expression requires assay, units, and versioned provenance"
            )
        object.__setattr__(self, "gene_id", gene_id)
        object.__setattr__(
            self, "value", _finite_nonnegative(self.value, "Tumor expression")
        )


@dataclass(frozen=True)
class CTAOverrideEvidence(DataclassSerializable):
    evidence_source: str
    evidence_version: str
    rationale: str
    evidence_id: str = ""

    def __post_init__(self):
        if not self.evidence_source or not self.evidence_version or not self.rationale:
            raise ValueError("CTA override requires versioned evidence and rationale")


@dataclass(frozen=True)
class CTAReferenceEvidence(DataclassSerializable):
    gene_id: str
    gene_name: str
    canonical_default: bool
    unfiltered_candidate: bool
    passes_filters: bool
    never_expressed: bool
    specificity_status: str
    specificity_action: str
    specificity_source_anchor: str
    specificity_rationale: str
    restriction: str
    restriction_confidence: str
    safety_flags: str
    source_databases: str
    canonical_transcript_id: str
    oncoref_version: str
    oncoref_data_version: str
    oncoref_source_matrix_version: str
    canonical_gene_ids_sha256: str
    unfiltered_gene_ids_sha256: str
    source_row_sha256: str

    def __post_init__(self):
        object.__setattr__(
            self, "gene_id", normalize_ensembl_gene_id(self.gene_id)
        )
        for digest in (
            self.canonical_gene_ids_sha256,
            self.unfiltered_gene_ids_sha256,
            self.source_row_sha256,
        ):
            if len(digest) != 64:
                raise ValueError("CTA reference evidence requires SHA-256 provenance")


@dataclass(frozen=True)
class CTAAdmissionAssessment(DataclassSerializable):
    antigen: VaccineAntigen
    policy: CTAAdmissionPolicy
    tumor_expression: PatientTumorExpressionEvidence
    reference_evidence: CTAReferenceEvidence
    override_evidence: Optional[CTAOverrideEvidence] = None

    def __post_init__(self):
        if self.antigen.kind != ANTIGEN_KIND_CTA:
            raise ValueError("CTA assessment requires a CTA vaccine antigen")
        if self.antigen.gene_id != self.reference_evidence.gene_id:
            raise ValueError("CTA antigen and reference evidence genes differ")

    def to_report_dict(self) -> dict:
        return asdict(self)


def _resolve_reference_evidence(gene_id: str):
    import oncoref
    from oncoref.cta import cta_evidence, cta_gene_ids, cta_unfiltered_gene_ids
    from oncoref.version import DATA_VERSION, SOURCE_MATRIX_VERSION

    canonical_ids = frozenset(
        normalize_ensembl_gene_id(value) for value in cta_gene_ids()
    )
    unfiltered_ids = frozenset(
        normalize_ensembl_gene_id(value)
        for value in cta_unfiltered_gene_ids()
    )
    if not canonical_ids or not unfiltered_ids:
        raise CTAAdmissionError("oncoref returned an empty CTA reference set")
    if gene_id not in unfiltered_ids:
        raise CTAAdmissionError(
            f"Gene {gene_id} is not in oncoref's CTA candidate universe"
        )
    frame = cta_evidence()
    required = {
        "Ensembl_Gene_ID",
        "Symbol",
        "passes_filters",
        "never_expressed",
        "specificity_status",
        "specificity_action",
        "specificity_source_anchor",
        "specificity_rationale",
        "restriction",
        "restriction_confidence",
        "safety_flags",
        "source_databases",
        "Canonical_Transcript_ID",
    }
    missing = required - set(frame.columns)
    if missing:
        raise CTAAdmissionError(
            f"oncoref CTA evidence is missing columns: {sorted(missing)}"
        )
    normalized_ids = frame["Ensembl_Gene_ID"].map(normalize_ensembl_gene_id)
    rows = frame.loc[normalized_ids == gene_id]
    if len(rows) != 1:
        raise CTAAdmissionError(
            f"Expected one oncoref CTA evidence row for {gene_id}; found {len(rows)}"
        )
    row = rows.iloc[0]
    row_payload = {
        column: _clean(row[column]) for column in sorted(required)
    }
    evidence = CTAReferenceEvidence(
        gene_id=gene_id,
        gene_name=_clean(row["Symbol"]),
        canonical_default=gene_id in canonical_ids,
        unfiltered_candidate=True,
        passes_filters=_required_bool(row["passes_filters"], "passes_filters"),
        never_expressed=_required_bool(
            row["never_expressed"], "never_expressed"
        ),
        specificity_status=_clean(row["specificity_status"]),
        specificity_action=_clean(row["specificity_action"]),
        specificity_source_anchor=_clean(row["specificity_source_anchor"]),
        specificity_rationale=_clean(row["specificity_rationale"]),
        restriction=_clean(row["restriction"]),
        restriction_confidence=_clean(row["restriction_confidence"]),
        safety_flags=_clean(row["safety_flags"]),
        source_databases=_clean(row["source_databases"]),
        canonical_transcript_id=_clean(row["Canonical_Transcript_ID"]),
        oncoref_version=oncoref.__version__,
        oncoref_data_version=DATA_VERSION,
        oncoref_source_matrix_version=SOURCE_MATRIX_VERSION,
        canonical_gene_ids_sha256=_fingerprint(sorted(canonical_ids)),
        unfiltered_gene_ids_sha256=_fingerprint(sorted(unfiltered_ids)),
        source_row_sha256=_fingerprint(row_payload),
    )
    return evidence, tuple(sorted(unfiltered_ids))


def assess_cta_antigen(
    *,
    amino_acids: str,
    gene_id: str,
    tumor_expression: PatientTumorExpressionEvidence,
    policy: CTAAdmissionPolicy,
    override_evidence: Optional[CTAOverrideEvidence] = None,
    transcript_ids: tuple[str, ...] = (),
    protein_ids: tuple[str, ...] = (),
    species: str = "Homo sapiens",
    source_identifier: str = "",
) -> CTAAdmissionAssessment:
    """Resolve oncoref CTA status and patient-expression construct admission."""
    gene_id = normalize_ensembl_gene_id(gene_id)
    if tumor_expression.gene_id != gene_id:
        raise ValueError("CTA and tumor-expression gene IDs differ")
    if tumor_expression.unit != policy.expression_unit:
        raise ValueError(
            "CTA tumor-expression units do not match the admission policy"
        )
    reference, self_excluded_gene_ids = _resolve_reference_evidence(gene_id)
    expression_passes = tumor_expression.value >= policy.min_tumor_expression

    if reference.canonical_default and expression_passes:
        status = ATTESTATION_ADMITTED
        rationale_code = CTA_REASON_CANONICAL_EXPRESSION_PASS
        requires_review = False
        override_reason = ""
    elif not expression_passes:
        status = ATTESTATION_HELD_OUT
        rationale_code = CTA_REASON_EXPRESSION_BELOW_THRESHOLD
        requires_review = True
        override_reason = ""
    elif override_evidence is None:
        status = ATTESTATION_HELD_OUT
        rationale_code = CTA_REASON_NONCANONICAL_HELD_OUT
        requires_review = True
        override_reason = ""
    else:
        status = ATTESTATION_OVERRIDDEN
        rationale_code = CTA_REASON_EXPLICIT_OVERRIDE
        requires_review = True
        override_reason = override_evidence.rationale

    reference_record = TumorSpecificityEvidence(
        evidence_kind="oncoref_cta_membership",
        evidence_source="oncoref.cta",
        evidence_version=reference.oncoref_version,
        subject_id=gene_id,
        passed=reference.canonical_default,
        details=(
            ("specificity_status", reference.specificity_status),
            ("specificity_action", reference.specificity_action),
            ("restriction", reference.restriction),
            ("restriction_confidence", reference.restriction_confidence),
            ("source_row_sha256", reference.source_row_sha256),
        ),
    )
    expression_record = TumorSpecificityEvidence(
        evidence_kind="patient_tumor_expression",
        evidence_source=tumor_expression.evidence_source,
        evidence_version=tumor_expression.evidence_version,
        subject_id=gene_id,
        sample_id=tumor_expression.sample_id,
        patient_specific=True,
        passed=expression_passes,
        numeric_value=tumor_expression.value,
        unit=tumor_expression.unit,
        threshold=policy.min_tumor_expression,
        details=(("assay", tumor_expression.assay),),
    )
    records = [reference_record, expression_record]
    if status == ATTESTATION_OVERRIDDEN:
        records.append(TumorSpecificityEvidence(
            evidence_kind="cta_admission_override",
            evidence_source=override_evidence.evidence_source,
            evidence_version=override_evidence.evidence_version,
            subject_id=gene_id,
            passed=True,
            details=(
                ("evidence_id", override_evidence.evidence_id),
                ("rationale", override_evidence.rationale),
            ),
        ))

    antigen = VaccineAntigen(
        kind=ANTIGEN_KIND_CTA,
        amino_acids=amino_acids,
        targetable_mask=TargetableMask((
            AminoAcidInterval(0, len(amino_acids)),
        )),
        tumor_specificity=TumorSpecificityAttestation(
            status=status,
            evidence_kind="oncoref_cta_and_patient_tumor_expression",
            evidence_source="oncoref.cta plus patient tumor expression",
            evidence_version=reference.oncoref_version,
            patient_specific=True,
            rationale_code=rationale_code,
            requires_review=requires_review,
            override_reason=override_reason,
            evidence_records=tuple(records),
        ),
        self_reference_excluded_gene_ids=self_excluded_gene_ids,
        gene_name=reference.gene_name,
        gene_id=gene_id,
        transcript_ids=transcript_ids or (
            (reference.canonical_transcript_id,)
            if reference.canonical_transcript_id else ()
        ),
        protein_ids=protein_ids,
        species=species,
        source_identifier=source_identifier or gene_id,
    )
    return CTAAdmissionAssessment(
        antigen=antigen,
        policy=policy,
        tumor_expression=tumor_expression,
        reference_evidence=reference,
        override_evidence=override_evidence,
    )


def build_cta_vaccine_antigen(
    *,
    amino_acids: str,
    gene_id: str,
    tumor_expression: PatientTumorExpressionEvidence,
    policy: CTAAdmissionPolicy,
    override_evidence: Optional[CTAOverrideEvidence] = None,
    transcript_ids: tuple[str, ...] = (),
    protein_ids: tuple[str, ...] = (),
    species: str = "Homo sapiens",
    source_identifier: str = "",
) -> VaccineAntigen:
    """Build a CTA antigen while retaining evidence in its attestation."""
    return assess_cta_antigen(
        amino_acids=amino_acids,
        gene_id=gene_id,
        tumor_expression=tumor_expression,
        policy=policy,
        override_evidence=override_evidence,
        transcript_ids=transcript_ids,
        protein_ids=protein_ids,
        species=species,
        source_identifier=source_identifier,
    ).antigen
