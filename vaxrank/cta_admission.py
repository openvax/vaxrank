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


@dataclass(frozen=True)
class CTAAdmissionPolicy(DataclassSerializable):
    """Run-specific quantitative tumor-expression admission rule."""

    min_tumor_expression: float
    expression_unit: str

    def __post_init__(self):
        try:
            minimum = float(self.min_tumor_expression)
        except (TypeError, ValueError) as error:
            raise ValueError("Minimum tumor expression must be numeric") from error
        if not math.isfinite(minimum) or minimum < 0:
            raise ValueError(
                "Minimum tumor expression must be finite and non-negative"
            )
        if minimum <= 0:
            raise ValueError("Minimum tumor expression must be greater than zero")
        object.__setattr__(self, "min_tumor_expression", minimum)
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
        try:
            value = float(self.value)
        except (TypeError, ValueError) as error:
            raise ValueError("Tumor expression must be numeric") from error
        if not math.isfinite(value) or value < 0:
            raise ValueError("Tumor expression must be finite and non-negative")
        object.__setattr__(self, "value", value)


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
class CTAReferenceResolution(DataclassSerializable):
    """Pinned oncoref evidence and its CTA self-reference exclusions."""

    evidence: CTAReferenceEvidence
    self_reference_excluded_gene_ids: tuple[str, ...]

    def __post_init__(self):
        excluded_gene_ids = tuple(sorted({
            normalize_ensembl_gene_id(gene_id)
            for gene_id in self.self_reference_excluded_gene_ids
        }))
        if not excluded_gene_ids or self.evidence.gene_id not in excluded_gene_ids:
            raise ValueError(
                "CTA reference resolution must include its source gene"
            )
        excluded_digest = hashlib.sha256(json.dumps(
            list(excluded_gene_ids),
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")).hexdigest()
        if excluded_digest != self.evidence.unfiltered_gene_ids_sha256:
            raise ValueError(
                "CTA reference exclusions disagree with oncoref provenance"
            )
        object.__setattr__(
            self, "self_reference_excluded_gene_ids", excluded_gene_ids
        )


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


def resolve_cta_reference_evidence(gene_id: str) -> CTAReferenceResolution:
    """Resolve one CTA against canonical and unfiltered oncoref facts.

    ``cta_gene_ids()`` determines default target admission, while the broader
    ``cta_unfiltered_gene_ids()`` universe determines which CTA source genes
    are omitted from the negative self reference. Both sets and the selected
    evidence row are fingerprinted in the returned immutable result.
    """
    import oncoref
    from oncoref.cta import cta_evidence, cta_gene_ids, cta_unfiltered_gene_ids
    from oncoref.version import DATA_VERSION, SOURCE_MATRIX_VERSION

    gene_id = normalize_ensembl_gene_id(gene_id)
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
        column: "" if pd.isna(row[column]) else str(row[column])
        for column in sorted(required)
    }
    boolean_values = {}
    for column in ("passes_filters", "never_expressed"):
        value = row[column]
        if isinstance(value, bool):
            boolean_values[column] = value
        elif value in (0, 1):
            boolean_values[column] = bool(value)
        elif row_payload[column].strip().casefold() in {"true", "false"}:
            boolean_values[column] = (
                row_payload[column].strip().casefold() == "true"
            )
        else:
            raise CTAAdmissionError(
                f"oncoref CTA {column} is not boolean"
            )
    fingerprint_payloads = {
        "canonical": sorted(canonical_ids),
        "unfiltered": sorted(unfiltered_ids),
        "source_row": row_payload,
    }
    fingerprints = {
        name: hashlib.sha256(json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")).hexdigest()
        for name, payload in fingerprint_payloads.items()
    }
    evidence = CTAReferenceEvidence(
        gene_id=gene_id,
        gene_name=row_payload["Symbol"],
        canonical_default=gene_id in canonical_ids,
        unfiltered_candidate=True,
        passes_filters=boolean_values["passes_filters"],
        never_expressed=boolean_values["never_expressed"],
        specificity_status=row_payload["specificity_status"],
        specificity_action=row_payload["specificity_action"],
        specificity_source_anchor=row_payload["specificity_source_anchor"],
        specificity_rationale=row_payload["specificity_rationale"],
        restriction=row_payload["restriction"],
        restriction_confidence=row_payload["restriction_confidence"],
        safety_flags=row_payload["safety_flags"],
        source_databases=row_payload["source_databases"],
        canonical_transcript_id=row_payload["Canonical_Transcript_ID"],
        oncoref_version=oncoref.__version__,
        oncoref_data_version=DATA_VERSION,
        oncoref_source_matrix_version=SOURCE_MATRIX_VERSION,
        canonical_gene_ids_sha256=fingerprints["canonical"],
        unfiltered_gene_ids_sha256=fingerprints["unfiltered"],
        source_row_sha256=fingerprints["source_row"],
    )
    return CTAReferenceResolution(
        evidence=evidence,
        self_reference_excluded_gene_ids=tuple(unfiltered_ids),
    )


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
    reference_resolution = resolve_cta_reference_evidence(gene_id)
    reference = reference_resolution.evidence
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
        self_reference_excluded_gene_ids=(
            reference_resolution.self_reference_excluded_gene_ids
        ),
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
