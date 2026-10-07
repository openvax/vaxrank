"""Auditable expression-table admission for patient CTA vaccine targets.

This adapter selects reference proteins; RNA measurements remain gene or
transcript measurements. It does not infer protein abundance or read counts.
"""

from dataclasses import dataclass, replace
import csv
from fnmatch import fnmatchcase
import hashlib
import math
from pathlib import Path
import re

import pandas as pd
from serializable import DataclassSerializable
from topiary.rna import load_expression

from .cta_admission import (
    CTAAdmissionAssessment, CTAAdmissionPolicy, CTAReferenceResolution,
    PatientTumorExpressionEvidence,
    assess_cta_antigen, resolve_cta_reference_evidence,
)
from .native_serialization import from_native_json, to_native_json


@dataclass(frozen=True)
class CTAExpressionInput(DataclassSerializable):
    """Explicit contract for one patient expression column (no unit inference)."""

    measurement_level: str
    id_column: str
    value_column: str
    expression_unit: str
    sample_id: str
    evidence_source: str
    evidence_version: str
    assay: str

    def __post_init__(self):
        if self.measurement_level not in {"gene", "transcript"}:
            raise ValueError("Expression measurement level must be gene or transcript")
        for name, value in self.to_dict().items():
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"Expression input requires {name}")
        if self.id_column == self.value_column:
            raise ValueError("Expression identifier and value columns must differ")

    @property
    def id_namespace(self):
        return f"ensembl_{self.measurement_level}_id"


@dataclass(frozen=True)
class CTATargetExclusionPolicy(DataclassSerializable):
    """Gene-symbol eligibility only; does not change the CTA self reference."""

    gene_patterns: tuple[str, ...] = ("MAGE*",)
    gene_exceptions: tuple[str, ...] = ("MAGEA4",)

    def __post_init__(self):
        for name in ("gene_patterns", "gene_exceptions"):
            if isinstance(getattr(self, name), str):
                raise ValueError(f"{name} must be a sequence, not a string")
            values = tuple(getattr(self, name))
            if any(not isinstance(v, str) or not v or v != v.strip() for v in values):
                raise ValueError(f"{name} requires nonempty gene symbols or patterns")
            object.__setattr__(self, name, tuple(sorted(set(values))))
        if any(any(c in v for c in "*?[]") for v in self.gene_exceptions):
            raise ValueError("Gene exceptions must be exact symbols")

    def excludes(self, gene_name):
        return (gene_name not in self.gene_exceptions
                and any(fnmatchcase(gene_name, p) for p in self.gene_patterns))


@dataclass(frozen=True)
class CTAExpressionDecision(DataclassSerializable):
    """One input feature, including rejected/unknown evidence and its reason."""

    input_identifier: str
    value: float | None
    gene_id: str = ""
    gene_name: str = ""
    transcript_id: str = ""
    status: str = ""
    assessment: CTAAdmissionAssessment | None = None
    reference_resolution: CTAReferenceResolution | None = None

    def __post_init__(self):
        if self.value is not None and (not math.isfinite(self.value) or self.value < 0):
            raise ValueError("Expression decisions require finite non-negative values or None")
        if self.status == "admitted" and (
            self.assessment is None or not self.assessment.antigen.tumor_specificity.admits_construct
        ):
            raise ValueError("Admitted expression decisions require an admitted CTA assessment")


@dataclass(frozen=True)
class CTAExpressionResult(DataclassSerializable):
    """Persisted admission decisions with no expression/reference replay reads."""

    input_contract: CTAExpressionInput
    admission_policy: CTAAdmissionPolicy
    exclusion_policy: CTATargetExclusionPolicy
    input_sha256: str
    reference_assembly: str
    annotation_name: str
    annotation_version: str
    decisions: tuple[CTAExpressionDecision, ...]

    @property
    def admitted_antigens(self):
        return tuple(d.assessment.antigen for d in self.decisions
                     if d.status == "admitted")

    def save(self, path):
        Path(path).write_text(to_native_json(self) + "\n")

    @classmethod
    def load(cls, path):
        return from_native_json(Path(path).read_text(), cls)


def _feature_id(identifier, level):
    prefix = "ENSG" if level == "gene" else "ENST"
    if not isinstance(identifier, str) or not re.fullmatch(
            prefix + r"\d{11}(?:\.\d+)?", identifier):
        raise ValueError(f"Expected a human Ensembl {level} ID, got {identifier!r}")
    return identifier.split(".")[0]


def _expression_rows(path, contract):
    """Preflight omissions/invalid values before using Topiary's public loader."""
    if Path(path).suffix.lower() not in {".csv", ".tsv"}:
        raise ValueError("CTA expression input must be a CSV or TSV table")
    separator = "\t" if Path(path).suffix.lower() == ".tsv" else ","
    with Path(path).open(newline='') as stream:
        header = next(csv.reader((line for line in stream if line.strip() and not line.lstrip().startswith('#')),
                                 delimiter=separator), [])
    if len(set(header)) != len(header):
        raise ValueError('Expression input requires unique column names')
    raw = pd.read_csv(path, sep=separator, comment="#", dtype=str)
    missing = {contract.id_column, contract.value_column} - set(raw.columns)
    if missing:
        raise ValueError(f"Expression input is missing columns: {sorted(missing)}")
    # Topiary drops missing identifiers and coerces invalid numeric cells to
    # missing. Reject those cases here so no input feature disappears or an
    # invalid measurement becomes an unknown measurement silently.
    identifiers = raw[contract.id_column]
    normalized = [_feature_id(v, contract.measurement_level) for v in identifiers]
    if len(set(normalized)) != len(normalized):
        raise ValueError("Duplicate normalized expression feature IDs require explicit aggregation")
    for value in raw[contract.value_column].dropna():
        try:
            number = float(value)
        except ValueError as error:
            raise ValueError(f"Invalid expression measurement {value!r}") from error
        if not math.isfinite(number) or number < 0:
            raise ValueError("Expression measurements must be finite and non-negative")
    loaded = load_expression(path, id_col=contract.id_column,
                             val_cols=contract.value_column)
    if len(loaded) != len(raw):
        raise ValueError("Expression loader omitted input rows")
    return tuple(zip(identifiers, normalized, loaded[contract.value_column]))


def admit_cta_expression(path, *, input_contract, admission_policy, genome,
                         exclusion_policy=None):
    """Intersect patient gene/transcript expression with direct OncoRef CTAs.

    Gene inputs choose the OncoRef canonical transcript as a reference sequence;
    transcript inputs choose only the measured transcript. Measurements are
    never pooled across features, including identical protein sequences.
    Unresolvable transcripts are retained as held-out decisions. Reference
    failures (e.g. missing genome caches) propagate rather than becoming empty
    CTA categories. Only human Ensembl identifiers are supported here.
    """
    from oncoref.cta import cta_unfiltered_gene_ids

    if input_contract.expression_unit != admission_policy.expression_unit:
        raise ValueError("Expression input units do not match the CTA admission policy")
    if genome.species.latin_name != "homo_sapiens":
        raise ValueError("OncoRef CTA expression admission requires a human genome")
    reference = (str(genome.reference_name), str(genome.annotation_name),
                 str(genome.annotation_version))
    if any(v in {"", "None"} for v in reference):
        raise ValueError("CTA expression admission requires a versioned genome reference")
    exclusion_policy = exclusion_policy or CTATargetExclusionPolicy()
    fingerprint = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    rows = _expression_rows(path, input_contract)
    candidate_ids = {value.split(".")[0] for value in cta_unfiltered_gene_ids()}
    if not candidate_ids:
        raise ValueError("OncoRef returned an empty CTA candidate universe")
    # Establish annotation availability before resolving individual features.
    # A missing cache or malformed annotation must fail, not masquerade as
    # unknown transcripts in an otherwise successful, empty admission run.
    transcript_ids = frozenset(genome.transcript_ids())
    decisions = []
    for identifier, feature_id, abundance in rows:
        value = None if pd.isna(abundance) else float(abundance)
        decision = CTAExpressionDecision(identifier, value)
        transcript = None
        gene_id = feature_id
        if input_contract.measurement_level == "transcript":
            if feature_id not in transcript_ids:
                decisions.append(replace(decision, transcript_id=feature_id,
                                         status="unresolved_transcript"))
                continue
            transcript = genome.transcript_by_id(feature_id)
            gene_id = transcript.gene_id
        decision = replace(decision, gene_id=gene_id,
                           transcript_id=feature_id if transcript else "")
        if gene_id not in candidate_ids:
            decisions.append(replace(decision, status="not_cta"))
            continue
        cta = resolve_cta_reference_evidence(gene_id)
        decision = replace(decision, gene_name=cta.evidence.gene_name,
                           reference_resolution=cta)
        if exclusion_policy.excludes(decision.gene_name):
            decisions.append(replace(decision, status="excluded_target"))
            continue
        if value is None:
            decisions.append(replace(decision, status="unknown_expression"))
            continue
        transcript_id = (feature_id if transcript else cta.evidence.canonical_transcript_id)
        if not transcript_id:
            decisions.append(replace(decision, status="missing_canonical_transcript"))
            continue
        decision = replace(decision, transcript_id=transcript_id)
        if transcript is None:
            if transcript_id.split(".")[0] not in transcript_ids:
                decisions.append(replace(decision, status="unresolved_transcript"))
                continue
            transcript = genome.transcript_by_id(transcript_id.split(".")[0])
        if transcript.gene_id != gene_id:
            raise ValueError(f"Selected transcript {transcript_id} belongs to a different gene")
        versioned_id = identifier if input_contract.measurement_level == "transcript" else transcript_id
        if "." in versioned_id:
            version = versioned_id.split(".")[1]
            if str(transcript.transcript_version) != version:
                decisions.append(replace(decision, status="transcript_version_mismatch"))
                continue
        sequence = transcript.protein_sequence
        if not sequence:
            decisions.append(replace(decision, status="no_protein_sequence"))
            continue
        measurement = PatientTumorExpressionEvidence(
            gene_id=gene_id, sample_id=input_contract.sample_id, value=value,
            unit=input_contract.expression_unit,
            evidence_source=input_contract.evidence_source,
            evidence_version=input_contract.evidence_version, assay=input_contract.assay,
            measurement_level=input_contract.measurement_level,
            transcript_id=identifier if input_contract.measurement_level == "transcript" else "",
            input_identifier=identifier, input_sha256=fingerprint)
        assessment = assess_cta_antigen(
            amino_acids=sequence, gene_id=gene_id, tumor_expression=measurement,
            policy=admission_policy, transcript_ids=(transcript_id,),
            protein_ids=(transcript.protein_id,) if transcript.protein_id else (),
            source_identifier=f"{input_contract.sample_id}:{identifier}")
        metadata = (
            ("expression_id_namespace", input_contract.id_namespace),
            ("reference_assembly", reference[0]), ("annotation_name", reference[1]),
            ("annotation_version", reference[2]),
            ("sequence_selection", "oncoref_canonical" if input_contract.measurement_level == "gene"
             else "measured_transcript"),
        )
        assessment = replace(assessment, antigen=replace(
            assessment.antigen, source_metadata=metadata))
        status = "admitted" if assessment.antigen.tumor_specificity.admits_construct else "held_out"
        decisions.append(replace(decision, status=status, assessment=assessment))
    if fingerprint != hashlib.sha256(Path(path).read_bytes()).hexdigest():
        raise ValueError("Expression input changed during admission")
    return CTAExpressionResult(input_contract, admission_policy, exclusion_policy,
                               fingerprint, *reference, tuple(decisions))
