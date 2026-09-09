# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Non-mutating assessment and lossless reporting of explicit constructs."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import jinja2
from serializable import DataclassSerializable

from .construct_sequence import ConstructSequence
from .native_serialization import from_native_json, to_native_json
from .safety_assessment import (
    ConstructBoundary,
    SafetyAssessmentError,
    WindowSafetyAssessment,
    assess_vaccine_antigen_window,
)


@dataclass(frozen=True)
class ConstructAudit(DataclassSerializable):
    """Observed prediction inventory, never a clinical safety verdict.

    ``predictions_returned`` does not promise complete model/HLA/length coverage.
    Processing compartments remain separately unassessed until site-resolved
    evidence is provided; MHC binding and half-life numbers cannot fill that gap.
    """

    construct: ConstructSequence
    mhc_status: str
    mhc_assessment: Optional[WindowSafetyAssessment] = None
    reason_codes: tuple[str, ...] = field(default_factory=tuple)
    error_message: str = ""
    unassessed_processing_compartments: tuple[str, ...] = (
        "proteasomal", "extracellular_serum", "endolysosomal")

    def __post_init__(self):
        object.__setattr__(self, "reason_codes", tuple(self.reason_codes))
        object.__setattr__(self, "unassessed_processing_compartments",
                           tuple(self.unassessed_processing_compartments))
        if not isinstance(self.construct, ConstructSequence):
            raise ValueError("Construct audit requires a ConstructSequence")
        if self.mhc_status not in {"unassessed", "failed", "no_predictions", "predictions_returned"}:
            raise ValueError("Unknown construct MHC assessment status")
        if self.mhc_status in {"unassessed", "failed"}:
            if self.mhc_assessment is not None or not self.reason_codes:
                raise ValueError("Unassessed/failed audit needs reasons and no prediction inventory")
        else:
            if not isinstance(self.mhc_assessment, WindowSafetyAssessment):
                raise ValueError("Construct audit requires its prediction inventory")
            if (self.mhc_assessment.antigen != self.construct.assessment_antigen()
                    or self.mhc_assessment.window_start_offset != 0
                    or self.mhc_assessment.window_end_offset != len(self.construct.sequence)):
                raise ValueError("Audit predictions do not describe the complete final construct")
            if bool(self.mhc_assessment.ligands) != (self.mhc_status == "predictions_returned"):
                raise ValueError("Audit MHC status disagrees with its prediction inventory")
        if self.mhc_status == "failed" and not self.error_message:
            raise ValueError("A failed construct audit requires an error message")
        # This version has no site-evidence field. Do not allow a serialized
        # record to remove a missing compartment and thereby imply assessment.
        if set(self.unassessed_processing_compartments) != {
                "proteasomal", "extracellular_serum", "endolysosomal"}:
            raise ValueError("Processing remains unassessed without site-resolved evidence")


def audit_construct_sequence(construct, mhc_predictor=None, *, genome=None,
                            reference_proteome=None):
    """Assess the actual final context through Topiary and existing self policy.

    Native-context predictions are deliberately not accepted as an argument.
    Genomes retain full self-source provenance; membership-only reference inputs
    stay explicitly incomplete. No ranking, filtering or manufacturing occurs.
    """
    reason = construct.sequence_prediction_unassessed_reason
    if reason is None and mhc_predictor is None:
        reason = "mhc_predictor_not_requested"
    if reason is None and genome is None and reference_proteome is None:
        reason = "self_reference_unavailable"
    if reason:
        return ConstructAudit(construct, "unassessed", reason_codes=(reason,))
    antigen = construct.assessment_antigen()
    boundaries = tuple(ConstructBoundary(offset, "construct_edit", "construct_edit")
                       for offset in construct.modification_boundaries)
    try:
        assessment = assess_vaccine_antigen_window(
            mhc_predictor, antigen, genome=genome,
            reference_proteome=reference_proteome, construct_boundaries=boundaries)
    except SafetyAssessmentError as error:
        return ConstructAudit(
            construct, "failed", reason_codes=("mhc_prediction_failed",),
            error_message=str(error))
    return ConstructAudit(
        construct, "predictions_returned" if assessment.ligands else "no_predictions",
        mhc_assessment=assessment, reason_codes=assessment.reason_codes)


def save_construct_sequences(constructs, path):
    """Write a lossless allowlisted native file for explicit construct inputs."""
    records = list(constructs)
    if any(not isinstance(record, ConstructSequence) for record in records):
        raise ValueError("Construct input file requires ConstructSequence records")
    Path(path).write_text(to_native_json(records), encoding="utf8")


def load_construct_sequences(path):
    """Load attributed constructs without allowing input-selected imports."""
    try:
        records = from_native_json(Path(path).read_text(encoding="utf8"), list)
        for i, record in enumerate(records):
            if not isinstance(record, ConstructSequence):
                raise ValueError(f"Record {i + 1} is not a ConstructSequence")
        return records
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid construct input {str(path)!r}: {error}") from error


def write_construct_audits(audits, *, json_path, html_path=None):
    """Write the full native graph and an escaped, source-aware review report."""
    audits = list(audits)
    if any(not isinstance(audit, ConstructAudit) for audit in audits):
        raise ValueError("Construct report requires ConstructAudit records")
    payload = to_native_json(audits)
    rendered = None
    if html_path is not None:
        environment = jinja2.Environment(
            loader=jinja2.PackageLoader("vaxrank", "templates"),
            autoescape=jinja2.select_autoescape(("html", "xml")),
        )
        rendered = environment.get_template("construct_audit.html").render(audits=audits)
    Path(json_path).write_text(payload, encoding="utf8")
    if rendered is not None:
        Path(html_path).write_text(rendered, encoding="utf8")
