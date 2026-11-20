# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Matched-normal DNA and normal/blood RNA observations alongside attestations.

Preserves per-library, per-pipeline matched-normal observations separately by
assay (``DNA`` vs. ``RNA``), treats repeated pipelines on one library as
dependent analyses (never summing their counts), and surfaces conflicts between
asserted somatic status and observed normal alternate support without labelling
alternate reads alone as germline or contamination.
"""

from collections import defaultdict
from dataclasses import dataclass, replace
from typing import Optional

from serializable import DataclassSerializable


NORMAL_ASSAYS = frozenset({"DNA", "RNA"})
NORMAL_CONFLICT_POLICIES = frozenset({"hold_for_review", "flag_for_review"})

STATE_MISSING = "missing"
STATE_ZERO_ALT_ADEQUATE = "zero_alt_adequate_coverage"
STATE_LOW_COVERAGE = "low_coverage"
STATE_OBSERVED_ALT = "observed_alt_support"

ADEQUACY_MISSING = "missing"
ADEQUACY_LOW = "low_coverage"
ADEQUACY_ADEQUATE = "adequate_coverage"

NORMAL_ALT_INTERPRETATION_NOTE = (
    "Observed alternate support in matched normal or blood conflicts with an "
    "unqualified somatic assumption and requires review; alternate reads alone "
    "do not establish germline origin, clonal hematopoiesis, tumor-in-normal "
    "carryover, or sequencing contamination."
)


@dataclass(frozen=True)
class MatchedNormalObservation(DataclassSerializable):
    """One pipeline observation of a locus/feature in a normal or blood library."""

    assay: str
    library_id: str
    pipeline: str
    pipeline_version: str = ""
    source_identity: str = ""
    tissue: str = "normal"
    depth: Optional[int] = None
    alt_count: Optional[int] = None
    min_adequate_depth: int = 20

    def __post_init__(self):
        assay = str(self.assay or "").strip().upper()
        if assay not in NORMAL_ASSAYS:
            raise ValueError(
                f"Matched-normal assay must be one of {sorted(NORMAL_ASSAYS)}, got {self.assay!r}"
            )
        object.__setattr__(self, "assay", assay)
        if not isinstance(self.library_id, str) or not self.library_id.strip():
            raise ValueError("Matched-normal observation requires a nonempty library_id")
        if not isinstance(self.pipeline, str) or not self.pipeline.strip():
            raise ValueError("Matched-normal observation requires a nonempty pipeline")
        if not isinstance(self.tissue, str) or not self.tissue.strip():
            raise ValueError("Matched-normal observation requires a nonempty tissue label")
        if type(self.min_adequate_depth) is not int or self.min_adequate_depth < 1:
            raise ValueError("min_adequate_depth must be a positive integer")
        for label in ("depth", "alt_count"):
            value = getattr(self, label)
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError(f"Matched-normal {label} must be a nonnegative integer or None")
        if (
            self.depth is not None
            and self.alt_count is not None
            and self.alt_count > self.depth
        ):
            raise ValueError("Matched-normal alt_count cannot exceed depth")

    @property
    def coverage_adequacy(self) -> str:
        if self.depth is None:
            return ADEQUACY_MISSING
        if self.depth < self.min_adequate_depth:
            return ADEQUACY_LOW
        return ADEQUACY_ADEQUATE

    @property
    def observation_state(self) -> str:
        if self.depth is None or self.alt_count is None:
            return STATE_MISSING
        if self.alt_count > 0:
            return STATE_OBSERVED_ALT
        if self.depth < self.min_adequate_depth:
            return STATE_LOW_COVERAGE
        return STATE_ZERO_ALT_ADEQUATE


def validate_normal_observations(observations) -> tuple[MatchedNormalObservation, ...]:
    records = tuple(observations or ())
    seen = {}
    for record in records:
        if not isinstance(record, MatchedNormalObservation):
            raise TypeError("normal_observations requires MatchedNormalObservation records")
        key = (
            record.assay,
            record.library_id,
            record.pipeline,
            record.pipeline_version,
            record.source_identity,
            record.tissue,
        )
        prior = seen.get(key)
        if prior is not None and prior != record:
            raise ValueError(
                f"Conflicting matched-normal observation for {key!r}"
            )
        seen[key] = record
    return tuple(seen.values())


def group_normal_observations_by_library(observations) -> list[dict]:
    """Group dependent pipeline runs by ``(assay, library_id)`` without summing."""
    grouped = defaultdict(list)
    for obs in validate_normal_observations(observations):
        grouped[(obs.assay, obs.library_id)].append(obs)
    result = []
    for (assay, library_id) in sorted(grouped):
        runs = grouped[(assay, library_id)]
        states = [r.observation_state for r in runs]
        if STATE_OBSERVED_ALT in states:
            summary_state = STATE_OBSERVED_ALT
        elif STATE_ZERO_ALT_ADEQUATE in states:
            summary_state = STATE_ZERO_ALT_ADEQUATE
        elif STATE_LOW_COVERAGE in states:
            summary_state = STATE_LOW_COVERAGE
        else:
            summary_state = STATE_MISSING
        result.append(dict(
            assay=assay,
            library_id=library_id,
            tissues=list(dict.fromkeys(r.tissue for r in runs)),
            dependent_pipeline_count=len(runs),
            summary_state=summary_state,
            low_coverage_present=any(r.coverage_adequacy == ADEQUACY_LOW for r in runs),
            max_observed_alt=max((r.alt_count for r in runs if r.alt_count is not None), default=None),
            pipelines=[
                dict(
                    pipeline=r.pipeline,
                    pipeline_version=r.pipeline_version,
                    source_identity=r.source_identity,
                    tissue=r.tissue,
                    depth=r.depth,
                    alt_count=r.alt_count,
                    min_adequate_depth=r.min_adequate_depth,
                    coverage_adequacy=r.coverage_adequacy,
                    observation_state=r.observation_state,
                )
                for r in runs
            ],
        ))
    return result


def summarize_normal_assessment(attestation, observations) -> dict:
    """Summarize normal DNA/RNA states and somatic-attestation conflicts."""
    libraries = group_normal_observations_by_library(observations)
    dna_libs = [lib for lib in libraries if lib["assay"] == "DNA"]
    rna_libs = [lib for lib in libraries if lib["assay"] == "RNA"]

    def _assay_state(libs):
        if not libs:
            return STATE_MISSING
        states = {lib["summary_state"] for lib in libs}
        if STATE_OBSERVED_ALT in states:
            return STATE_OBSERVED_ALT
        if STATE_ZERO_ALT_ADEQUATE in states:
            return STATE_ZERO_ALT_ADEQUATE
        if STATE_LOW_COVERAGE in states:
            return STATE_LOW_COVERAGE
        return STATE_MISSING

    dna_state = _assay_state(dna_libs)
    rna_state = _assay_state(rna_libs)
    asserts_somatic = (
        "somatic" in (attestation.evidence_kind or "").lower()
        or "somatic" in (attestation.rationale_code or "").lower()
    )
    conflict_reasons = []
    if asserts_somatic and dna_state == STATE_OBSERVED_ALT:
        conflict_reasons.append("observed_normal_dna_alt_support")
    if asserts_somatic and rna_state == STATE_OBSERVED_ALT:
        conflict_reasons.append("observed_normal_rna_alt_support")

    review_reasons = list(conflict_reasons)
    if asserts_somatic and dna_state == STATE_MISSING:
        review_reasons.append("matched_normal_dna_missing")
    elif asserts_somatic and dna_state == STATE_LOW_COVERAGE:
        review_reasons.append("matched_normal_dna_low_coverage")
    if any(lib["low_coverage_present"] and lib["summary_state"] == STATE_OBSERVED_ALT for lib in libraries):
        review_reasons.append("observed_normal_alt_at_low_coverage")

    return dict(
        asserts_somatic=asserts_somatic,
        normal_dna_state=dna_state,
        normal_rna_state=rna_state,
        conflict=bool(conflict_reasons),
        conflict_reasons=conflict_reasons,
        review_reasons=review_reasons,
        libraries=libraries,
        note=NORMAL_ALT_INTERPRETATION_NOTE if conflict_reasons else "",
    )


def adjudicate_somatic_attestation_with_normal(
    attestation,
    observations,
    *,
    conflict_policy: str = "hold_for_review",
):
    """Apply explicit conflict/review policy to ``attestation`` without erasing history."""
    if conflict_policy not in NORMAL_CONFLICT_POLICIES:
        raise ValueError(
            f"Unknown normal conflict_policy {conflict_policy!r}; "
            f"expected one of {sorted(NORMAL_CONFLICT_POLICIES)}"
        )
    records = validate_normal_observations(observations)
    if not records:
        return attestation
    summary = summarize_normal_assessment(attestation, records)
    if not summary["review_reasons"]:
        return attestation

    base_rationale = attestation.rationale_code or attestation.evidence_kind
    added = ",".join(summary["review_reasons"])
    rationale = f"{base_rationale};normal_review:{added}" if base_rationale else f"normal_review:{added}"

    if summary["conflict"] and conflict_policy == "hold_for_review" and attestation.status == "admitted":
        return replace(
            attestation,
            status="held_out",
            rationale_code=rationale,
            requires_review=True,
            override_reason="",
        )
    return replace(
        attestation,
        rationale_code=rationale if summary["conflict"] or not attestation.rationale_code else rationale,
        requires_review=True,
    )


_STATE_LABELS = {
    STATE_ZERO_ALT_ADEQUATE: "0 alt",
    STATE_OBSERVED_ALT: "alt reads observed",
    STATE_LOW_COVERAGE: "low coverage",
    STATE_MISSING: "unmeasured",
}


def _dedup_pipeline_runs(pipelines):
    seen = set()
    unique = []
    for run in pipelines:
        key = (
            run["pipeline"],
            run.get("pipeline_version", ""),
            run["alt_count"],
            run["depth"],
            run["coverage_adequacy"],
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(run)
    return unique


def _format_run_compact(run):
    if run["depth"] is None or run["alt_count"] is None:
        return f"{run['pipeline']} unmeasured"
    suffix = " (low cov)" if run["coverage_adequacy"] == ADEQUACY_LOW else ""
    return f"{run['pipeline']} {run['alt_count']}/{run['depth']}{suffix}"


def _format_assay_compact(assay_name, libs, *, inline_libraries=False):
    if not libs:
        return f"{assay_name} unmeasured"
    if inline_libraries:
        parts = []
        for lib in libs:
            runs = ", ".join(_format_run_compact(r) for r in lib["dedup_pipelines"])
            parts.append(f"{lib['library_id']}: {runs}")
        states = {lib["summary_state"] for lib in libs}
        state_label = (
            _STATE_LABELS[STATE_OBSERVED_ALT] if STATE_OBSERVED_ALT in states
            else _STATE_LABELS[STATE_ZERO_ALT_ADEQUATE] if STATE_ZERO_ALT_ADEQUATE in states
            else _STATE_LABELS[STATE_LOW_COVERAGE] if STATE_LOW_COVERAGE in states
            else _STATE_LABELS[STATE_MISSING]
        )
        return f"{assay_name} {state_label} ({'; '.join(parts)})"

    measured_depths = [
        r["depth"] for lib in libs for r in lib["dedup_pipelines"] if r["depth"] is not None
    ]
    alt_libs = [lib for lib in libs if lib["summary_state"] == STATE_OBSERVED_ALT]
    if alt_libs:
        max_alt = max(
            (lib["max_observed_alt"] for lib in alt_libs if lib["max_observed_alt"] is not None),
            default=0,
        )
        # Show the worst (highest alt / lowest depth at max_alt) ratio for context
        worst_depth = min(
            (
                r["depth"] for lib in alt_libs for r in lib["dedup_pipelines"]
                if r["alt_count"] == max_alt and r["depth"] is not None
            ),
            default=None,
        )
        ratio = f"{max_alt}/{worst_depth}" if worst_depth is not None else str(max_alt)
        return f"{assay_name} alt reads in {len(alt_libs)}/{len(libs)} libraries (max {ratio})"

    adequate_zero = [lib for lib in libs if lib["summary_state"] == STATE_ZERO_ALT_ADEQUATE]
    if adequate_zero:
        if measured_depths:
            lo, hi = min(measured_depths), max(measured_depths)
            depth_span = f"{lo}x" if lo == hi else f"{lo}-{hi}x"
            return f"{assay_name} 0 alt across {len(adequate_zero)}/{len(libs)} libraries ({depth_span})"
        return f"{assay_name} 0 alt across {len(adequate_zero)}/{len(libs)} libraries"

    low_cov = [lib for lib in libs if lib["summary_state"] == STATE_LOW_COVERAGE]
    if low_cov:
        max_d = max(measured_depths) if measured_depths else 0
        return f"{assay_name} low coverage ({len(low_cov)}/{len(libs)} libraries, max {max_d}x)"
    return f"{assay_name} unmeasured ({len(libs)} libraries)"


def format_normal_product_summary(summary: dict, *, inline_limit: int = 4) -> dict:
    """Format a matched-normal assessment into concise product-view prose + optional library rows."""
    libraries = []
    for lib in summary["libraries"]:
        dedup = _dedup_pipeline_runs(lib["pipelines"])
        libraries.append(dict(
            assay=lib["assay"],
            library_id=lib["library_id"],
            tissue=", ".join(lib["tissues"]),
            state_label=_STATE_LABELS.get(lib["summary_state"], lib["summary_state"]),
            summary_state=lib["summary_state"],
            max_observed_alt=lib["max_observed_alt"],
            dedup_pipelines=dedup,
            runs_text=", ".join(_format_run_compact(r) for r in dedup),
        ))
    dna_libs = [lib for lib in libraries if lib["assay"] == "DNA"]
    rna_libs = [lib for lib in libraries if lib["assay"] == "RNA"]
    inline = len(libraries) <= inline_limit
    dna_text = _format_assay_compact("DNA", dna_libs, inline_libraries=inline)
    rna_text = _format_assay_compact("RNA", rna_libs, inline_libraries=inline)
    headline = f"Matched normal: {dna_text}; {rna_text}."
    if summary["conflict"]:
        headline += " Review required."
    text_line = f"{headline} {summary['note']}" if summary["note"] else headline
    return dict(
        headline=headline,
        note=summary["note"],
        text_line=text_line,
        collapsible=not inline,
        library_count=len(libraries),
        libraries=libraries,
    )
