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

"""Within-library Isovar read-set unions for shared candidate vaccine sequences.

When multiple source genes or transcripts encode the exact same candidate
vaccine amino-acid sequence, adding per-gene counts can double-count shared
reads, and taking one source's count omits distinct reads supporting other
origins. This module combines scoped Isovar ``RNAFeatureEvidence`` records
within each ``evidence_scope = (sample_id, source_id)`` via
``isovar.union_rna_support``, preserving:

- separate libraries and timepoints (never pooled across ``evidence_scope``);
- individually spanning reads vs. fragments supported jointly by paired mates;
- verified common whole-library denominators for ``per_million`` quantities;
- per-source overlap and exclusive counts alongside supporting set identities;
- explicit unavailability of effective-length-normalized group quantities,
  because a raw read-set union across alternative transcript paths does not
  define a shared effective sampling length.
"""

from collections import defaultdict
from dataclasses import replace
import json

from isovar import union_rna_support
from isovar.rna_evidence import content_identifier

from .config.defaults import (
    SHARED_SEQUENCE_RNA_POLICIES,
    SHARED_SEQUENCE_RNA_POLICY_SOURCE,
    SHARED_SEQUENCE_RNA_POLICY_UNION as SHARED_SEQUENCE_RNA_POLICY_UNION,
)
from .rna_support import applicable_rna_support


SHARED_RNA_UNION_SCHEMA = "vaxrank.shared_sequence_rna_union.v1"

SHARED_RNA_AVAILABLE_UNITS = frozenset({
    "rna_shared_reads",
    "rna_shared_fragments",
    "rna_shared_library_reads",
    "rna_shared_library_fragments",
    "rna_shared_reads_per_million",
    "rna_shared_fragments_per_million",
})
SHARED_RNA_UNDEFINED_LENGTH_BINDINGS = frozenset({
    "rna_shared_effective_length_reads",
    "rna_shared_effective_length_fragments",
    "rna_shared_reads_per_effective_kb",
    "rna_shared_fragments_per_effective_kb",
    "rna_shared_reads_per_effective_kb_per_million",
    "rna_shared_fragments_per_effective_kb_per_million",
})
SHARED_RNA_BINDINGS = SHARED_RNA_AVAILABLE_UNITS | SHARED_RNA_UNDEFINED_LENGTH_BINDINGS


def validate_shared_sequence_rna_policy(policy: str) -> str:
    if not isinstance(policy, str) or policy not in SHARED_SEQUENCE_RNA_POLICIES:
        raise ValueError(
            f"Unknown shared_sequence_rna_policy {policy!r}; "
            f"expected one of {sorted(SHARED_SEQUENCE_RNA_POLICIES)}"
        )
    return policy


def compute_shared_sequence_rna_unions(antigens) -> tuple[dict, ...]:
    """Compute within-library RNA set unions across origins of one sequence.

    Parameters
    ----------
    antigens : iterable of VaccineAntigen
        Candidate-window ``VaccineAntigen`` objects sharing one exact
        ``amino_acids`` sequence.

    Returns
    -------
    tuple of dict
        One validated ``vaxrank.shared_sequence_rna_union.v1`` record per
        ``evidence_scope`` present in applicable ``RNAFeatureEvidence`` records,
        ordered deterministically by ``evidence_scope``.
    """
    antigens = list(antigens)
    if not antigens:
        return ()
    sequences = {antigen.amino_acids for antigen in antigens}
    if len(sequences) != 1:
        raise ValueError(
            "Shared-sequence RNA union requires candidate antigens with the exact same sequence"
        )
    sequence = next(iter(sequences))

    by_scope = defaultdict(list)
    seen_member_keys = defaultdict(set)
    for antigen in antigens:
        for evidence in applicable_rna_support(antigen):
            measurement = evidence.measurement
            scope = tuple(measurement["evidence_scope"])
            dedup_key = (
                antigen.source_identifier,
                antigen.gene_id,
                antigen.gene_name,
                antigen.transcript_ids,
                antigen.protein_ids,
                evidence.scope,
                evidence.protein_interval,
                measurement["measurement_id"],
            )
            if dedup_key in seen_member_keys[scope]:
                continue
            seen_member_keys[scope].add(dedup_key)
            by_scope[scope].append((antigen, evidence))

    unions = []
    for scope in sorted(by_scope):
        members = by_scope[scope]
        fragment_modes = {
            ev.measurement["feature"].get("fragment_mode", "union")
            for _, ev in members
        }
        if len(fragment_modes) > 1:
            raise ValueError(
                f"Cannot union shared-sequence RNA features with conflicting fragment_mode "
                f"{sorted(fragment_modes)} in evidence_scope {list(scope)!r}"
            )
        fragment_mode = next(iter(fragment_modes))

        for unit in ("reads", "fragments"):
            denominators = {
                ev.measurement["values"]["library_" + unit]
                for _, ev in members
            }
            if len(denominators) > 1:
                raise ValueError(
                    f"Conflicting library_{unit} denominators {sorted(denominators, key=str)} "
                    f"within evidence_scope {list(scope)!r}"
                )

        measurements = [ev.measurement for _, ev in members]
        united = union_rna_support(measurements)

        # A fragment is individually read-supported in the union iff at least
        # one contributing feature observed an individually spanning read for
        # that fragment. Joint-mate-only fragments remain distinct.
        single_read_fragments = set().union(*(
            set(m["fragment_ids"]) - set(m["joint_only_fragment_ids"])
            for m in measurements
        ))
        joint_only_fragment_ids = sorted(set(united["fragment_ids"]) - single_read_fragments)
        single_count = united["fragments"] - len(joint_only_fragment_ids)
        if not single_count <= united["reads"] <= 2 * single_count:
            raise ValueError(
                "Shared-sequence RNA union read and fragment cardinalities disagree"
            )

        member_summaries = []
        any_shared_units = False
        for index, (antigen, ev) in enumerate(members):
            m = ev.measurement
            own_reads = set(m["read_ids"])
            own_frags = set(m["fragment_ids"])
            other_reads = set().union(*(
                set(m2["read_ids"]) for j, m2 in enumerate(measurements) if j != index
            )) if len(measurements) > 1 else set()
            other_frags = set().union(*(
                set(m2["fragment_ids"]) for j, m2 in enumerate(measurements) if j != index
            )) if len(measurements) > 1 else set()
            shared_reads = len(own_reads & other_reads)
            shared_frags = len(own_frags & other_frags)
            if shared_reads or shared_frags:
                any_shared_units = True
            member_summaries.append(dict(
                source_identifier=antigen.source_identifier,
                gene_name=antigen.gene_name,
                gene_id=antigen.gene_id,
                transcript_ids=list(antigen.transcript_ids),
                protein_ids=list(antigen.protein_ids),
                applicability=ev.scope,
                protein_interval=(
                    list(ev.protein_interval) if ev.protein_interval is not None else None
                ),
                feature=m["feature"],
                measurement_id=m["measurement_id"],
                reads=len(own_reads),
                fragments=len(own_frags),
                joint_only_fragments=len(m["joint_only_fragment_ids"]),
                reads_shared_with_other_sources=shared_reads,
                reads_exclusive_to_source=len(own_reads - other_reads),
                fragments_shared_with_other_sources=shared_frags,
                fragments_exclusive_to_source=len(own_frags - other_frags),
                effective_length_reads=m["values"].get("effective_length_reads"),
                effective_length_fragments=m["values"].get("effective_length_fragments"),
            ))

        values = {}
        unavailable = {}
        for unit in ("reads", "fragments"):
            denominator = measurements[0]["values"]["library_" + unit]
            count = united[unit]
            sum_of_sources = sum(m["values"][unit] for m in measurements)
            per_million = (
                None if denominator is None or denominator == 0
                else 1e6 * count / denominator
            )
            if per_million is None:
                unavailable[unit + "_per_million"] = (
                    "whole_library_denominator_unavailable"
                    if denominator is None else "feature_unobservable"
                )
            values[unit] = count
            values["sum_of_source_" + unit] = sum_of_sources
            values["library_" + unit] = denominator
            values[unit + "_per_million"] = per_million
            for suffix in (
                "effective_length_" + unit,
                unit + "_per_effective_kb",
                unit + "_per_effective_kb_per_million",
            ):
                values[suffix] = None
                unavailable[suffix] = "shared_sequence_effective_length_undefined"

        distinct_origins = {
            (ag.source_identifier, ag.gene_id, ag.gene_name, ag.transcript_ids)
            for ag, _ in members
        }
        limitations = []
        if len(distinct_origins) > 1:
            limitations.append("shared_sequence_ambiguous_origin_attribution")
        if any_shared_units:
            limitations.append("multi_source_shared_read_or_fragment_identities")
        if any(ev.scope == "event" for _, ev in members):
            limitations.append("event_scope_does_not_assert_full_window_coverage")
        limitations.append("shared_sequence_effective_length_undefined")

        record = dict(
            schema=SHARED_RNA_UNION_SCHEMA,
            sequence=sequence,
            evidence_scope=list(scope),
            fragment_mode=fragment_mode,
            values=values,
            unavailable=unavailable,
            read_ids=united["read_ids"],
            fragment_ids=united["fragment_ids"],
            joint_only_fragment_ids=joint_only_fragment_ids,
            evidence_set_id=united["evidence_set_id"],
            fragment_evidence_set_id=content_identifier(
                "fragment_evidence", united["fragment_ids"]
            ),
            members=member_summaries,
            attribution_limitations=limitations,
        )
        record["union_id"] = content_identifier("shared_rna_union", record)
        unions.append(json.loads(json.dumps(record, sort_keys=True, allow_nan=False)))
    return tuple(unions)


def attach_shared_sequence_rna_support(peptides, *, policy: str = SHARED_SEQUENCE_RNA_POLICY_SOURCE):
    """Attach within-library shared-sequence RNA unions across ``VaccinePeptide``s.

    When ``policy == 'source'``, returns ``peptides`` unchanged so default
    per-source scoring and historical replay remain untouched.
    """
    validate_shared_sequence_rna_policy(policy)
    peptides = list(peptides)
    if policy == SHARED_SEQUENCE_RNA_POLICY_SOURCE or not peptides:
        return peptides
    by_sequence = defaultdict(list)
    for vp in peptides:
        by_sequence[vp.amino_acids].append(vp.antigen)
    unions_by_sequence = {
        seq: compute_shared_sequence_rna_unions(antigens)
        for seq, antigens in by_sequence.items()
    }
    return [
        replace(vp, shared_sequence_rna_support=unions_by_sequence[vp.amino_acids])
        if unions_by_sequence[vp.amino_acids] != getattr(vp, "shared_sequence_rna_support", ())
        else vp
        for vp in peptides
    ]


def attach_shared_sequence_rna_to_groups(grouped_peptides, *, policy: str = SHARED_SEQUENCE_RNA_POLICY_SOURCE):
    """Attach within-library shared-sequence RNA unions across ``[(key, list[VaccinePeptide])]``."""
    validate_shared_sequence_rna_policy(policy)
    grouped = [(key, list(vps)) for key, vps in grouped_peptides]
    if policy == SHARED_SEQUENCE_RNA_POLICY_SOURCE or not grouped:
        return grouped
    flat = [vp for _, vps in grouped for vp in vps]
    updated = attach_shared_sequence_rna_support(flat, policy=policy)
    offset = 0
    result = []
    for key, vps in grouped:
        count = len(vps)
        result.append((key, updated[offset:offset + count]))
        offset += count
    return result


def shared_sequence_rna_bindings(vaccine_peptide) -> dict:
    """Return ``rna_shared_*`` scalar bindings for the combined-score DSL."""
    records = getattr(vaccine_peptide, "shared_sequence_rna_support", ()) or ()
    values = records[0]["values"] if len(records) == 1 else {}
    return {
        name: values.get(name.removeprefix("rna_shared_"))
        for name in SHARED_RNA_BINDINGS
    }
