# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Explicit native-antigen versus manufactured-sequence provenance.

Here, native means the unmodified antigen (including any tumor variant), not
the wild-type protein. Mapping a residue to that antigen does not itself prove
RNA coverage. Added and substituted residues never acquire a native coordinate.
No alignment, solubility heuristic or automatic sequence design is performed.
"""

from dataclasses import dataclass, field, replace
import hashlib
import json
from typing import Optional

from serializable import DataclassSerializable

from .amino_acids import validate_amino_acid_sequence
from .vaccine_antigen import (
    AminoAcidInterval,
    TargetableMask,
    VaccineAntigen,
)


def _interval(start, end, label):
    if type(start) is not int or type(end) is not int:
        raise ValueError(f"{label} coordinates must be integers")
    if start < 0 or end < start:
        raise ValueError(f"{label} must be a nonnegative half-open interval")


def _evidence(evidence, rationale=None):
    if not isinstance(evidence, ConstructEvidence):
        raise ValueError("Construct evidence must be an attributed ConstructEvidence record")
    if rationale is not None and not isinstance(rationale, ConstructEvidence):
        raise ValueError("Modification rationale must be an attributed ConstructEvidence record")


@dataclass(frozen=True)
class ConstructEvidence(DataclassSerializable):
    """One attributed observation or explanation, not a safety attestation.

    ``documented`` identifies a published/recorded observation; ``reported``
    can attribute a user's explanation; ``inferred`` marks a hypothesis.
    Separate instances describe sequence evidence and modification rationale.
    """

    source: str
    evidence_level: str
    description: str
    provider: str = ""
    vaccine_version: str = ""

    def __post_init__(self):
        if any(not isinstance(value, str) for value in (
                self.source, self.evidence_level, self.description,
                self.provider, self.vaccine_version)):
            raise ValueError("Construct evidence fields must be strings")
        if not self.source.strip() or not self.description.strip():
            raise ValueError("Construct evidence requires a source and description")
        if self.evidence_level not in {"documented", "reported", "inferred"}:
            raise ValueError("Unknown construct evidence level")


@dataclass(frozen=True)
class ConstructSequenceEdit(DataclassSerializable):
    """Replace native antigen [start, end) by replacement amino acids.

    Zero width denotes insertion; an empty replacement denotes deletion.
    Coordinates are absolute within the recorded native antigen, not relative
    to either the selected window or a previously applied edit. Rationale is
    independently attributed and may be unavailable.
    """

    start: int
    end: int
    replacement: str
    evidence: ConstructEvidence
    rationale: Optional[ConstructEvidence] = None

    def __post_init__(self):
        _interval(self.start, self.end, "Sequence edit")
        _evidence(self.evidence, self.rationale)
        if not isinstance(self.replacement, str):
            raise ValueError("Sequence edit replacement must be a string")
        if self.replacement:
            validate_amino_acid_sequence(self.replacement, "Sequence edit")
        elif self.start == self.end:
            raise ValueError("A sequence edit must change at least one residue")

    @property
    def kind(self):
        if self.start == self.end:
            return "insertion"
        return "substitution" if self.replacement else "deletion"


@dataclass(frozen=True)
class ConstructChemicalModification(DataclassSerializable):
    """A named chemical change in final-sequence coordinates.

    [0, 0) and [length, length) denote the N and C termini. Nonempty intervals
    identify modified residues. This does not rewrite the canonical sequence
    or imply that a sequence-only predictor supports the chemical form.
    """

    start: int
    end: int
    name: str
    evidence: ConstructEvidence
    rationale: Optional[ConstructEvidence] = None

    def __post_init__(self):
        _interval(self.start, self.end, "Chemical modification")
        _evidence(self.evidence, self.rationale)
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("Chemical modification name is required")


@dataclass(frozen=True)
class ConstructSequence(DataclassSerializable):
    """One final antigen construct with an explicit native-window mapping.

    For mRNA/DNA, ``sequence`` is the encoded amino-acid sequence, not nucleotide
    sequence. Chemically synthesized modifications are only valid for peptide
    products. An unresolved mapping is reportable but cannot enter sequence
    prediction as though its target coordinates were established.
    """

    name: str
    sequence: str
    modality: str
    evidence: ConstructEvidence
    native_antigen: Optional[VaccineAntigen] = None
    native_start: Optional[int] = None
    native_end: Optional[int] = None
    sequence_edits: tuple[ConstructSequenceEdit, ...] = field(default_factory=tuple)
    chemical_modifications: tuple[ConstructChemicalModification, ...] = field(
        default_factory=tuple)
    mapping_status: str = "resolved"
    mapping_reason: str = ""

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("Construct name is required")
        if not isinstance(self.sequence, str):
            raise ValueError("Construct sequence must be a string")
        _evidence(self.evidence)
        if self.native_antigen is not None and not isinstance(self.native_antigen, VaccineAntigen):
            raise ValueError("Construct native antigen must be a VaccineAntigen")
        validate_amino_acid_sequence(self.sequence, "Construct")
        if self.modality not in {"peptide", "mrna", "dna"}:
            raise ValueError("Construct modality must be peptide, mrna or dna")
        object.__setattr__(self, "sequence_edits", tuple(self.sequence_edits))
        object.__setattr__(
            self, "chemical_modifications", tuple(self.chemical_modifications))
        if any(not isinstance(edit, ConstructSequenceEdit) for edit in self.sequence_edits):
            raise ValueError("Construct sequence edits must be typed records")
        if any(not isinstance(modification, ConstructChemicalModification)
               for modification in self.chemical_modifications):
            raise ValueError("Construct chemical modifications must be typed records")
        if self.chemical_modifications and self.modality != "peptide":
            raise ValueError("Chemical synthesis changes are not encoded mRNA/DNA edits")
        for modification in self.chemical_modifications:
            if modification.end > len(self.sequence):
                raise ValueError("Chemical modification lies outside the final sequence")
            if (modification.start == modification.end
                    and modification.start not in {0, len(self.sequence)}):
                raise ValueError("A zero-width chemical modification must name a terminus")
        if len(set(self.chemical_modifications)) != len(self.chemical_modifications):
            raise ValueError("Duplicate chemical modification")
        if self.mapping_status == "unresolved":
            if not self.mapping_reason:
                raise ValueError("An unresolved mapping requires a reason")
            if (self.native_start is not None or self.native_end is not None
                    or self.sequence_edits):
                raise ValueError("Unresolved mappings cannot assert native coordinates/edits")
            return
        if self.mapping_status != "resolved":
            raise ValueError("Unknown construct mapping status")
        if self.native_antigen is None:
            raise ValueError("A resolved construct requires its native antigen")
        _interval(self.native_start, self.native_end, "Native window")
        if self.native_start == self.native_end:
            raise ValueError("Native window must contain residues")
        if self.native_end > len(self.native_antigen.amino_acids):
            raise ValueError("Native window lies outside its antigen")
        edited, _ = self._apply_sequence_edits()
        if edited != self.sequence:
            raise ValueError("Recorded native window and edits do not reproduce final sequence")

    def _apply_sequence_edits(self):
        native = self.native_antigen.amino_acids
        cursor = self.native_start
        parts, offsets = [], []
        previous = None
        for edit in self.sequence_edits:
            if edit.start < cursor or edit.end > self.native_end:
                raise ValueError("Sequence edits must be ordered, disjoint and inside the window")
            if previous is not None and edit.start == previous.start:
                raise ValueError("Multiple edits at one native start are ambiguous")
            if native[edit.start:edit.end] == edit.replacement:
                raise ValueError("A sequence edit must change the native sequence")
            parts.extend((native[cursor:edit.start], edit.replacement))
            offsets.extend(range(cursor, edit.start))
            offsets.extend([None] * len(edit.replacement))
            cursor = edit.end
            previous = edit
        parts.append(native[cursor:self.native_end])
        offsets.extend(range(cursor, self.native_end))
        return "".join(parts), tuple(offsets)

    @property
    def native_sequence(self):
        if self.mapping_status != "resolved":
            return None
        return self.native_antigen.amino_acids[self.native_start:self.native_end]

    @property
    def native_residue_offsets(self):
        """Native coordinates, or None for added/replaced/unresolved residues."""
        if self.mapping_status != "resolved":
            return (None,) * len(self.sequence)
        return self._apply_sequence_edits()[1]

    @property
    def targetable_mask(self):
        """Map retained target residues and intact native deletion junctions."""
        if self.mapping_status != "resolved":
            return TargetableMask()
        offsets = self.native_residue_offsets
        native_mask = self.native_antigen.targetable_mask
        positions = [i for i, native in enumerate(offsets)
                     if native is not None and any(
                         interval.start <= native < interval.end
                         for interval in native_mask.intervals)]
        intervals = []
        for position in positions:
            if intervals and position == intervals[-1].end:
                intervals[-1] = AminoAcidInterval(intervals[-1].start, position + 1)
            else:
                intervals.append(AminoAcidInterval(position, position + 1))
        for interval in native_mask.intervals:
            if interval.start == interval.end:
                for i in range(1, len(offsets)):
                    if offsets[i - 1:i + 1] == (interval.start - 1, interval.start):
                        # An existing nonempty target interval can already cover
                        # this junction; do not create overlapping mask entries.
                        if not any(x.start <= i <= x.end for x in intervals):
                            intervals.append(AminoAcidInterval(i, i))
        return TargetableMask(tuple(sorted(intervals, key=lambda x: (x.start, x.end))))

    @property
    def removed_target_intervals(self):
        """Native target intervals not fully retained, including cropped targets."""
        if self.mapping_status != "resolved":
            return None
        if self.native_antigen is None:
            return ()
        offsets = self.native_residue_offsets
        retained = set(offsets) - {None}
        removed = []
        for interval in self.native_antigen.targetable_mask.intervals:
            if interval.start == interval.end:
                preserved = any(
                    offsets[i - 1:i + 1] == (interval.start - 1, interval.start)
                    for i in range(1, len(offsets)))
            else:
                preserved = set(range(interval.start, interval.end)) <= retained
            if not preserved:
                removed.append(interval)
        return tuple(removed)

    @property
    def modification_boundaries(self):
        """Internal final-sequence boundaries created by the recorded edits."""
        if self.mapping_status != "resolved":
            return ()
        boundaries = set()
        shift = -self.native_start
        for edit in self.sequence_edits:
            start = edit.start + shift
            boundaries.update((start, start + len(edit.replacement)))
            shift += len(edit.replacement) - (edit.end - edit.start)
        return tuple(sorted(b for b in boundaries if 0 < b < len(self.sequence)))

    @property
    def cache_identity(self):
        """Content identity for this full record, not merely its bare sequence.

        Predictor caches must additionally include model/settings provenance.
        """
        from .native_serialization import to_native_json

        canonical = json.dumps(json.loads(to_native_json(self)), sort_keys=True,
                               separators=(",", ":"), allow_nan=False)
        return "construct-v1:" + hashlib.sha256(canonical.encode("utf8")).hexdigest()

    @property
    def sequence_prediction_unassessed_reason(self):
        if self.mapping_status != "resolved":
            return "unresolved_native_mapping"
        if self.chemical_modifications:
            return "unsupported_chemical_modifications"
        return None

    def assessment_antigen(self):
        """Create final-context target/self input without changing the native antigen.

        This is an assessment view, not permission to manufacture a changed
        construct. Lost target content cannot silently inherit construct admission.
        """
        reason = self.sequence_prediction_unassessed_reason
        if reason:
            raise ValueError(f"Construct sequence prediction is unassessed: {reason}")
        attestation = self.native_antigen.tumor_specificity
        if self.removed_target_intervals:
            attestation = replace(
                attestation, status="held_out", requires_review=True,
                rationale_code="construct_removed_targetable_content")
        return replace(
            self.native_antigen,
            amino_acids=self.sequence,
            targetable_mask=self.targetable_mask,
            tumor_specificity=attestation,
            source_identifier=self.cache_identity,
        )


@dataclass(frozen=True)
class ConstructPlacement(DataclassSerializable):
    """An explicit occurrence of an antigen construct in a larger final product.

    Offset is in final product amino acids, also for an encoded mRNA/DNA ORF.
    Linkers and other product elements are not assigned to the native antigen.
    """

    construct: ConstructSequence
    offset: int

    def __post_init__(self):
        if not isinstance(self.construct, ConstructSequence):
            raise ValueError("Construct placement requires a ConstructSequence")
        _interval(self.offset, self.offset, "Construct placement")


def validate_construct_placements(placements, sequence, modality):
    """Validate exact emitted occurrences; never resolve repeated sequences by find."""
    previous_end = 0
    for placement in placements:
        if not isinstance(placement, ConstructPlacement):
            raise ValueError("Construct placements must be typed records")
        record = placement.construct
        if record.modality != modality:
            raise ValueError("Construct provenance modality differs from emitted product")
        end = placement.offset + len(record.sequence)
        if placement.offset < previous_end:
            raise ValueError("Construct placements must be ordered and nonoverlapping")
        if sequence[placement.offset:end] != record.sequence:
            raise ValueError("Construct provenance does not match the emitted sequence")
        for chemical in record.chemical_modifications:
            if chemical.start == chemical.end:
                terminal_offset = placement.offset + chemical.start
                if terminal_offset not in {0, len(sequence)}:
                    raise ValueError("A declared chemical terminus is internal to the final product")
        previous_end = end
