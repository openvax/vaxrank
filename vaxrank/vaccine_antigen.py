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

"""Source-agnostic antigen targetability and exact-self result types."""

from dataclasses import dataclass, field
from typing import Any

from serializable import DataclassSerializable


ANTIGEN_KIND_MUTATION = "mutation"
ANTIGEN_KIND_CTA = "CTA"
ANTIGEN_KIND_ERV = "ERV"
ANTIGEN_KIND_SPLICE = "splice"
ANTIGEN_KIND_VIRAL = "viral"
ANTIGEN_KINDS = frozenset({
    ANTIGEN_KIND_MUTATION,
    ANTIGEN_KIND_CTA,
    ANTIGEN_KIND_ERV,
    ANTIGEN_KIND_SPLICE,
    ANTIGEN_KIND_VIRAL,
})

ATTESTATION_ADMITTED = "admitted"
ATTESTATION_HELD_OUT = "held_out"
ATTESTATION_OVERRIDDEN = "overridden"
ATTESTATION_STATUSES = frozenset({
    ATTESTATION_ADMITTED,
    ATTESTATION_HELD_OUT,
    ATTESTATION_OVERRIDDEN,
})
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _normalized_gene_id(gene_id: str) -> str:
    return str(gene_id).split(".")[0]


def _validate_amino_acids(sequence: str, label: str) -> None:
    if not sequence:
        raise ValueError(f"{label} amino-acid sequence is required")
    invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
    if invalid:
        raise ValueError(
            f"{label} contains non-canonical amino acids: {', '.join(invalid)}"
        )


@dataclass(frozen=True)
class AminoAcidInterval(DataclassSerializable):
    """Zero-based, half-open targetable interval.

    ``start == end`` is allowed for deletion junctions and retains the
    historical ``MutantProteinFragment.interval_overlaps_mutation`` behavior.
    """

    start: int
    end: int

    def __post_init__(self):
        if self.start < 0:
            raise ValueError("Amino-acid interval start must be non-negative")
        if self.end < self.start:
            raise ValueError("Amino-acid interval end must not precede start")

    def overlaps(self, start: int, end: int) -> bool:
        if start < 0 or end < start:
            raise ValueError("Query interval must be a valid half-open interval")
        return start < self.end and end > self.start


@dataclass(frozen=True)
class TargetableMask(DataclassSerializable):
    """Explicit amino-acid regions carrying tumor-targetable content."""

    intervals: tuple[AminoAcidInterval, ...] = field(default_factory=tuple)

    def __post_init__(self):
        intervals = tuple(self.intervals)
        if intervals != self.intervals:
            object.__setattr__(self, "intervals", intervals)
        sorted_intervals = tuple(sorted(
            intervals,
            key=lambda interval: (interval.start, interval.end),
        ))
        if sorted_intervals != intervals:
            raise ValueError("Targetable intervals must be sorted")
        for previous, current in zip(intervals, intervals[1:]):
            if current.start < previous.end:
                raise ValueError("Targetable intervals must not overlap")

    def overlaps(self, start: int, end: int) -> bool:
        return any(interval.overlaps(start, end) for interval in self.intervals)

    def validate_for_sequence(self, amino_acids: str) -> None:
        sequence_length = len(amino_acids)
        for interval in self.intervals:
            if interval.end > sequence_length:
                raise ValueError(
                    "Targetable interval extends beyond antigen amino-acid sequence"
                )


@dataclass(frozen=True)
class TumorSpecificityAttestation(DataclassSerializable):
    """Evidence-bearing decision about construct admission."""

    status: str
    evidence_kind: str
    evidence_source: str
    evidence_version: str = ""
    patient_specific: bool = False
    rationale_code: str = ""
    requires_review: bool = False
    override_reason: str = ""

    def __post_init__(self):
        if self.status not in ATTESTATION_STATUSES:
            raise ValueError(
                f"Unknown tumor-specificity status {self.status!r}; "
                f"expected one of {sorted(ATTESTATION_STATUSES)}"
            )
        if not self.evidence_kind or not self.evidence_source:
            raise ValueError("Tumor-specificity evidence kind and source are required")
        if self.status == ATTESTATION_OVERRIDDEN and not self.override_reason:
            raise ValueError("An overridden attestation requires an override reason")

    @property
    def admits_construct(self) -> bool:
        return self.status in {ATTESTATION_ADMITTED, ATTESTATION_OVERRIDDEN}


@dataclass(frozen=True)
class SelfReferenceSource(DataclassSerializable):
    """One non-excluded protein source for an exact peptide match."""

    gene_id: str
    transcript_id: str = ""
    protein_id: str = ""
    gene_name: str = ""

    def __post_init__(self):
        if not self.gene_id:
            raise ValueError("Self-reference source requires a gene ID")
        object.__setattr__(self, "gene_id", _normalized_gene_id(self.gene_id))


@dataclass(frozen=True)
class SelfReferenceMatch(DataclassSerializable):
    """Antigen-aware exact-self classification for one peptide.

    ``source_provenance_complete`` explicitly distinguishes the transitional
    set-membership index from the later provenance-complete index. Missing
    provenance is never represented as evidence that no source exists.
    """

    peptide: str
    occurs: bool
    antigen_kind: str
    excluded_gene_ids: tuple[str, ...] = field(default_factory=tuple)
    sources: tuple[SelfReferenceSource, ...] = field(default_factory=tuple)
    source_provenance_complete: bool = False
    genome_release: str = ""

    def __post_init__(self):
        if self.antigen_kind not in ANTIGEN_KINDS:
            raise ValueError(f"Unknown antigen kind {self.antigen_kind!r}")
        _validate_amino_acids(self.peptide, "Self-reference peptide")
        excluded_gene_ids = tuple(sorted({
            _normalized_gene_id(gene_id) for gene_id in self.excluded_gene_ids
        }))
        sources = tuple(self.sources)
        object.__setattr__(self, "excluded_gene_ids", excluded_gene_ids)
        object.__setattr__(self, "sources", sources)
        if sources and not self.occurs:
            raise ValueError("A non-occurring self match cannot have match sources")
        if self.source_provenance_complete and self.occurs != bool(sources):
            raise ValueError(
                "A provenance-complete match must agree with its source list"
            )
        excluded = set(excluded_gene_ids)
        if any(source.gene_id in excluded for source in sources):
            raise ValueError("Self-match sources cannot include excluded genes")


@dataclass(frozen=True)
class VaccineAntigen(DataclassSerializable):
    """Source-agnostic antigen content, targetability, and admission evidence."""

    kind: str
    amino_acids: str
    targetable_mask: TargetableMask
    tumor_specificity: TumorSpecificityAttestation
    self_reference_excluded_gene_ids: tuple[str, ...] = field(default_factory=tuple)
    gene_name: str = ""
    gene_id: str = ""
    transcript_ids: tuple[str, ...] = field(default_factory=tuple)
    protein_ids: tuple[str, ...] = field(default_factory=tuple)
    species: str = ""
    source_identifier: str = ""

    def __post_init__(self):
        if self.kind not in ANTIGEN_KINDS:
            raise ValueError(
                f"Unknown antigen kind {self.kind!r}; expected one of "
                f"{sorted(ANTIGEN_KINDS)}"
            )
        _validate_amino_acids(self.amino_acids, "Vaccine antigen")
        if (
            self.tumor_specificity.admits_construct
            and not self.targetable_mask.intervals
        ):
            raise ValueError("An admitted antigen requires targetable content")
        self.targetable_mask.validate_for_sequence(self.amino_acids)
        excluded_gene_ids = tuple(sorted({
            _normalized_gene_id(gene_id)
            for gene_id in self.self_reference_excluded_gene_ids
        }))
        object.__setattr__(
            self, "self_reference_excluded_gene_ids", excluded_gene_ids
        )
        if self.gene_id:
            object.__setattr__(self, "gene_id", _normalized_gene_id(self.gene_id))
        object.__setattr__(self, "transcript_ids", tuple(self.transcript_ids))
        object.__setattr__(self, "protein_ids", tuple(self.protein_ids))

    def interval_is_targetable(self, start: int, end: int) -> bool:
        return self.targetable_mask.overlaps(start, end)

    @classmethod
    def from_mutant_protein_fragment(cls, fragment: Any) -> "VaccineAntigen":
        """Represent a legacy mutation fragment without changing its policy."""
        transcripts = tuple(
            getattr(fragment, "supporting_reference_transcripts", ()) or ()
        )
        gene_ids = tuple(sorted({
            _normalized_gene_id(getattr(transcript, "gene_id", ""))
            for transcript in transcripts
            if getattr(transcript, "gene_id", "")
        }))
        transcript_ids = tuple(sorted({
            str(getattr(transcript, "transcript_id", ""))
            for transcript in transcripts
            if getattr(transcript, "transcript_id", "")
        }))
        protein_ids = tuple(sorted({
            str(getattr(transcript, "protein_id", ""))
            for transcript in transcripts
            if getattr(transcript, "protein_id", "")
        }))
        variant = getattr(fragment, "variant", None)
        source_identifier = str(
            getattr(variant, "short_description", "") or variant or ""
        )
        return cls(
            kind=ANTIGEN_KIND_MUTATION,
            amino_acids=fragment.amino_acids,
            targetable_mask=TargetableMask((AminoAcidInterval(
                fragment.mutant_amino_acid_start_offset,
                fragment.mutant_amino_acid_end_offset,
            ),)),
            tumor_specificity=TumorSpecificityAttestation(
                status=ATTESTATION_ADMITTED,
                evidence_kind="somatic_variant",
                evidence_source="vaxrank mutation input",
                patient_specific=True,
                rationale_code="legacy_mutation_admission",
            ),
            gene_name=str(getattr(fragment, "gene_name", "") or ""),
            gene_id=gene_ids[0] if len(gene_ids) == 1 else "",
            transcript_ids=transcript_ids,
            protein_ids=protein_ids,
            source_identifier=source_identifier,
        )

    def self_reference_match(
        self,
        peptide: str,
        occurs: bool,
        *,
        genome_release: str = "",
    ) -> SelfReferenceMatch:
        """Build an explicit membership-only result for current indexes."""
        return SelfReferenceMatch(
            peptide=peptide,
            occurs=occurs,
            antigen_kind=self.kind,
            excluded_gene_ids=self.self_reference_excluded_gene_ids,
            source_provenance_complete=False,
            genome_release=genome_release,
        )
