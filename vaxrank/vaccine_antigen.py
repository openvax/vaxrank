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

import hashlib
import math
from dataclasses import dataclass, field
from typing import Any, Optional

from serializable import DataclassSerializable

from .amino_acids import validate_amino_acid_sequence
from .antigen_origin import AntigenOrigin, CodingReference, VariantOrigin
from .file_digest import is_sha256_digest
from .identifiers import normalize_ensembl_gene_id
from .normal_observations import MatchedNormalObservation, validate_normal_observations
from .rna_support import RNAFeatureEvidence, validate_rna_support, slice_rna_support


from .config.defaults import (
    ANTIGEN_KIND_CTA as ANTIGEN_KIND_CTA,
    ANTIGEN_KIND_ERV as ANTIGEN_KIND_ERV,
    ANTIGEN_KIND_FUSION as ANTIGEN_KIND_FUSION,
    ANTIGEN_KIND_ISOFORM as ANTIGEN_KIND_ISOFORM,
    ANTIGEN_KIND_MUTATION as ANTIGEN_KIND_MUTATION,
    ANTIGEN_KIND_SPLICE as ANTIGEN_KIND_SPLICE,
    ANTIGEN_KIND_SV as ANTIGEN_KIND_SV,
    ANTIGEN_KIND_VIRAL as ANTIGEN_KIND_VIRAL,
    ANTIGEN_KINDS as ANTIGEN_KINDS,
)

ATTESTATION_ADMITTED = "admitted"
ATTESTATION_HELD_OUT = "held_out"
ATTESTATION_OVERRIDDEN = "overridden"
ATTESTATION_STATUSES = frozenset({
    ATTESTATION_ADMITTED,
    ATTESTATION_HELD_OUT,
    ATTESTATION_OVERRIDDEN,
})

SOURCE_PROTEIN_COMPLETE = "complete"
SOURCE_PROTEIN_PARTIAL = "partial"
SOURCE_PROTEIN_UNASSESSED = "unassessed"
SOURCE_PROTEIN_COMPLETENESS_VALUES = frozenset({
    SOURCE_PROTEIN_COMPLETE,
    SOURCE_PROTEIN_PARTIAL,
    SOURCE_PROTEIN_UNASSESSED,
})
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

    def sliced(self, start: int, end: int) -> "TargetableMask":
        """Clip targetable intervals to a sequence slice and rebase them."""
        if start < 0 or end < start:
            raise ValueError("Sequence slice must be a valid half-open interval")
        intervals = []
        for interval in self.intervals:
            if not interval.overlaps(start, end):
                continue
            clipped_start = max(interval.start, start) - start
            clipped_end = min(interval.end, end) - start
            intervals.append(AminoAcidInterval(clipped_start, clipped_end))
        return TargetableMask(tuple(intervals))


@dataclass(frozen=True)
class TumorSpecificityEvidence(DataclassSerializable):
    """One machine-readable input to a tumor-specificity attestation."""

    evidence_kind: str
    evidence_source: str
    evidence_version: str = ""
    subject_id: str = ""
    sample_id: str = ""
    patient_specific: bool = False
    passed: Optional[bool] = None
    numeric_value: Optional[float] = None
    unit: str = ""
    threshold: Optional[float] = None
    details: tuple[tuple[str, str], ...] = field(default_factory=tuple)

    def __post_init__(self):
        if not self.evidence_kind or not self.evidence_source:
            raise ValueError("Tumor-specificity evidence kind and source are required")
        for label in ("numeric_value", "threshold"):
            value = getattr(self, label)
            if value is not None and not math.isfinite(value):
                raise ValueError(f"Tumor-specificity {label} must be finite")
        if (self.numeric_value is not None or self.threshold is not None) and not self.unit:
            raise ValueError("Quantitative tumor-specificity evidence requires units")
        details = tuple(sorted({
            (str(name), str(value)) for name, value in self.details
            if str(name)
        }))
        object.__setattr__(self, "details", details)


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
    evidence_records: tuple[TumorSpecificityEvidence, ...] = field(
        default_factory=tuple
    )

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
        records = tuple(self.evidence_records)
        if records and self.patient_specific != any(
            record.patient_specific for record in records
        ):
            raise ValueError(
                "Attestation patient-specific flag disagrees with evidence records"
            )
        object.__setattr__(self, "evidence_records", records)

    @property
    def admits_construct(self) -> bool:
        return self.status in {ATTESTATION_ADMITTED, ATTESTATION_OVERRIDDEN}


@dataclass(frozen=True)
class SelfReferenceSource(DataclassSerializable):
    """One annotated sequence source, not proof of expression/presentation.

    Biotype preserves distinctions such as predicted NMD and reference LoF.
    An empty value denotes missing annotation, including older saved sources.
    """

    gene_id: str
    transcript_id: str = ""
    protein_id: str = ""
    gene_name: str = ""
    species: str = ""
    transcript_biotype: str = ""
    orf_id: str = ""
    reference_entry: str = ""
    reference_occurrence: int | None = None

    def __post_init__(self):
        if not self.gene_id:
            raise ValueError("Self-reference source requires a gene ID")
        if self.reference_occurrence is not None and (
                type(self.reference_occurrence) is not int or self.reference_occurrence < 0
                or not self.reference_entry):
            raise ValueError('Self-reference occurrence requires an entry and nonnegative index')
        object.__setattr__(
            self, "gene_id", normalize_ensembl_gene_id(self.gene_id)
        )


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
    reference_sha256: str = ""

    def __post_init__(self):
        if self.antigen_kind not in ANTIGEN_KINDS:
            raise ValueError(f"Unknown antigen kind {self.antigen_kind!r}")
        validate_amino_acid_sequence(self.peptide, "Self-reference peptide")
        excluded_gene_ids = tuple(sorted({
            normalize_ensembl_gene_id(gene_id)
            for gene_id in self.excluded_gene_ids
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
        if any(source.gene_id in excluded and not (self.reference_sha256 and source.orf_id)
               for source in sources):
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
    source_metadata: tuple[tuple[str, str], ...] = field(default_factory=tuple)
    rna_support: tuple[RNAFeatureEvidence, ...] = field(default_factory=tuple)
    source_protein_offset: int = 0
    origin: AntigenOrigin | None = None
    source_protein_completeness: str = SOURCE_PROTEIN_UNASSESSED
    source_protein_sha256: str = ""
    source_protein_length: Optional[int] = None
    normal_observations: tuple[MatchedNormalObservation, ...] = field(default_factory=tuple)

    def __post_init__(self):
        normalized = next((kind for kind in ANTIGEN_KINDS if kind.lower() == self.kind.lower()), self.kind)
        object.__setattr__(self, 'kind', normalized)
        object.__setattr__(
            self, "normal_observations", validate_normal_observations(self.normal_observations)
        )
        if self.origin is not None and not isinstance(self.origin, AntigenOrigin):
            raise TypeError('Antigen origin must be an AntigenOrigin record')
        if type(self.source_protein_offset) is not int or self.source_protein_offset < 0:
            raise ValueError('Source protein offset must be a nonnegative integer')
        if self.source_protein_completeness not in SOURCE_PROTEIN_COMPLETENESS_VALUES:
            raise ValueError(
                f"Unknown source_protein_completeness {self.source_protein_completeness!r}; "
                f"expected one of {sorted(SOURCE_PROTEIN_COMPLETENESS_VALUES)}"
            )
        object.__setattr__(self, "rna_support", validate_rna_support(self.rna_support, len(self.amino_acids)))
        if self.kind not in ANTIGEN_KINDS:
            raise ValueError(
                f"Unknown antigen kind {self.kind!r}; expected one of "
                f"{sorted(ANTIGEN_KINDS)}"
            )
        validate_amino_acid_sequence(self.amino_acids, "Vaccine antigen")
        if self.source_protein_completeness == SOURCE_PROTEIN_COMPLETE:
            if not self.source_protein_sha256 and self.source_protein_length is None:
                if self.source_protein_offset != 0:
                    raise ValueError(
                        "Complete source protein at nonzero offset requires explicit "
                        "source_protein_sha256 and source_protein_length"
                    )
                digest = hashlib.sha256(self.amino_acids.encode("ascii")).hexdigest()
                object.__setattr__(self, "source_protein_sha256", digest)
                object.__setattr__(self, "source_protein_length", len(self.amino_acids))
            else:
                if not is_sha256_digest(self.source_protein_sha256):
                    raise ValueError("Complete source protein requires a SHA-256 digest")
                if (
                    type(self.source_protein_length) is not int
                    or self.source_protein_length < self.source_protein_offset + len(self.amino_acids)
                ):
                    raise ValueError(
                        "Complete source protein length must be an integer covering the window"
                    )
        elif self.source_protein_sha256 or self.source_protein_length is not None:
            raise ValueError(
                "Only complete source proteins carry source_protein_sha256 and source_protein_length"
            )
        if (
            self.tumor_specificity.admits_construct
            and not self.targetable_mask.intervals
        ):
            raise ValueError("An admitted antigen requires targetable content")
        self.targetable_mask.validate_for_sequence(self.amino_acids)
        excluded_gene_ids = tuple(sorted({
            normalize_ensembl_gene_id(gene_id)
            for gene_id in self.self_reference_excluded_gene_ids
        }))
        object.__setattr__(
            self, "self_reference_excluded_gene_ids", excluded_gene_ids
        )
        if self.gene_id:
            object.__setattr__(
                self, "gene_id", normalize_ensembl_gene_id(self.gene_id)
            )
        object.__setattr__(self, "transcript_ids", tuple(self.transcript_ids))
        object.__setattr__(self, "protein_ids", tuple(self.protein_ids))
        source_metadata = tuple(sorted({
            (str(name), str(value))
            for name, value in self.source_metadata
            if str(name) and str(value)
        }))
        object.__setattr__(self, "source_metadata", source_metadata)

    def interval_is_targetable(self, start: int, end: int) -> bool:
        return self.targetable_mask.overlaps(start, end)

    def targetable_span(self) -> AminoAcidInterval:
        """Smallest contiguous interval containing all targetable residues."""
        intervals = self.targetable_mask.intervals
        if not intervals:
            raise ValueError("Vaccine antigen has no targetable content")
        return AminoAcidInterval(intervals[0].start, intervals[-1].end)

    def sliced(self, start: int, end: int) -> "VaccineAntigen":
        """Return a peptide-window antigen with target coordinates rebased."""
        if start < 0 or end > len(self.amino_acids) or end <= start:
            raise ValueError("Antigen slice must be non-empty and in bounds")
        targetable_mask = self.targetable_mask.sliced(start, end)
        if self.tumor_specificity.admits_construct and not targetable_mask.intervals:
            raise ValueError("Antigen slice does not overlap targetable content")
        return VaccineAntigen(
            kind=self.kind,
            amino_acids=self.amino_acids[start:end],
            targetable_mask=targetable_mask,
            tumor_specificity=self.tumor_specificity,
            self_reference_excluded_gene_ids=self.self_reference_excluded_gene_ids,
            gene_name=self.gene_name,
            gene_id=self.gene_id,
            transcript_ids=self.transcript_ids,
            protein_ids=self.protein_ids,
            species=self.species,
            source_identifier=self.source_identifier,
            source_metadata=self.source_metadata,
            rna_support=slice_rna_support(self.rna_support, start, end),
            source_protein_offset=self.source_protein_offset + start,
            origin=self.origin,
            source_protein_completeness=self.source_protein_completeness,
            source_protein_sha256=self.source_protein_sha256,
            source_protein_length=self.source_protein_length,
            normal_observations=self.normal_observations,
        )

    @classmethod
    def from_assembled_mutation_sequence(
            cls,
            amino_acids: str,
            targetable_intervals: tuple[AminoAcidInterval, ...],
            tumor_specificity: TumorSpecificityAttestation,
            gene_name: str = "",
            gene_id: str = "",
            transcript_ids: tuple[str, ...] = (),
            protein_ids: tuple[str, ...] = (),
            species: str = "",
            source_identifier: str = "",
            source_metadata: tuple[tuple[str, str], ...] = (),
            source_protein_completeness: str = SOURCE_PROTEIN_UNASSESSED,
    ) -> "VaccineAntigen":
        """Represent a sample-specific RNA translation of mutation(s).

        This constructor is for compound, phased, or otherwise assembled
        mutation products whose protein coordinates cannot be represented by
        one Varcode effect. It intentionally carries no inferred wild-type
        comparator; callers must provide the translated sequence and every
        targetable amino-acid interval explicitly.
        """
        return cls(
            kind=ANTIGEN_KIND_MUTATION,
            amino_acids=amino_acids,
            targetable_mask=TargetableMask(tuple(targetable_intervals)),
            tumor_specificity=tumor_specificity,
            gene_name=gene_name,
            gene_id=gene_id,
            transcript_ids=transcript_ids,
            protein_ids=protein_ids,
            species=species,
            source_identifier=source_identifier,
            source_metadata=source_metadata,
            source_protein_completeness=source_protein_completeness,
        )

    @property
    def display_identifier(self) -> str:
        """Best stable source identifier for reports and construct labels."""
        identifiers = (
            (self.source_identifier,)
            + self.protein_ids
            + self.transcript_ids
            + ((self.gene_id,) if self.gene_id else ())
        )
        return next((value for value in identifiers if value), self.kind)

    @property
    def display_gene_name(self) -> str:
        """Best available gene label for reports and construct labels."""
        return self.gene_name or self.gene_id or "unknown"

    @classmethod
    def from_fusion_sequence(
            cls,
            amino_acids: str,
            junction_offset: int,
            tumor_specificity: TumorSpecificityAttestation,
            gene_name: str = "",
            transcript_ids: tuple[str, ...] = (),
            species: str = "",
            source_identifier: str = "",
            source_metadata: tuple[tuple[str, str], ...] = (),
            source_protein_completeness: str = SOURCE_PROTEIN_UNASSESSED,
    ) -> "VaccineAntigen":
        """Represent an assembled fusion translation and its exact junction.

        ``junction_offset`` is the zero-based boundary between the last amino
        acid from the left partner and the first amino acid from the right
        partner. A zero-width target interval makes an epitope targetable only
        when it contains residues from both sides of that boundary.
        """
        if not isinstance(junction_offset, int):
            raise TypeError("Fusion junction offset must be an integer")
        if junction_offset <= 0 or junction_offset >= len(amino_acids):
            raise ValueError(
                "Fusion junction offset must have amino acids on both sides"
            )
        metadata = tuple(source_metadata) + (
            ("junction_offset_aa", str(junction_offset)),
        )
        return cls(
            kind=ANTIGEN_KIND_FUSION,
            origin=AntigenOrigin(rna_product="fusion"),
            amino_acids=amino_acids,
            targetable_mask=TargetableMask((AminoAcidInterval(
                junction_offset, junction_offset,
            ),)),
            tumor_specificity=tumor_specificity,
            gene_name=gene_name,
            transcript_ids=transcript_ids,
            species=species,
            source_identifier=source_identifier,
            source_metadata=metadata,
            source_protein_completeness=source_protein_completeness,
        )

    @classmethod
    def from_mutant_protein_fragment(cls, fragment: Any) -> "VaccineAntigen":
        """Represent a legacy mutation fragment without changing its policy."""
        from varcode import StructuralVariant, Variant
        transcripts = tuple(
            getattr(fragment, "supporting_reference_transcripts", ()) or ()
        )
        gene_ids = tuple(sorted({
            normalize_ensembl_gene_id(getattr(transcript, "gene_id", ""))
            for transcript in transcripts
            if getattr(transcript, "gene_id", "")
        }))
        transcript_ids = tuple(sorted({
            str(getattr(transcript, "transcript_id", ""))
            for transcript in transcripts
            if getattr(transcript, "transcript_id", "")
        }))
        protein_ids = getattr(fragment, "supporting_reference_protein_ids", None)
        if protein_ids is None:
            protein_ids = tuple(sorted({
                str(getattr(transcript, "protein_id", ""))
                for transcript in transcripts
                if getattr(transcript, "protein_id", "")
            }))
        variant = getattr(fragment, "variant", None)
        source_identifier = str(
            getattr(variant, "short_description", "") or variant or ""
        )
        proxy_meta = ()
        proxy_code = getattr(fragment, "proxy_expression_cancer_code", "") or ""
        proxy_tpm = getattr(fragment, "proxy_expression_tpm", None)
        if proxy_code or proxy_tpm is not None:
            entries = []
            if proxy_code:
                entries.append(("proxy_expression_cancer_code", str(proxy_code)))
            if proxy_tpm is not None:
                entries.append(("proxy_expression_tpm", f"{float(proxy_tpm):.4f}"))
            proxy_source = getattr(fragment, "proxy_expression_source", "") or ""
            if proxy_source:
                entries.append(("proxy_expression_source", str(proxy_source)))
            proxy_meta = tuple(entries)
        return cls(
            kind=ANTIGEN_KIND_SV if isinstance(variant, StructuralVariant) else ANTIGEN_KIND_MUTATION,
            origin=AntigenOrigin(variants=(VariantOrigin.from_variant(variant),) if isinstance(variant, Variant) else (),
                coding_references=tuple(CodingReference(kind="annotated", source_transcript_id=t)
                                        for t in transcript_ids)),
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
            source_metadata=proxy_meta,
            rna_support=getattr(fragment, "rna_support", ()),
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
