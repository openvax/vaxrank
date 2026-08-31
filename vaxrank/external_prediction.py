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

"""Stable identities for prediction records imported from external tools.

External reports routinely contain the same peptide and flanking sequence for
different variants, genes, or transcripts. Sequence position is therefore not
an identity. :class:`ExternalPredictionKey` is the public join key shared by
loading, Topiary evaluation, audit reports, native exports, and construct
selection. Its JSON representation is deterministic, human-readable, and
round-trippable without access to the original report.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json

from . import cells

# Re-exported under their historical names. :mod:`vaxrank.cells` is the one
# normalization vocabulary; these aliases exist so the many external-report
# call sites keep reading in the domain's language.
external_text = cells.text
external_values = cells.values


# pVACseq's two flavors spell the same fact differently; topiary normalizes
# both onto one of these names. Rows only ever reach this module already
# normalized, so the raw pVACseq spellings are deliberately absent — listing
# them would imply a second vocabulary is still in play here.
_PVACSEQ_RNA_DEPTH_COLUMNS = ("rna_depth", "tumor_rna_depth")
_PVACSEQ_RNA_VAF_COLUMNS = ("rna_vaf", "tumor_rna_vaf")
_PVACSEQ_DNA_VAF_COLUMNS = ("dna_vaf", "tumor_dna_vaf")


def pvacseq_variant_id(row) -> str:
    """Return one stable variant identifier from either pVACseq flavor.

    Precedence is **not** arbitrary: it mirrors ``topiary.read_pvacseq``,
    which prefers pVACseq's own ``Index`` when the column exists and only
    falls back to composing ``chrom-start-ref-alt`` when it does not. A reader
    that chose coordinates first would build identities that never join
    against topiary-normalized ones for any file carrying both.

    ``variant`` / ``ID`` come first because those are the already-normalized
    spellings (topiary's output column and pVACseq's aggregated-flavor column
    respectively); they carry the same value ``Index`` would have produced.
    """
    normalized = cells.first_text(row, "variant", "ID", "Index")
    if normalized:
        return normalized
    chrom = cells.text(row.get("Chromosome"))
    start = cells.text(row.get("Start"))
    ref = cells.text(row.get("Reference"))
    alt = cells.text(row.get("Variant"))
    if all((chrom, start, ref, alt)):
        return f"{chrom}-{start}-{ref}-{alt}"
    return ""


def lens_variant_id(row) -> str:
    """Return a LENS variant identity including alleles when available."""
    coords = external_text(row.get("variant"))
    ref = (
        external_text(row.get("snv_ref_allele"))
        or external_text(row.get("indel_ref_allele"))
    )
    alt = (
        external_text(row.get("snv_alt_allele"))
        or external_text(row.get("indel_alt_allele"))
    )
    if ref or alt:
        return f"{coords}:{ref}>{alt}"
    return coords


@dataclass(frozen=True)
class ExternalPredictionKey:
    """Complete provenance and peptide position for one score group.

    Allele and predictor are deliberately absent: they are leaves within one
    candidate score group. Variant, gene, transcript, source context, peptide,
    and offset are present because none of those fields is individually unique.
    """

    source_format: str
    schema_version: int = 1
    variant_id: str = ""
    antigen_source: str = ""
    gene_ids: tuple[str, ...] = ()
    gene_names: tuple[str, ...] = ()
    transcript_ids: tuple[str, ...] = ()
    # The row's own gene / transcript, before the multi-value columns are
    # folded in. The sets above are sorted so that two rows listing the same
    # genes in different orders resolve to one identity — which makes
    # ``gene_names[0]`` alphabetical, not primary. Anything labelling a
    # construct must use these instead.
    primary_gene_name: str = ""
    primary_transcript_id: str = ""
    species: str = ""
    peptide: str = ""
    source_sequence: str = ""
    offset: int = 0

    def __post_init__(self):
        if not self.source_format:
            raise ValueError("External prediction source format is required")
        if self.schema_version != 1:
            raise ValueError(
                f"Unsupported external prediction schema {self.schema_version}")
        if not self.peptide:
            raise ValueError("External prediction peptide is required")
        if self.offset < 0:
            raise ValueError("External prediction offset cannot be negative")
        # The key is only a usable join target if its position fields describe
        # a real position. Enforcing it here means every producer is held to
        # the same invariant instead of each one re-checking (or forgetting).
        if self.source_sequence:
            window = self.source_sequence[
                self.offset:self.offset + len(self.peptide)]
            if window != self.peptide:
                raise ValueError(
                    "External prediction peptide %r is not at offset %d of "
                    "its source sequence" % (self.peptide, self.offset))
        for name in ("gene_ids", "gene_names", "transcript_ids"):
            values = tuple(sorted(set(getattr(self, name))))
            object.__setattr__(self, name, values)

    # Annotation carried alongside the key but deliberately outside its
    # identity. Which gene a report *named first* for a row is presentation,
    # not biology: two rows listing the same genes in different orders
    # describe the same candidate, and folding order into the identity would
    # split them — the exact fragility the sorted sets above prevent.
    NON_IDENTITY_FIELDS = ("primary_gene_name", "primary_transcript_id")

    def payload(self, include_position=True) -> dict:
        """Return the canonical JSON-compatible identity fields."""
        result = asdict(self)
        for name in self.NON_IDENTITY_FIELDS:
            result.pop(name, None)
        if not include_position:
            result.pop("peptide")
            result.pop("offset")
        return result

    @property
    def ordered_transcript_ids(self) -> tuple[str, ...]:
        """Transcript IDs with the row's own first.

        ``transcript_ids`` is sorted for identity stability, so its first
        element is alphabetical. Anything that treats one transcript as
        representative — construct annotation, varcode cross-checks — wants
        the one the report actually named for this row.
        """
        if not self.primary_transcript_id:
            return self.transcript_ids
        return (self.primary_transcript_id,) + tuple(
            t for t in self.transcript_ids if t != self.primary_transcript_id)

    @property
    def identifier(self) -> str:
        """Canonical, externally legible prediction-group identifier."""
        return json.dumps(
            self.payload(), sort_keys=True, separators=(",", ":"))

    @property
    def construct_identifier(self) -> str:
        """Identity of the source context used to build one construct.

        Multiple candidate epitopes from the same source window share this
        value. The peptide and offset are omitted while all biological source
        provenance remains.
        """
        return json.dumps(
            self.payload(include_position=False),
            sort_keys=True,
            separators=(",", ":"),
        )

    @classmethod
    def from_identifier(cls, identifier: str) -> "ExternalPredictionKey":
        """Reconstruct a key from :attr:`identifier`.

        The annotation fields in :data:`NON_IDENTITY_FIELDS` are absent from
        the identifier by design, so a reconstructed key carries the identity
        but not the report's own naming order.
        """
        payload = json.loads(identifier)
        for name in ("gene_ids", "gene_names", "transcript_ids"):
            payload[name] = tuple(payload.get(name) or ())
        return cls(**payload)

    @classmethod
    def from_lens_row(cls, row) -> "ExternalPredictionKey | None":
        """Build the exact score-group key for one raw LENS row."""
        peptide = external_text(row.get("peptide"))
        source_sequence = external_text(row.get("pep_context"))
        if not peptide:
            return None
        if source_sequence:
            offset = source_sequence.find(peptide)
            if offset < 0:
                return None
        else:
            source_sequence = peptide
            offset = 0
        return cls(
            source_format="lens",
            variant_id=lens_variant_id(row),
            antigen_source=external_text(row.get("antigen_source")),
            gene_ids=external_values(
                row.get("gene_id"),
                row.get("all_gene_ids_encoding_peptide"),
            ),
            gene_names=external_values(
                row.get("gene"),
                row.get("all_gene_names_encoding_peptide"),
            ),
            transcript_ids=external_values(
                row.get("transcript_id"),
                row.get("all_transcript_ids_encoding_peptide"),
            ),
            primary_gene_name=cells.text(row.get("gene")),
            primary_transcript_id=cells.text(row.get("transcript_id")),
            species=external_text(row.get("species")),
            peptide=peptide,
            source_sequence=source_sequence,
            offset=offset,
        )

    @classmethod
    def from_pvacseq_row(cls, row) -> "ExternalPredictionKey | None":
        """Build the exact score-group key for either pVACseq flavor."""
        peptide = (
            external_text(row.get("peptide"))
            or external_text(row.get("Best Peptide"))
            or external_text(row.get("MT Epitope Seq"))
        )
        if not peptide:
            return None
        return cls(
            source_format="pvacseq",
            variant_id=pvacseq_variant_id(row),
            # The aggregated and all-epitopes readers do not expose effect
            # type consistently. Variant/gene/transcript carry the stable
            # provenance; effect remains a report annotation, not identity.
            antigen_source="",
            gene_names=external_values(
                row.get("gene"), row.get("Gene"), row.get("Gene Name")),
            transcript_ids=external_values(
                row.get("transcript"),
                row.get("Best Transcript"),
                row.get("Transcript"),
            ),
            primary_gene_name=cells.first_text(
                row, "gene", "Gene", "Gene Name"),
            primary_transcript_id=cells.first_text(
                row, "transcript", "Best Transcript", "Transcript"),
            species=external_text(row.get("species")),
            peptide=peptide,
            # pVACseq ships no source-protein context, so the peptide is its
            # own window and its offset within that window is 0 by
            # definition. Carrying a source-protein offset here instead would
            # describe the peptide as sitting partway inside itself, which
            # the position invariant rejects — and rejecting it mid-load
            # would take down the whole run, where every other malformed-row
            # case skips the row.
            source_sequence=peptide,
            offset=0,
        )
