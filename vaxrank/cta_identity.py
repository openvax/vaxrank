"""Persisted annotation-scoped CTA identities from the public OncoRef API."""

from dataclasses import asdict, dataclass
import hashlib
import json

from serializable import DataclassSerializable


def identity_digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@dataclass(frozen=True)
class CTAAnnotationReference(DataclassSerializable):
    """Full candidate exclusion evidence, including rejected source mappings.

    Canonical membership and annotation source exclusions have separate hashes.
    All validation uses saved facts so decoding never consults a live reference.
    """

    identity_contract_version: int
    annotation_species: str
    annotation_assembly: str
    annotation_release: int
    annotation_sha256: str
    annotation_hash_kind: str
    canonical_candidate_gene_ids: tuple[str, ...]
    canonical_candidate_gene_ids_sha256: str
    identities: tuple[dict, ...]
    self_reference_excluded_gene_ids: tuple[str, ...]
    self_reference_excluded_gene_ids_sha256: str

    def __post_init__(self):
        if self.identity_contract_version != 1:
            raise ValueError("Unsupported CTA identity contract")
        if (self.annotation_species != "homo_sapiens"
                or self.annotation_assembly not in {"GRCh37", "GRCh38"}
                or self.annotation_release < 1):
            raise ValueError("CTA identity requires an explicit human annotation")
        if len(self.annotation_sha256) != 64:
            raise ValueError("CTA annotation requires SHA-256 provenance")
        candidates = tuple(sorted(set(self.canonical_candidate_gene_ids)))
        excluded = tuple(sorted(set(self.self_reference_excluded_gene_ids)))
        if identity_digest(list(candidates)) != self.canonical_candidate_gene_ids_sha256:
            raise ValueError("CTA canonical membership hash disagrees with saved facts")
        if identity_digest(list(excluded)) != self.self_reference_excluded_gene_ids_sha256:
            raise ValueError("CTA source exclusions hash disagrees with saved facts")
        verified = set()
        seen = set()
        for record in self.identities:
            if record["source_gene_id"] in seen:
                raise ValueError("Duplicate CTA source identity")
            seen.add(record["source_gene_id"])
            for key in ("identity_contract_version", "annotation_species",
                        "annotation_assembly", "annotation_release"):
                if record[key] != getattr(self, key):
                    raise ValueError("CTA identity and annotation contracts disagree")
            if (record["status"] in {"canonical", "alias"}
                    and record["canonical_gene_id"] in candidates):
                verified.add(record["source_gene_id"])
        if verified != set(excluded):
            raise ValueError("CTA exclusions disagree with verified annotation identities")
        object.__setattr__(self, "canonical_candidate_gene_ids", candidates)
        object.__setattr__(self, "self_reference_excluded_gene_ids", excluded)

    @property
    def fingerprint(self):
        return identity_digest(asdict(self))


def snapshot_cta_annotation(genome):
    """Resolve once for admission, exclusions and the persisted audit."""
    from oncoref import cta_annotation_gene_identities
    from oncoref.cta import cta_unfiltered_gene_ids
    from oncoref.gene_identity import GENE_IDENTITY_CONTRACT_VERSION
    from .reference_proteome import ensembl_dataset_cache_identity, _protein_source_snapshot

    records = tuple(record.as_dict() for record in cta_annotation_gene_identities(genome))
    candidates = tuple(sorted(cta_unfiltered_gene_ids()))
    excluded = tuple(sorted(record["source_gene_id"] for record in records
                            if record["status"] in {"canonical", "alias"}
                            and record["canonical_gene_id"] in candidates))
    annotation_hash = ensembl_dataset_cache_identity(genome)
    hash_kind = "pyensembl-source-files"
    if annotation_hash is None:
        # Custom in-memory annotations have no source files. Fingerprint the
        # actual annotated proteins and loci, explicitly labeling that scope.
        annotation_hash = identity_digest({
            "identities": records,
            "proteins": [(sequence, [asdict(source) for source in sources])
                         for sequence, sources in _protein_source_snapshot(genome)],
        })
        hash_kind = "annotated-cta-loci-and-proteins"
    return CTAAnnotationReference(
        identity_contract_version=GENE_IDENTITY_CONTRACT_VERSION,
        annotation_species=genome.species.latin_name,
        annotation_assembly=genome.reference_name,
        annotation_release=genome.release,
        annotation_sha256=annotation_hash, annotation_hash_kind=hash_kind,
        canonical_candidate_gene_ids=candidates,
        canonical_candidate_gene_ids_sha256=identity_digest(list(candidates)),
        identities=records, self_reference_excluded_gene_ids=excluded,
        self_reference_excluded_gene_ids_sha256=identity_digest(list(excluded)),
    )
