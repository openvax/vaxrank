"""Experimental consumption of checksummed ORFdata host proteomes."""

from collections import defaultdict
from functools import cached_property
import hashlib
import json
from pathlib import Path
import re

from .input_scope import InputScope, combine_scopes
from .reference_coverage import DISCOVER_CATALOGUES, reference_coverage, validate_catalogues
from .vaccine_antigen import SelfReferenceMatch, SelfReferenceSource


class HostReference:
    """One frozen host inventory, shared by candidate and product assessment.

    Build and select ORFs with ORFdata, then pass its Proteome or saved snapshot.
    Viral entries remain in the snapshot but do not count as host self.
    ``catalogues`` optionally asserts the exact source release labels recorded
    in its host ORF layers; an empty list asserts annotated-only coverage.
    ``coverage`` restores a saved coverage record (None means historical
    availability was not recorded); by default availability is captured offline.
    """

    def __init__(self, proteome, *, catalogues=None, coverage=DISCOVER_CATALOGUES):
        from orfdata import Proteome
        if not isinstance(proteome, Proteome):
            raise TypeError('host_reference requires an ORFdata Proteome or snapshot')
        self.proteome = proteome
        self.identity = proteome.identity
        validate_catalogues(catalogues, proteome)
        self._coverage = reference_coverage(proteome, coverage)
        reference = proteome.reference
        self.scope = InputScope(
            host_taxon_id=reference.taxon_id, reference_assembly=reference.genome_build,
            annotation=(f'{reference.annotation_name}:{reference.annotation_version}'
                        if reference.annotation_version is not None else None),
            host_reference_sha256=self.identity)
        sources = defaultdict(list)
        for entry in proteome:
            if entry.kind not in ('annotated', 'extra_orf', 'viral'):
                raise ValueError(f'Unsupported host reference entry kind: {entry.kind!r}')
            if entry.kind == 'viral':
                continue
            if not entry.occurrences or any(not occurrence.gene_id for occurrence in entry.occurrences):
                raise ValueError(f'Host reference entry {entry.accession!r} lacks source occurrences')
            for index, occurrence in enumerate(entry.occurrences):
                sources[entry.sequence].append(SelfReferenceSource(
                    gene_id=occurrence.gene_id, transcript_id=occurrence.transcript_id,
                    protein_id=entry.accession if entry.kind == 'annotated' else '',
                    gene_name=occurrence.gene_name or '', species=reference.species,
                    transcript_biotype=occurrence.transcript_biotype or '',
                    orf_id=entry.accession if entry.kind == 'extra_orf' else '',
                    reference_entry=entry.accession, reference_occurrence=index))
        self._sources = tuple((sequence, tuple(records)) for sequence, records in sources.items())

    @classmethod
    def load(cls, value, *, coverage=DISCOVER_CATALOGUES):
        """Accept a path, a pinned path mapping, a Proteome, or this wrapper."""
        if isinstance(value, cls):
            if coverage is not DISCOVER_CATALOGUES and coverage != value.coverage:
                raise ValueError('Cannot replace frozen host catalogue coverage')
            return value
        from orfdata import Proteome
        if isinstance(value, Proteome):
            return cls(value, coverage=coverage)
        expected, catalogues = None, None
        if isinstance(value, dict):
            if not {'path', 'identity'} <= set(value) or set(value) - {'path', 'identity', 'catalogues'}:
                raise ValueError('host_reference mapping requires path and identity, with optional catalogues')
            catalogues = value.get('catalogues')
            if 'catalogues' in value and catalogues is None:
                raise ValueError('host_reference catalogues must be an explicit list')
            expected, value = value['identity'], value['path']
            if not isinstance(expected, str) or not re.fullmatch('[0-9a-f]{64}', expected):
                raise ValueError('host_reference identity must be a SHA-256')
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise ValueError('host_reference requires an ORFdata snapshot path')
        result = cls(Proteome.load(value), catalogues=catalogues, coverage=coverage)
        if expected is not None and result.identity != expected:
            raise ValueError('Host reference snapshot differs from its declared identity')
        return result

    @property
    def coverage(self):
        """Detached, frozen catalogue availability and actual inventory coverage."""
        from copy import deepcopy
        return deepcopy(self._coverage)

    def validate(self, scope, *, genome=None, location='Host reference'):
        """Reject incompatible host, annotation or previously saved identity."""
        scope = combine_scopes(scope, self.scope, location)
        if genome is not None:
            self.proteome.reference.validate_genome(genome)
        return scope

    def matches(self, peptides, antigen):
        """Keep all mapped occurrences, including extra ORFs at excluded genes."""
        self.validate_antigen(antigen, native_sequence=False)
        peptides = tuple(dict.fromkeys(peptides))
        if not peptides:
            return {}
        sources = {peptide: [] for peptide in peptides}
        excluded = set(antigen.self_reference_excluded_gene_ids)
        for sequence, records in self._sources:
            retained = tuple(r for r in records if r.orf_id or r.gene_id not in excluded)
            for peptide in peptides:
                if peptide in sequence:
                    sources[peptide].extend(retained)
        return {peptide: SelfReferenceMatch(
            peptide=peptide, occurs=bool(records), antigen_kind=antigen.kind,
            excluded_gene_ids=antigen.self_reference_excluded_gene_ids,
            sources=tuple(records), source_provenance_complete=True,
            genome_release=str(self.proteome.reference.annotation_version or ''),
            reference_sha256=self.identity) for peptide, records in sources.items()}

    @cached_property
    def _sequence_hashes(self):
        sequences = defaultdict(set)
        for sequence, _ in self._sources:
            sequences[hashlib.sha256(sequence.encode()).hexdigest()].add(sequence)
        return sequences

    def validate_antigen(self, antigen, *, native_sequence=True):
        """Check frozen restriction evidence against the actual complete host."""
        records = [r for r in antigen.tumor_specificity.evidence_records if r.evidence_kind == 'cta_rna_restriction']
        expected = dict(antigen.source_metadata).get('cta_restriction_sha256')
        if not records and expected is None:
            return
        if len(records) != 1 or expected is None:
            raise ValueError('CTA antigen lost its frozen restriction evidence')
        from .cta_restriction import CTARestrictionAssessment, host_occurrences
        from .native_serialization import from_native_json
        details = dict(records[0].details)
        if set(details) != {'assessment'} or hashlib.sha256(details['assessment'].encode()).hexdigest() != expected:
            raise ValueError('CTA antigen restriction evidence checksum differs')
        assessment = from_native_json(details['assessment'], CTARestrictionAssessment)
        native = json.loads(assessment.native_evidence_json)
        # Final-product audit views can include leaders/linkers and edited flanks.
        # Their restriction evidence remains bound to the original source protein;
        # ConstructSequence separately validates the actual placement and edits.
        digest = dict(antigen.source_metadata).get('protein_sequence_sha256')
        sequences = self._sequence_hashes.get(digest, ())
        source_sequence = next(iter(sequences)) if len(sequences) == 1 else None
        if (native['taxon_id'] != self.scope.host_taxon_id
                or source_sequence is None
                or (native_sequence and antigen.amino_acids not in source_sequence)
                or set(assessment.occurrences) != set(host_occurrences(self, source_sequence))
                or (assessment.evidence is not None and assessment.evidence.host_reference_sha256 != self.identity)
                or (antigen.tumor_specificity.admits_construct and not assessment.passed)
                or antigen.tumor_specificity.status == 'overridden'):
            raise ValueError('CTA antigen restriction differs from its complete host reference')

    def save(self, path):
        """Use ORFdata's lossless snapshot, including layers and exclusions."""
        identity = self.proteome.save(path)
        if identity != self.identity:
            raise ValueError('Host reference changed after validation')
        return identity


def combine_host_references(references):
    """Resolve one shared reference without input-order overrides."""
    result = None
    seen = set()
    for reference in references:
        if reference is None:
            continue
        if isinstance(reference, (str, Path)):
            key = str(Path(reference).resolve())
        elif (isinstance(reference, dict) and set(reference) == {'path', 'identity'}
              and isinstance(reference['path'], (str, Path)) and isinstance(reference['identity'], str)):
            key = (str(Path(reference['path']).resolve()), reference['identity'])
        else:
            key = id(reference)
        if key in seen:
            continue
        seen.add(key)
        reference = HostReference.load(reference)
        if result is not None and result.identity != reference.identity:
            raise ValueError('Inputs contain conflicting host_reference snapshots')
        if result is not None and result.coverage != reference.coverage:
            raise ValueError('Inputs contain conflicting frozen host catalogue coverage; reuse the same saved reference')
        result = reference
    return result
