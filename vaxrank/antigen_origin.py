"""Shared biological origin, independent of input format and targetability.

Coordinates are the producer's one-based genomic allele coordinates. RNA
products do not imply a causal DNA event. Missing facts remain unspecified.
"""
from dataclasses import asdict, dataclass

from serializable import DataclassSerializable


VARIANT_TYPES = frozenset({'SNV', 'MNV', 'INDEL', 'SV'})
SV_TYPES = frozenset({'DEL', 'DUP', 'INV', 'INS', 'BND', 'CNV'})


@dataclass(frozen=True)
class VariantOrigin(DataclassSerializable):
    variant_type: str
    sv_type: str = ''
    reference_assembly: str = ''
    contig: str = ''
    start: int | None = None
    ref: str = ''
    alt: str = ''
    end: int | None = None
    mate_contig: str = ''
    mate_start: int | None = None
    mate_orientation: str = ''
    event_id: str = ''

    def __post_init__(self):
        object.__setattr__(self, 'variant_type', self.variant_type.upper())
        object.__setattr__(self, 'sv_type', self.sv_type.upper())
        if self.variant_type not in VARIANT_TYPES:
            raise ValueError('Unknown variant type: ' + self.variant_type)
        if self.sv_type and self.variant_type != 'SV':
            raise ValueError('SV subtype requires variant_type=SV')
        for name in ('start', 'end', 'mate_start'):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or value < 1):
                raise ValueError(name + ' must be a positive one-based genomic coordinate')

    @classmethod
    def from_variant(cls, variant):
        """Use Varcode's normalized allele classification and coordinates."""
        structural = getattr(variant, 'is_structural', False)
        return cls(
            variant_type='SV' if structural else 'SNV' if variant.is_snv else
                         'INDEL' if variant.is_indel else 'MNV',
            sv_type=variant.sv_type if structural else '',
            reference_assembly=variant.genome.reference_name,
            contig=variant.contig, start=variant.start, ref=variant.ref, alt=variant.alt,
            **({name: getattr(variant, name) for name in
                ('end', 'mate_contig', 'mate_start', 'mate_orientation')
                if getattr(variant, name) is not None} if structural else {}))

    @classmethod
    def from_alleles(cls, value):
        from varcode import Variant
        # Isovar's contract uses normalized genomic alleles. Let Varcode own
        # trimming/classification; reference identity stays the producer's.
        variant = Variant(value['contig'], value['start'], value['ref'], value['alt'])
        return cls(variant_type='SNV' if variant.is_snv else 'INDEL' if variant.is_indel else 'MNV',
            reference_assembly=value.get('reference_name') or value.get('reference_assembly') or '',
            contig=variant.contig, start=variant.start, ref=variant.ref, alt=variant.alt)


@dataclass(frozen=True)
class CodingReference(DataclassSerializable):
    kind: str | None = None
    orf_type: str | None = None
    orf_id: str | None = None
    source: str | None = None
    reference_key: str | None = None
    source_transcript_id: str | None = None

    def __post_init__(self):
        if self.kind is not None:
            object.__setattr__(self, 'kind', self.kind.lower())


@dataclass(frozen=True)
class GenomicEndpoint(DataclassSerializable):
    contig: str = ''
    position: int | None = None
    strand: str = ''
    description: str = ''

    def __post_init__(self):
        if self.position is not None and (type(self.position) is not int or self.position < 1):
            raise ValueError('Endpoint position must be a positive one-based coordinate')
        if self.strand not in ('', '+', '-'):
            raise ValueError('Endpoint strand must be +, - or unspecified')

    @classmethod
    def from_record(cls, value):
        if isinstance(value, (list, tuple)) and len(value) == 3:
            return cls(str(value[0]), value[1], str(value[2]))
        if isinstance(value, dict):
            return cls(value.get('contig') or '', value.get('position'), value.get('strand') or '')
        text = str(value or '')
        parts = text.rsplit(':', 2)
        if len(parts) == 3 and parts[1].isdigit() and parts[2] in ('+', '-'):
            return cls(parts[0], int(parts[1]), parts[2])
        return cls(description=text)

    def label(self):
        return ('%s:%s (%s)' % (self.contig, self.position if self.position is not None else '?',
                              self.strand or '?')) if self.contig else self.description or 'Not recorded'


@dataclass(frozen=True)
class JunctionOrigin(DataclassSerializable):
    left: GenomicEndpoint
    right: GenomicEndpoint
    # Full source protein coordinates, zero-based half-open. A boundary can
    # be zero-width; a split codon spans one residue. Display, not admission.
    protein_intervals: tuple[tuple[int, int], ...] = ()
    description: str = ''

    def __post_init__(self):
        if not isinstance(self.left, GenomicEndpoint) or not isinstance(self.right, GenomicEndpoint):
            raise TypeError('Junction endpoints must be GenomicEndpoint records')
        intervals = tuple(tuple(interval) for interval in self.protein_intervals)
        if any(len(interval) != 2 or any(type(i) is not int for i in interval)
               or not 0 <= interval[0] <= interval[1] for interval in intervals):
            raise ValueError('Junction protein intervals must be nonnegative ordered coordinate pairs')
        object.__setattr__(self, 'protein_intervals', intervals)


@dataclass(frozen=True)
class AntigenOrigin(DataclassSerializable):
    variants: tuple[VariantOrigin, ...] = ()
    coding_references: tuple[CodingReference, ...] = ()
    rna_product: str = ''
    junctions: tuple[JunctionOrigin, ...] = ()

    def __post_init__(self):
        object.__setattr__(self, 'rna_product', self.rna_product.lower())
        if self.rna_product not in ('', 'fusion', 'splice', 'isoform'):
            raise ValueError('Unknown RNA product: ' + self.rna_product)
        for name, cls in (('variants', VariantOrigin), ('coding_references', CodingReference),
                          ('junctions', JunctionOrigin)):
            values = tuple(getattr(self, name))
            if any(not isinstance(value, cls) for value in values):
                raise TypeError(name + ' requires typed origin records')
            object.__setattr__(self, name, values)


def origin_from_isovar(result, family, candidate):
    variants = ()
    variant = result.get('variant')
    if family == 'protein_hypotheses':
        event = next((event for event in result['events'] if event['event_id'] == candidate['event_id']), {})
        variant = event.get('variant')
    if variant:
        variants = (VariantOrigin.from_alleles(variant),)
    elif family == 'sv':
        event = result.get('event_provenance') or {}
        variants = (VariantOrigin('SV', sv_type=str(result.get('event_type') or
            event.get('sv_type') or event.get('SVTYPE') or ''),
            reference_assembly=result.get('reference_name') or '', event_id=result.get('event_id') or ''),)
    junctions = []
    for junction in candidate.get('junctions', ()):
        if junction.get('annotated') or junction.get('termination_only'):
            continue
        start = candidate.get('coding_start')
        intervals = tuple(sorted({((boundary - start) // 3, (boundary - start + 2) // 3)
            for boundary in junction.get('query_interval', ()) if start is not None and boundary >= start}))
        junctions.append(JunctionOrigin(GenomicEndpoint.from_record(junction.get('left')),
            GenomicEndpoint.from_record(junction.get('right')), intervals))
    if not junctions and (result.get('donor') or result.get('acceptor')):
        junctions.append(JunctionOrigin(GenomicEndpoint.from_record(result.get('donor')),
            GenomicEndpoint.from_record(result.get('acceptor'))))
    reference = candidate.get('coding_reference')
    return AntigenOrigin(variants=variants,
        coding_references=(CodingReference(**reference),) if reference is not None else (),
        rna_product={'fusion': 'fusion', 'isoforms': 'isoform'}.get(family, ''), junctions=tuple(junctions))


def has_variant(antigen):
    """Historical labels remain interpretable without rewriting saved records."""
    return bool(antigen.origin and antigen.origin.variants) or antigen.kind in {'mutation', 'SV'}


def source_class(antigen):
    from .candidate_epitope import SOURCE_CLASS_MUTATION, SOURCE_CLASS_SELF, SOURCE_CLASS_VIRUS
    if antigen.kind == 'viral':
        return SOURCE_CLASS_VIRUS
    if has_variant(antigen) or antigen.kind in {'fusion', 'splice'}:
        return SOURCE_CLASS_MUTATION
    return SOURCE_CLASS_SELF


def origin_dict(antigen):
    return asdict(antigen.origin) if antigen.origin is not None else None


def origin_from_lens(kind, metadata):
    """Use caller-supplied loci; never infer translated boundaries from them."""
    metadata = dict(metadata)
    junctions = ()
    if kind == 'fusion':
        junctions = (JunctionOrigin(GenomicEndpoint.from_record(metadata.get('left_breakpoint')),
            GenomicEndpoint.from_record(metadata.get('right_breakpoint')),
            description=metadata.get('fusion_annotation', '')),)
    elif kind == 'splice':
        junctions = (JunctionOrigin(GenomicEndpoint(), GenomicEndpoint(),
            description='; '.join(metadata[k] for k in ('splice_coords', 'splice_description') if metadata.get(k))),)
    return AntigenOrigin(rna_product=kind if kind in {'fusion', 'splice'} else '', junctions=junctions)


def variant_labels(origin):
    """Display the same explicit event labels in tables, products and reports."""
    if origin is None:
        return ()
    events = origin.get('variants', ()) if isinstance(origin, dict) else map(asdict, origin.variants)
    return tuple(sorted({event['variant_type'] + (':' + event['sv_type'] if event.get('sv_type') else '')
                         for event in events}))
