"""Select distinct emitted sequences while retaining their separate sources."""

from dataclasses import dataclass, field
import json
import logging

from varcode import Variant

from .vaccine_library import iter_named_antigens, select_antigen_window, top_target_epitopes
from .coding_context import coding_contexts

logger = logging.getLogger(__name__)


@dataclass
class AssemblyAntigen:
    """One sequence to pack, with all contributing candidate mappings.

    Sources remain separate observations: their scores and RNA support are not
    added. The first occurrence in selection order supplies the display name.
    """
    name: str
    amino_acids: str
    sources: list = field(default_factory=list)


def _source_record(source, name, antigen, peptide, source_evidence):
    metadata = dict(antigen.source_metadata)
    evidence_id = metadata.get('isovar_evidence_id')
    result = source_evidence.get(evidence_id, {}).get('isovar_result', {})
    if not result and metadata.get('isovar_fusion_result'):
        result = json.loads(metadata['isovar_fusion_result'])
    event = None
    if isinstance(source, Variant):
        event = dict(reference_assembly=source.ensembl.reference_name,
                     contig=source.contig, start=source.start, ref=source.ref, alt=source.alt)
    elif result.get('variant'):
        event = result['variant']
    elif result.get('event_id'):
        event = dict(event_id=result['event_id'], reference_assembly=result.get('reference_name'))
    genes = set(json.loads(metadata.get('gene_ids', '[]')))
    if antigen.gene_id:
        genes.add(antigen.gene_id)
    fragment = peptide.mutant_protein_fragment
    genes.update(t.gene_id for t in getattr(fragment, 'supporting_reference_transcripts', ())
                 if getattr(t, 'gene_id', None))
    record = dict(
        candidate_name=name, candidate_sequence=peptide.amino_acids,
        source_identifier=antigen.source_identifier, kind=antigen.kind,
        gene_name=antigen.gene_name, gene_ids=sorted(genes),
        transcript_ids=list(antigen.transcript_ids), protein_ids=list(antigen.protein_ids),
        species=antigen.species,
        prediction_ids=sorted({e.prediction_group_source for e in peptide.epitopes}),
        input_ids=sorted({e.input_provenance.input_id for e in peptide.epitopes if e.input_provenance}),
        evidence_group=evidence_id, source_event=event)
    contexts = coding_contexts(antigen, source_evidence)
    if contexts or ('isovar_family', 'variant_orfs') in antigen.source_metadata:
        record['coding_contexts'] = contexts
    return record


def antigen_occurrences(ranked, options, modality, *, source_evidence=None):
    """Yield exact windows/ligands using the same content rules in both modalities."""
    if options.antigen_content not in ('mutation_spanning', 'minimal_epitope'):
        raise ValueError("antigen_content must be 'mutation_spanning' or 'minimal_epitope'")
    sources = {id(peptide): source for source, peptides in ranked for peptide in peptides}
    for name, antigen, peptide in iter_named_antigens(ranked, options.candidates_per_slot):
        source = _source_record(sources[id(peptide)], name, antigen, peptide, source_evidence or {})
        if options.antigen_content == 'minimal_epitope':
            targets = top_target_epitopes(peptide, n=options.epitopes_per_antigen)
            if not targets:
                logger.info('Skipping %s: no target epitope predictions available.', name)
            for index, epitope in enumerate(targets):
                suffix = 'epitope' if len(targets) == 1 else 'epitope%d' % (index + 1)
                # Preserve the established modality labels; the source mapping
                # records the actual emitted name rather than guessing a suffix.
                label = (name[:-1] + ', ' + suffix + ')' if modality == 'peptide' and name.endswith(')')
                         else name + ' (' + suffix + ')' if modality == 'peptide' else name + '_' + suffix)
                yield AssemblyAntigen(label, epitope.sequence, [source])
        else:
            sequence = select_antigen_window(antigen, name, options.max_antigen_length_aa)
            if len(sequence) < options.min_antigen_length_aa:
                logger.warning('Antigen %s emitted at %d aa, below --%s-min-antigen-length-aa (%d).',
                               name, len(sequence), modality, options.min_antigen_length_aa)
            yield AssemblyAntigen(name, sequence, [source])


def antigens_for_assembly(ranked, options, modality, *, source_evidence=None):
    """Group identical final antigen sequences before any packing/product cap.

    Distinct overlapping windows remain distinct sequences. Shared event and
    evidence-group identities make their related provenance explicit; neither
    sequence count nor candidate count establishes independent target coverage.
    """
    by_sequence, used_names = {}, set()
    for occurrence in antigen_occurrences(ranked, options, modality, source_evidence=source_evidence):
        existing = by_sequence.get(occurrence.amino_acids)
        if existing is not None:
            existing.sources.extend(s for s in occurrence.sources if s not in existing.sources)
            continue
        name, index = occurrence.name, 2
        while occurrence.name in used_names:
            occurrence.name = name + ' [sequence %d]' % index
            index += 1
        used_names.add(occurrence.name)
        by_sequence[occurrence.amino_acids] = occurrence
    return list(by_sequence.values())


def summarize_products(products):
    """Count actual products, sequences and recorded events without inferring independence."""
    antigens = [a for product in products for a in product.antigens]
    sources = {json.dumps(s, sort_keys=True) for a in antigens for s in a.get('sources', ())}
    events, unassigned, evidence_groups = set(), 0, set()
    for serialized in sources:
        source = json.loads(serialized)
        if source['source_event'] is None:
            unassigned += 1
        else:
            events.add(json.dumps(source['source_event'], sort_keys=True))
        if source['evidence_group']:
            evidence_groups.add(source['evidence_group'])
    return dict(constructs=len(products), unique_construct_sequences=len({p.sequence for p in products}),
                antigen_slots=sum(len(p.antigen_names) for p in products),
                unique_antigen_sequences=(len({a['aa'] for a in antigens})
                    if len(antigens) == sum(len(p.antigen_names) for p in products) else None),
                source_candidates=len(sources), source_events=len(events),
                source_candidates_without_event=unassigned, evidence_groups=len(evidence_groups))


def source_matches_candidate(source, antigen, sequence, prediction_ids, input_ids):
    """Match one retained product source to a candidate's exact evidence identity."""
    return (source['source_identifier'] == antigen.source_identifier
            and source['gene_name'] == antigen.gene_name
            and set(source['input_ids']) <= input_ids
            and set(source['prediction_ids']) <= prediction_ids
            and source['candidate_sequence'] in sequence)


def emitted_antigens_for_candidate(antigen, sequence, epitopes, products):
    """Map a ranked candidate to actual product/antigen names using retained evidence."""
    prediction_ids = {e.prediction_group_source for e in epitopes}
    input_ids = {e.input_provenance.input_id for e in epitopes if e.input_provenance}
    return {product.name: sorted({name for name, sources in product.antigen_sources.items()
        for source in sources
        if source_matches_candidate(source, antigen, sequence, prediction_ids, input_ids)}) for product in products}
