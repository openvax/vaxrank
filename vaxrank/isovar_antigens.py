"""Lossless adaptation of released Isovar sequence/evidence contracts.

No reference sequence fills a gap in an observed candidate. All protein
alternatives, including held-out ones, stay in the antigen catalog; the full
producer result is stored once and linked by a content identifier.
"""

from dataclasses import asdict, dataclass, field, replace
from hashlib import sha256
import json

import msgspec
from isovar.genetic_code import genetic_code_by_id
from isovar.rna_evidence import content_identifier, evidence_set
from orfdata.reference import canonical_species

from .antigen_input import AntigenDataset, _typed_record
from .coding_context import reference_context
from .input_scope import InputScope, combine_scopes
from .isovar_options import source_options
from .vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation,
    TumorSpecificityEvidence, VaccineAntigen,
)


RESULT_FAMILIES = {
    'isovar.variant_rna_orfs.v1': 'variant_orfs',
    'isovar.sv_rna_candidates.v5': 'sv',
    'isovar.sv_rna_orfs.v4': 'sv',
    'isovar.fusion_rna.v3': 'fusion',
    'isovar.viral_sequences.v1': 'viral',
    'isovar.antigen_isoforms.v1': 'isoforms',
    'isovar.protein_hypotheses.v2': 'protein_hypotheses',
}


def _json(value):
    return json.dumps(value, sort_keys=True, allow_nan=False)


def _interval(pair, length):
    if (len(pair) != 2 or any(type(x) is not int for x in pair)
            or not 0 <= pair[0] <= pair[1] <= length):
        raise ValueError('Invalid Isovar interval')
    return tuple(pair)


def _translation(dna, protein, stop, table=1, start=False):
    if not isinstance(dna, str) or not dna or any(b not in 'ACGT' for b in dna):
        raise ValueError('Isovar coding sequence must be unambiguous DNA')
    if type(stop) is not bool or type(start) is not bool:
        raise ValueError('Isovar coding boundary flags must be booleans')
    amino, observed_stop = genetic_code_by_id(table).translate(dna)
    if start and amino:
        amino = 'M' + amino[1:]
    if (amino, observed_stop) != (protein, stop):
        raise ValueError('Isovar protein disagrees with its coding sequence/genetic code')


def _junction_mask(junctions, start, protein_length):
    # A split codon contains both partners. A zero-width boundary requires
    # flanking residues; stop-only crossings never become target residues.
    pairs = set()
    for junction in junctions:
        if junction.get('termination_only'):
            continue
        for boundary in junction['query_interval']:
            boundary -= start
            if 0 < boundary < 3 * protein_length:
                pairs.add((boundary // 3, (boundary + 2) // 3))
    # Adjacent split-codon masks can coincide/overlap around short inserts.
    merged = []
    for low, high in sorted(pairs):
        if merged and low < merged[-1][1]:
            merged[-1] = (merged[-1][0], max(high, merged[-1][1]))
        else:
            merged.append((low, high))
    return tuple(merged)


@dataclass
class _Candidate:
    identifier: str
    kind: str
    protein: str
    intervals: tuple
    fragments: int | None
    partial: bool = False
    ambiguous: bool = False
    reasons: tuple = ()
    gene_names: tuple = ()
    gene_ids: tuple = ()
    transcripts: tuple = ()
    proteins: tuple = ()
    species: str = ''
    evidence: dict = field(default_factory=dict)


def _variant_frame(row, reference):
    """Compare frames at the focal allele, independent of window cropping."""
    length = len(row['nucleotide_sequence'])
    coding_start, _ = _interval(row['coding_interval'], length)
    variant_start, _ = _interval(row['variant_cdna_interval'], length)
    model = reference['model']
    resolved = reference['resolved_transcript']
    strand = resolved['strand'] if 'strand' in resolved else (model['strand'] if model is not None else None)
    if strand is not None and strand not in ('+', '-'):
        raise ValueError('Isovar reference strand must be + or -')
    return strand, (variant_start - coding_start) % 3, reference['translation_table']


def _variants(result):
    models = {m['model_id']: m['reference'] for m in result['models']}
    frames = {}
    for row in result['protein_candidates']:
        reference = models[row['model_id']]
        strand, phase, table = _variant_frame(row, reference)
        # Unknown carrier orientation (including historical exports without
        # strand) only establishes a shared frame within the same transcript.
        orientation = strand or (reference['reference_key'], reference['resolved_transcript']['id'])
        frames.setdefault(row['rna_candidate_id'], set()).add(
            (orientation, phase, table))
    for row in result['protein_candidates']:
        reference = models[row['model_id']]
        strand, phase, table = _variant_frame(row, reference)
        amino = row['protein_sequence']
        low, high = _interval(row['coding_interval'], len(row['nucleotide_sequence']))
        _translation(row['nucleotide_sequence'][low:high], amino, row['ends_with_stop_codon'],
                     reference['translation_table'], row['starts_at_reference_start'])
        interval = _interval(row['mutation_interval'], len(amino))
        changed = row['focal_amino_acid_change']
        if type(changed) is not bool:
            raise ValueError('Isovar focal_amino_acid_change must be boolean')
        intervals = (interval,) if changed and interval != (len(amino), len(amino)) else ()
        support = row['rna_support']
        yield _Candidate(row['candidate_id'], 'mutation', amino, intervals, support['fragments'],
            partial=row['sequence_extent'] != 'complete_coding_candidate',
            ambiguous=len(frames[row['rna_candidate_id']]) > 1,
            reasons=() if row['coding_support']['fully_tiled'] else ('coding_interval_not_observed',),
            gene_names=(row['gene_name'],),
            gene_ids=(reference.get('source_gene_id'),) if reference['on_source_transcript'] else (),
            transcripts=(row['transcript_id'],), proteins=(row['orf_id'],),
            species=reference['resolved_transcript'].get('species') or '',
            evidence=dict(model_id=row['model_id'], rna_candidate_id=row['rna_candidate_id'],
                          coding_reference=reference_context(reference),
                          rna_support=support, coding_support=row['coding_support'],
                          reference_orf_source_gene_id=reference.get('source_gene_id'),
                          reference_strand=strand, variant_codon_phase=phase,
                          sequence_extent=row['sequence_extent'], translation_table=table))


def _sv_support(support, result, scope):
    """Validate the source-qualified set carried by a measured SV placement."""
    evidence_scope = [result.get('sample_id'), result.get('source')]
    if (any(not isinstance(value, str) or not value.strip() for value in evidence_scope)
            or not isinstance(support, dict) or support.get('scope') != scope
            or support.get('evidence_scope') != evidence_scope):
        raise ValueError('Isovar SV RNA support has an inconsistent scope')
    for unit in ('read', 'fragment'):
        ids = support.get(unit + '_ids')
        if (not isinstance(ids, list) or any(not isinstance(key, str) or not key for key in ids)
                or len(set(ids)) != len(ids) or type(support.get(unit + 's')) is not int
                or support[unit + 's'] != len(ids)):
            raise ValueError('Isovar SV RNA support counts disagree with its identities')
    if (support.get('evidence_set_id') != content_identifier('evidence', sorted(support['read_ids']))
            or support['fragments'] > support['reads'] or bool(support['fragments']) != bool(support['reads'])):
        raise ValueError('Isovar SV RNA evidence identity is inconsistent')
    return support


def _sv_export(result, observations=None):
    for row in result['candidates']:
        _translation(row['nucleotide_sequence'], row['amino_acids'], row['ends_with_stop_codon'], start=True)
        measured = [o for o in row['occurrences'] if 'rna_support' in o]
        if measured:
            union = _sv_support(row['rna_support'], result, 'union_of_full_interval_junction_linked_witnesses_within_input')
            for occurrence in measured:
                scoped = _sv_support(occurrence['rna_support'], result,
                                     'full_interval_junction_linked_witnesses_for_occurrence')
                if bool(occurrence['witnesses']) != bool(scoped['reads']):
                    raise ValueError('Isovar SV witnesses disagree with measured RNA support')
                for ids in ('read_ids', 'fragment_ids'):
                    if not set(scoped[ids]) <= set(union[ids]):
                        raise ValueError('Isovar SV occurrence support is outside the candidate union')
                if observations is not None:
                    expected = evidence_set(scoped['evidence_scope'],
                        {tuple(observations[w['observation']]['identity']) for w in occurrence['witnesses']})
                    if any(set(scoped[ids]) != set(expected[ids]) for ids in ('read_ids', 'fragment_ids')):
                        raise ValueError('Isovar SV occurrence support disagrees with original observations')
            if len(measured) == len(row['occurrences']):
                for ids in ('read_ids', 'fragment_ids'):
                    if set(union[ids]) != {key for o in measured for key in o['rna_support'][ids]}:
                        raise ValueError('Isovar SV candidate union disagrees with measured occurrences')
        for index, occurrence in enumerate(row['occurrences']):
            low, high = occurrence['query_interval']
            if high - low != len(row['nucleotide_sequence']):
                raise ValueError('Isovar SV occurrence length disagrees with its ORF')
            transcripts = tuple(c['transcript_id'] for c in occurrence['reference_comparisons']
                                if c.get('transcript_id'))
            support = None
            if 'rna_support' in occurrence:
                support = occurrence['rna_support']
            elif observations is not None:
                from isovar.rna_evidence import rna_support
                support = rna_support({tuple(observations[w['observation']]['identity'])
                                       for w in occurrence['witnesses']})
            elif not occurrence['witnesses']:
                support = {'fragments': 0}
            elif len(row['occurrences']) == 1:
                support = row['rna_support']
            yield _Candidate(row['candidate_id'] + ':occurrence:' + str(index), 'SV', row['amino_acids'],
                _junction_mask(occurrence['junctions'], low, len(row['amino_acids'])),
                support['fragments'] if support is not None else None, partial=not row['ends_with_stop_codon'],
                ambiguous=occurrence['frame_status'] != 'translated' or not occurrence['annotated_start'],
                transcripts=transcripts,
                evidence=dict(candidate_id=row['candidate_id'], occurrence_index=index,
                              path_id=occurrence['path_id'], rna_support=row['rna_support'],
                              occurrence_rna_support=support,
                              uncertainty_flags=row['uncertainty_flags'], translation_table=1))


def _sv(result):
    from isovar import export_sv_rna_orfs
    # Keep annotated-frame products too: exploratory ATG exports alone omit
    # local protein windows whose initiation site lies outside observed RNA.
    for path in result['paths']:
        for index, row in enumerate(path['translations']):
            low, high = _interval([row['translation_start'], row['translation_end']], len(path['sequence']))
            _translation(path['sequence'][low:high], row['amino_acids'], row['ends_with_stop_codon'],
                         start=row['complete_5prime'])
            departures = {e['departure'] for e in row['frame_evidence']}
            junctions = [j for j in path['junctions'] if not j['annotated']
                         and departures.intersection(j['query_interval'])]
            supports = [j['direct_support']['fragments'] for j in junctions]
            yield _Candidate(path['path_id'] + ':translation:' + str(index), 'SV', row['amino_acids'],
                _junction_mask(junctions, low, len(row['amino_acids'])), min(supports) if supports else None,
                partial=not (row['complete_5prime'] and row['ends_with_stop_codon']),
                ambiguous=path['frame_status'] != 'translated', transcripts=tuple(row['transcript_ids']),
                evidence=dict(path_id=path['path_id'], translation_index=index, translation_table=1,
                              frame_evidence=row['frame_evidence'], support_scope='minimum_direct_junction_fragments'))
    yield from _sv_export(export_sv_rna_orfs(result), observations=result['observations'])


def _coding(result, family):
    viral = family == 'viral'
    sequences = result['sequences'] if viral else result['isoforms'] + result['partial_sequences']
    targets = {} if viral else {t['reference']['transcript_id']: t for t in result['targets']}
    for sequence in sequences:
        for index, row in enumerate(sequence['coding']):
            low, high = _interval([row['query_start'], row['query_end']], len(sequence['sequence']))
            dna = row['coding_sequence'] if viral else sequence['sequence'][low:high]
            _translation(dna, row['amino_acids'], row['observed_stop'], row.get('translation_table', 1), row['complete_start'])
            references = row['reference_orfs'] if viral else []
            transcripts = tuple(row.get('transcript_ids', ()))
            names = ([r.get('gene_name') for r in references] if viral else
                     [targets[t]['gene_name'] for t in transcripts if t in targets])
            taxons = {r['reference_context'].get('species') for r in references
                      if r.get('reference_context', {}).get('species')}
            if len(taxons) > 1:
                raise ValueError('Isovar coding prediction has conflicting source species')
            yield _Candidate(row.get('id', sequence['id'] + ':coding:' + str(index)),
                'viral' if viral else 'isoform', row['amino_acids'], ((0, len(row['amino_acids'])),),
                sequence['rna_support']['fragments'], partial=row['status'] != 'complete_CDS',
                reasons=() if row['all_bases_observed'] else ('coding_interval_not_observed',),
                gene_names=tuple(names), gene_ids=tuple(sequence.get('gene_ids', ())),
                transcripts=transcripts, proteins=tuple(r.get('protein_id') for r in references),
                species=next(iter(taxons), ''),
                evidence=dict(sequence_id=sequence['id'], coding_index=index,
                    sequence_status=sequence['sequence_status'], origin=sequence['origin'],
                    rna_support=sequence['rna_support'], all_bases_observed=row['all_bases_observed'],
                    whole_cds_segments=row['whole_cds_segments'], translation_table=row.get('translation_table', 1),
                    reference_orfs=references))


def _proteins(result):
    for event in result['events']:
        for row in event['protein_hypotheses']:
            for translation in row['translations']:
                low, high = _interval(translation['translated_interval'], len(translation['nucleotide_sequence']))
                # Isovar's legacy protein export translates the first codon as
                # ordinary sequence, even at an annotated non-ATG start.
                _translation(translation['nucleotide_sequence'][low:high], row['amino_acids'],
                             row['ends_with_stop_codon'], translation['reference_context']['translation_table'])
            interval = _interval(row['mutation_interval'], len(row['amino_acids']))
            intervals = (interval,) if row['contains_mutation'] and interval != (len(row['amino_acids']),) * 2 else ()
            yield _Candidate(row['hypothesis_id'], 'mutation', row['amino_acids'], intervals,
                row['rna_support']['fragments'],
                partial=row['n_terminus'] != 'annotated_start_codon' or not row['ends_with_stop_codon'],
                reasons=() if event['filters']['passes_all_filters'] else ('isovar_filters_failed',),
                gene_names=tuple(row['gene_names']), gene_ids=tuple(row['gene_ids']),
                transcripts=tuple(row['transcript_ids']),
                evidence=dict(event_id=event['event_id'], rna_support=row['rna_support'],
                              isovar_rank=row['isovar_rank'], filters=event['filters']))


def _fusion(result):
    from .isovar_fusion import fusion_antigens_from_isovar
    adapted = fusion_antigens_from_isovar(result)
    path, = result['paths']
    for index, (antigen, row) in enumerate(zip(adapted.antigens, path['translations'])):
        yield _Candidate(antigen.source_identifier, 'fusion', antigen.amino_acids,
            tuple((i.start, i.end) for i in antigen.targetable_mask.intervals),
            path['junctions'][0]['direct_support']['fragments'],
            partial=not (row['complete_5prime'] and row['ends_with_stop_codon']),
            ambiguous=result['status'] != 'translated', transcripts=antigen.transcript_ids,
            evidence=dict(translation_index=index, translation_table=1))


def _attestation(value):
    if value is None:
        return TumorSpecificityAttestation(status='held_out', evidence_kind='missing_tumor_specificity',
            evidence_source='Isovar RNA reconstruction', rationale_code='rna_is_not_tumor_specificity', requires_review=True)
    if isinstance(value, TumorSpecificityAttestation):
        return value
    value = dict(value)
    value['evidence_records'] = tuple(_typed_record(TumorSpecificityEvidence, r)
                                     for r in value.get('evidence_records', ()))
    return _typed_record(TumorSpecificityAttestation, value)


def antigens_from_isovar(result, *, scope=None, species='', tumor_specificity=None,
                         candidate_attestations=None, options=None):
    """Adapt all supported Isovar hypotheses into one native antigen dataset.

    ``tumor_specificity`` is an explicit source-wide assertion; use
    ``candidate_attestations`` to restrict it to reviewed candidate IDs instead.
    RNA support never supplies this assertion. Partial/ambiguous candidates
    require the corresponding admission opt-in in addition to an attestation.
    Full source results and decisions survive native antigen/prediction replay.
    """
    result = json.loads(_json(result))
    family = RESULT_FAMILIES.get(result.get('schema'))
    if family is None:
        raise ValueError('Unsupported Isovar result schema: ' + str(result.get('schema')))
    options = source_options(options)
    if options.reads or options.protein or options.reconstruction:
        raise ValueError('Reconstruction controls cannot be applied to a saved Isovar result')
    sample = result.get('sample_id') or result.get('provenance', {}).get('sample_id')
    if not sample:
        raise ValueError('Isovar result requires a sample_id')
    scope = combine_scopes(scope or InputScope(), InputScope(sample_id=sample), 'Isovar result')
    reference = result.get('reference_name') or result.get('variant', {}).get('reference_name')
    if reference and scope.reference_assembly and reference != scope.reference_assembly:
        raise ValueError('Isovar result and host scope reference assemblies disagree')
    identity = 'isovar:' + sha256(_json(result).encode()).hexdigest()
    candidates = (_variants(result) if family == 'variant_orfs' else _sv(result)
                  if result['schema'] == 'isovar.sv_rna_candidates.v5' else _sv_export(result)
                  if family == 'sv' else _fusion(result) if family == 'fusion' else _proteins(result)
                  if family == 'protein_hypotheses' else _coding(result, family))
    antigens, decisions = {}, []
    candidate_attestations = dict(candidate_attestations or {})
    default_attestation = _attestation(tumor_specificity)
    policy = options.admission
    for candidate in candidates:
        key = candidate.identifier
        if not isinstance(key, str) or not key or key in antigens:
            raise ValueError('Isovar candidate identifiers must be nonempty and unique')
        requested = _attestation(candidate_attestations.get(key, default_attestation))
        reasons = list(candidate.reasons)
        if not candidate.intervals:
            reasons.append('no_targetable_amino_acid_change')
        if candidate.fragments is None:
            reasons.append('rna_fragment_support_unavailable')
        elif type(candidate.fragments) is not int or candidate.fragments < 0:
            raise ValueError('Isovar fragment support must be a nonnegative integer')
        elif candidate.fragments < policy.min_fragments:
            reasons.append('insufficient_rna_fragments')
        if candidate.partial and not policy.allow_partial:
            reasons.append('partial_coding_candidate')
        if candidate.ambiguous and not policy.allow_ambiguous:
            reasons.append('ambiguous_coding_candidate')
        attestation = requested
        if reasons:
            attestation = replace(requested, status='held_out', rationale_code='isovar:' + ','.join(reasons),
                                  requires_review=True, override_reason='')
        names = tuple(sorted({s for s in candidate.gene_names if s}))
        genes = tuple(sorted({s for s in candidate.gene_ids if s}))
        if species and candidate.species and canonical_species(species) != canonical_species(candidate.species):
            raise ValueError('Declared species disagrees with Isovar source species')
        antigens[key] = VaccineAntigen(kind=candidate.kind, amino_acids=candidate.protein,
            targetable_mask=TargetableMask(tuple(AminoAcidInterval(*i) for i in candidate.intervals)),
            tumor_specificity=attestation, gene_name=';'.join(names), gene_id=genes[0] if len(genes) == 1 else '',
            transcript_ids=tuple(sorted(set(candidate.transcripts))),
            protein_ids=tuple(sorted({s for s in candidate.proteins if s})), species=candidate.species or species,
            source_identifier=key, source_metadata=(('isovar_evidence_id', identity), ('isovar_family', family),
                ('gene_ids', _json(genes)), ('isovar_candidate', _json(dict(
                    candidate.evidence, partial_coding_candidate=candidate.partial,
                    ambiguous_coding_candidate=candidate.ambiguous)))))
        decisions.append(dict(candidate_id=key, status=attestation.status, reasons=reasons,
                              requested_tumor_specificity=asdict(requested)))
    if set(candidate_attestations) - set(antigens):
        raise ValueError('Attestations name unknown Isovar candidates')
    evidence = dict(isovar_result=result, family=family, options=msgspec.to_builtins(options), decisions=decisions)
    return AntigenDataset(antigens, scope=scope, source_evidence={identity: evidence})


def isovar_evidence_summary(source_evidence):
    """Human-readable source outcomes, including inputs with no proteins."""
    from collections import Counter
    rows = []
    for evidence in source_evidence.values():
        if 'isovar_result' not in evidence:
            continue
        result, decisions = evidence['isovar_result'], evidence['decisions']
        admitted = sum(d['status'] in ('admitted', 'overridden') for d in decisions)
        reasons = Counter(reason for d in decisions for reason in d['reasons'])
        reasons.update(m['reason'] for m in result.get('models', ()) if m.get('reason'))
        for decision in decisions:
            if decision['status'] == 'held_out' and not decision['reasons']:
                reasons[decision['requested_tumor_specificity'].get('rationale_code') or 'source_attestation_held_out'] += 1
        state = result.get('assembly_state') or result.get('status')
        text = (f"{evidence['family']} / {result.get('sample_id') or result.get('provenance', {}).get('sample_id')}: "
                f'{len(decisions)} coding candidates, {admitted} admitted, {len(decisions) - admitted} held out')
        if state:
            text += '; producer outcome: ' + state
        if reasons:
            text += '; ' + ', '.join(f'{key} ({count})' for key, count in sorted(reasons.items()))
        rows.append(text)
    return rows
