"""Exact peptide containment against OncoRef's verified sequence roles."""

from dataclasses import replace
from collections import defaultdict
import hashlib
import math
from pathlib import Path

from .vaccine_antigen import SelfReferenceMatch, SelfReferenceSource


SCHEMA = 'vaxrank.cta_peptide_targeting.v1'


def load_bundle(path):
    from oncoref import sequence_bundle
    return sequence_bundle.load(None if path == 'current' else path)


def bundle_scope(bundle):
    from .cta_reference_adapter import _reference_context
    from .input_scope import InputScope
    from orfdata.reference import species_name_for_taxon
    reference = bundle.manifest['reference']
    context = _reference_context(reference, species_name_for_taxon(reference['taxon_id']))
    return InputScope(host_taxon_id=context.taxon_id, reference_assembly=context.genome_build,
                      annotation='%s:%s' % (context.annotation_name, context.annotation_version))


def admit_expression(path, contract, policy, bundle, scope):
    """Resolve current qualified representatives directly from the owner bundle."""
    from orfdata.reference import species_name_for_taxon
    from .cta_expression import (CTAExpressionDecision, CTAExpressionResult, CTASequenceAssessment,
                                 CTATargetExclusionPolicy, _expression_rows)
    from .input_scope import combine_scopes
    from .vaccine_antigen import (AminoAcidInterval, TargetableMask, TumorSpecificityAttestation,
                                 TumorSpecificityEvidence, VaccineAntigen)
    reference = bundle_scope(bundle)
    combine_scopes(scope, reference, 'CTA sequence bundle')
    if contract.expression_unit != policy.expression_unit:
        raise ValueError('Expression input units do not match the CTA admission policy')
    annotation_name, annotation_version = reference.annotation.split(':', 1)
    expression = list(_expression_rows(path, contract))
    wanted = {feature for _, feature, _ in expression}
    by_feature = {}
    for record in bundle.records():
        for occurrence in record['occurrences']:
            key = occurrence.get('gene_id' if contract.measurement_level == 'gene' else 'transcript_id')
            if key and key.split('.')[0] in wanted:
                by_feature.setdefault(key.split('.')[0], []).append((record, occurrence))
    decisions = []
    for identifier, feature, value in expression:
        measured = value is not None and math.isfinite(value)
        candidates = by_feature.get(feature, [])
        if contract.measurement_level == 'transcript' and '.' in identifier:
            candidates = [(r, o) for r, o in candidates if o['transcript_id'] == identifier]
        allowed = [(r, o) for r, o in candidates if r['direct_target_candidate']
                   and o['direct_target_candidate'] and r['shared_sequence'] == 'allow']
        if not allowed:
            decisions.append(CTAExpressionDecision(identifier, float(value) if measured else None,
                status='unsafe_reference' if any(r['shared_sequence'] == 'exclude' for r, _ in candidates)
                       else 'not_current_cta' if candidates else 'unresolved_reference'))
            continue
        groups = defaultdict(list)
        for record, occurrence in allowed:
            groups[(occurrence['gene_id'], record['sequence'])].append((record, occurrence))
        for (gene_id, sequence), origins in groups.items():
            record, occurrence = origins[0]
            transcripts = tuple(sorted({o['transcript_id'] for _, o in origins}))
            proteins = tuple(sorted({r['record_id'] for r, _ in origins}))
            sequence = record['sequence']
            # Qualified full representatives must be contiguous AA20 records.
            complete = (len(record['valid_segments']) == 1 and record['valid_segments'][0]['start'] == 0
                        and record['valid_segments'][0]['end'] == len(sequence))
            if not complete:
                raise ValueError('OncoRef allowed a representative with interrupted sequence')
            admitted = measured and value >= policy.min_tumor_expression
            antigen = VaccineAntigen(kind='CTA', amino_acids=sequence,
                targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
                gene_id=gene_id, gene_name=occurrence.get('gene_name') or '',
                transcript_ids=transcripts, protein_ids=proteins,
                species=species_name_for_taxon(reference.host_taxon_id),
                source_identifier=gene_id + ':' + (
                    identifier + ':' if contract.measurement_level == 'transcript' else ''
                ) + hashlib.sha256(sequence.encode()).hexdigest(),
                source_metadata=(('oncoref_sequence_bundle', bundle.identity),),
                source_protein_completeness='complete',
                tumor_specificity=TumorSpecificityAttestation(
                    status='admitted' if admitted else 'held_out', evidence_kind='cta_expression',
                    evidence_source='OncoRef current representative and patient expression',
                    patient_specific=True, rationale_code='expressed_current_cta' if admitted else 'insufficient_tumor_expression',
                    evidence_records=(TumorSpecificityEvidence(evidence_kind='patient_expression',
                        evidence_source=contract.evidence_source, evidence_version=contract.evidence_version,
                        subject_id=identifier, sample_id=contract.sample_id, patient_specific=True,
                        numeric_value=float(value) if measured else None, unit=contract.expression_unit,
                        threshold=policy.min_tumor_expression, passed=admitted),)))
            assessment = CTASequenceAssessment(antigen, bundle.identity, record['record_id'], occurrence,
                equivalent_occurrences=tuple(dict(record_id=r['record_id'], **o) for r, o in origins))
            decisions.append(CTAExpressionDecision(identifier, float(value) if measured else None,
                gene_id=occurrence['gene_id'], gene_name=antigen.gene_name,
                transcript_id=';'.join(transcripts), status='admitted' if admitted else 'held_out',
                assessment=assessment, gene_identity=dict(canonical_gene_id=occurrence['gene_id'],
                    source_gene_id=occurrence['gene_id'], status='verified_sequence_bundle')))
    return CTAExpressionResult(contract, policy, CTATargetExclusionPolicy((), ()),
        hashlib.sha256(Path(path).read_bytes()).hexdigest(), reference.reference_assembly,
        annotation_name, annotation_version, tuple(decisions), target_reference=dict(
            bundle_id=bundle.identity, manifest=bundle.manifest, coverage=bundle.coverage))


def peptide_containment(peptides, bundle):
    """Match all contiguous occurrences, preserving exclusion/review precedence."""
    import ahocorasick

    peptides = sorted(set(peptides))
    matches = {peptide: [] for peptide in peptides}
    if not peptides:
        return matches
    query = ahocorasick.Automaton()
    for peptide in peptides:
        query.add_word(peptide, peptide)
    query.make_automaton()
    for segment in bundle.segments(min_length=min(map(len, peptides))):
        for end, peptide in query.iter(segment['sequence']):
            matches[peptide].append(dict(
                record_id=segment['record_id'], role=segment['shared_sequence'],
                start=segment['start'] + end + 1 - len(peptide), end=segment['start'] + end + 1,
                occurrences=segment['occurrences']))
    return matches


def targeting_evidence(epitope):
    return next((record for record in epitope.input_evidence if record.get('schema') == SCHEMA), None)


def recorded_bundle_ids(dataset):
    """Identify frozen CTA policies from evidence, without resolving live aliases."""
    identities = set()
    for epitope in dataset.epitopes:
        antigen = dataset.antigens.get(epitope.prediction_group_source)
        if antigen is not None and antigen.kind == 'CTA':
            evidence = targeting_evidence(epitope)
            identities.add(evidence['bundle_id'] if evidence else None)
    for antigen in dataset.antigens.values():
        if antigen.kind == 'CTA':
            identities.update(value for key, value in antigen.source_metadata
                              if key == 'oncoref_sequence_bundle')
    # Admission-only sources can have no retained predictions.
    for identity, evidence in dataset.selection.get('source_evidence', {}).items():
        if evidence.get('schema') == 'vaxrank.oncoref_target_reference.v1':
            identities.add(identity)
    return identities


def matching_gene_summary(matches, roles=('exclude', 'review')):
    """Compact display only; exact protein occurrences stay in native evidence."""
    genes = defaultdict(set)
    for match in matches:
        if match['role'] in roles:
            for source in match['occurrences']:
                if source.get('shared_sequence', match['role']) in roles:
                    genes[source.get('gene_name') or source.get('gene_id') or '?'].add(match['record_id'])
    return '; '.join('%s (%d protein%s)' % (gene, len(records), '' if len(records) == 1 else 's')
                     for gene, records in sorted(genes.items()))


def _targeting_summary(epitopes, antigens, bundle_id):
    """Compact ranges of tested peptide starts, including unranked proteins."""
    sources = {}
    for epitope in epitopes:
        antigen = antigens[epitope.prediction_id]
        provenance = epitope.input_provenance.input_id if epitope.input_provenance else ''
        key = (provenance, antigen.gene_id or antigen.gene_name, antigen.amino_acids)
        summary, groups = sources.setdefault(key, (dict(
            schema='vaxrank.cta_targeting_summary.v1', gene=antigen.display_gene_name,
            gene_id=antigen.gene_id, identifier=antigen.source_identifier, identifiers=[],
            bundle_id=bundle_id), {}))
        if antigen.source_identifier not in summary['identifiers']:
            summary['identifiers'].append(antigen.source_identifier)
        evidence = targeting_evidence(epitope)
        starts, matches = groups.setdefault((len(epitope.sequence), evidence['status']), (set(), {}))
        starts.add(antigen.source_protein_offset + epitope.offset + 1)
        for match in evidence['matches']:
            matches[match['record_id']] = match
    for key, (summary, groups) in sources.items():
        rows = []
        for (length, status), (starts, matches) in sorted(groups.items()):
            ranges = []
            for position in sorted(starts):
                if ranges and position == ranges[-1][1] + 1:
                    ranges[-1][1] = position
                else:
                    ranges.append([position, position])
            rows.append(dict(length=length, status=status, count=len(starts),
                starts=', '.join(str(a) if a == b else '%d–%d' % (a, b) for a, b in ranges),
                matches=matching_gene_summary(matches.values(),
                    roles=('exclude',) if status == 'excluded' else ('review',))))
        summary['rows'] = rows
        identity = 'cta-targeting:' + hashlib.sha256(repr((bundle_id, key)).encode()).hexdigest()
        yield identity, summary


def apply_cta_targeting(dataset, bundle):
    """Attach shared targeting decisions before ranking; native replay is offline.

    Source admission remains separate. A peptide can be targeted only if it has
    a qualified allowed occurrence and no exclusion or unresolved occurrence.
    This is exact containment, not a near-self or clinical safety model.
    """
    cta = [epitope for epitope in dataset.epitopes
           if dataset.antigens[epitope.prediction_id].kind == 'CTA']
    if not cta:
        return
    manifest = bundle.manifest
    reference = manifest['reference']
    from .input_scope import combine_scopes
    reference_scope = bundle_scope(bundle)
    for provenance in dataset.provenance:
        combine_scopes(provenance.scope, reference_scope, 'OncoRef sequence bundle')
    matches = peptide_containment((e.sequence for e in cta), bundle)
    updated = {}
    for epitope in cta:
        occurrences = matches[epitope.sequence]
        roles = {m['role'] for m in occurrences}
        status = ('excluded' if 'exclude' in roles else 'review' if 'review' in roles
                  else 'eligible' if 'allow' in roles else 'unassessed')
        evidence = dict(schema=SCHEMA, peptide=epitope.sequence, status=status,
                        bundle_id=bundle.identity, matches=occurrences)
        antigen = dataset.antigens[epitope.prediction_id]
        sources = []
        for match in occurrences:
            if match['role'] != 'exclude':
                continue
            for index, source in enumerate(match['occurrences']):
                if source.get('gene_id') and source.get('shared_sequence', 'exclude') == 'exclude':
                    sources.append(SelfReferenceSource(gene_id=source['gene_id'],
                        gene_name=source.get('gene_name') or '', transcript_id=source.get('transcript_id') or '',
                        protein_id=source.get('protein_id') or match['record_id'],
                        reference_entry=match['record_id'], reference_occurrence=index))
        updated[epitope.prediction_id] = replace(epitope,
            occurs_in_reference=bool(occurrences),
            occurs_in_non_CTA_reference=status == 'excluded',
            overlaps_targetable=(antigen.interval_is_targetable(
                epitope.offset, epitope.offset + len(epitope.sequence)) and status == 'eligible'),
            self_reference_match=SelfReferenceMatch(epitope.sequence, status == 'excluded', antigen.kind,
                sources=tuple(dict.fromkeys(sources)), source_provenance_complete=True,
                reference_sha256=bundle.identity, genome_release=reference.get('reference_id', '')),
            input_evidence=tuple(r for r in epitope.input_evidence if r.get('schema') != SCHEMA) + (evidence,))
    dataset.epitopes = tuple(updated.get(e.prediction_id, e) for e in dataset.epitopes)
    dataset.selection.setdefault('source_evidence', {})[bundle.identity] = dict(
        schema='vaxrank.oncoref_target_reference.v1', manifest=manifest, coverage=bundle.coverage)
    dataset.selection['source_evidence'].update(_targeting_summary(updated.values(), dataset.antigens, bundle.identity))
    dataset.result.df['cta_target_status'] = dataset.result.df.prediction_id.map(
        {key: targeting_evidence(e)['status'] for key, e in updated.items()})
    dataset.result.df['cta_target_reference'] = dataset.result.df.prediction_id.map(
        {key: bundle.identity for key in updated})
    for field in ('occurs_in_reference', 'occurs_in_non_CTA_reference'):
        values = dataset.result.df.prediction_id.map({key: getattr(e, field) for key, e in updated.items()})
        dataset.result.df[field] = values.where(values.notna(), dataset.result.df.get(field))
