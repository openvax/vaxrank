"""Original input-stage observations, preserved separately from current ranking."""

from .selection_policy import canonical_digest
from .version import __version__

from collections import Counter
from dataclasses import asdict

from varcode import StructuralVariant, Variant

from .antigen_report import metadata_values
from .coding_context import coding_context_label, coding_contexts
from .native_serialization import to_native_json


FAMILIES = {
    'mutation': 'Small variants', 'variant_orfs': 'Additional variant ORFs',
    'sv': 'Structural variants', 'fusion': 'Fusions', 'isoforms': 'Isoforms',
    'cta': 'CTAs', 'viral': 'Viral sequences', 'erv': 'ERVs', 'splice': 'Splice antigens',
}


def report_source_key(source):
    # Same identity as the shared ranker's observation grouping. Input scope is
    # validated before pooling; duplicate observations do not become new events.
    return canonical_digest([type(source).__name__, str(source)])


def _source_row(source):
    if isinstance(source, Variant):
        return dict(identity=source.short_description, gene='', context='',
                    family='sv' if isinstance(source, StructuralVariant) else 'mutation',
                    admitted=None, reasons=[])
    families = metadata_values(source).get('isovar_family', ())
    family = families[0] if len(set(families)) == 1 else source.kind.lower()
    family = {'isoform': 'isoforms'}.get(family, family)
    return dict(identity=source.display_identifier, gene=source.display_gene_name,
                context=' / '.join((*source.protein_ids, *source.transcript_ids)),
                family=family,
                admitted=source.tumor_specificity.admits_construct,
                reasons=[source.tumor_specificity.rationale_code]
                    if not source.tumor_specificity.admits_construct else [])


def ranking_coverage(dataset, reports, outcomes, ranked):
    """Describe recorded sources and current ranking without re-running a filter.

    Catalogs include held-out/no-prediction antigens; ranker entries include
    direct variants and the existing external adapters' resolved sources. The
    inventory is explicitly not the unknown upstream discovery universe.
    """
    from .antigen_input import AntigenDataset, _antigen_from_dict

    rows = {}
    retained = Counter()
    target_overlap = Counter()
    ranked_keys = {report_source_key(source): i for i, (source, _) in enumerate(ranked, 1)}

    def add(source, input_ids=()):
        key = report_source_key(source)
        row = rows.setdefault(key, dict(_source_row(source), inputs=set(), has_window=False))
        if not isinstance(source, Variant):
            contexts = coding_contexts(source, dataset.selection.get('source_evidence'))
            if contexts or ('isovar_family', 'variant_orfs') in source.source_metadata:
                row['context'] = ' / '.join((coding_context_label(contexts),
                    *source.protein_ids, *source.transcript_ids))
        row['inputs'].update(input_ids)
        return row

    # Catalogs retain distinct assertions, even when sequences happen to match.
    for report, outcome in zip(reports, outcomes):
        original = report.dataset.provenance if report.dataset is not None else (report.input_provenance,)
        input_ids = [p.input_id for p in original]
        for exclusion in outcome.excluded_sources:
            family = exclusion['family'].lower()
            family = {'snv': 'mutation', 'indel': 'mutation', 'cta/self': 'cta'}.get(family, family) or 'unresolved'
            key = 'excluded:' + canonical_digest(exclusion)
            row = rows.setdefault(key, dict(identity=exclusion['identity'], family=family,
                gene='', context='', inputs=set(), has_window=False, admitted=False,
                recorded_stage=exclusion['stage'], reasons=[exclusion['reason']]))
            row['inputs'].update(input_ids)
        if report.dataset is not None:
            native = report.dataset
            origins = {}
            for input_id, records in native.selection.get('antigen_inputs', {}).items():
                for value in records.values():
                    origins.setdefault(to_native_json(_antigen_from_dict(value)), set()).add(input_id)
            fragment_antigens = {to_native_json(native.antigens[k]) for k in native.mutation_fragments
                                if k in native.antigens}
            for antigen in AntigenDataset.from_epitope_dataset(native).antigens.values():
                if to_native_json(antigen) not in fragment_antigens:
                    add(antigen, origins.get(to_native_json(antigen), ()))
                    retained.setdefault(report_source_key(antigen), 0)
            for epitope in native.epitopes:
                group = epitope.prediction_group_source
                source = (native.mutation_fragments[group].variant if group in native.mutation_fragments
                          else native.antigens.get(group))
                if source is not None:
                    key = report_source_key(source)
                    row = add(source, [epitope.input_provenance.input_id] if epitope.input_provenance else ())
                    if group in native.antigens:
                        row['admitted'] = native.antigens[group].tumor_specificity.admits_construct
                    if group in native.mutation_fragments:
                        row['gene'] = native.mutation_fragments[group].gene_name
                    retained[key] += 1
                    target_overlap[key] += bool(epitope.is_targetable)
                elif not native.construction_reports:
                    key = 'unresolved:' + group
                    rows.setdefault(key, dict(identity=group, gene='', context='', family='unresolved',
                        admitted=None, inputs=set(input_ids), has_window=False, recorded_stage='Source resolution',
                        reasons=['Retained prediction group has no resolved source antigen']))
                    retained[key] += 1
        for entry in outcome.entries:
            row = add(entry.ranking_source, input_ids if report.dataset is None else ())
            row['has_window'] |= bool(entry.ranking_peptides)
            if entry.ranking_peptides:
                row['gene'] = entry.ranking_peptides[0].antigen.display_gene_name
    for key, row in rows.items():
        row['source_key'] = key
        row['inputs'] = sorted(row['inputs'])
        row['rank'] = ranked_keys.get(key)
        # Zero retained rows is a fact only for a known native source catalog.
        row['retained_epitopes'] = retained.get(key)
        if row.get('recorded_stage'):
            row['stage'] = row['recorded_stage']
        elif row['admitted'] is False:
            row['stage'] = 'Construct admission'
        elif row['rank'] is not None:
            row['stage'] = 'Ranked'
        elif row['has_window']:
            row['stage'] = 'Shared ranking'
            row['reasons'] = ['Window available but absent from the final ranked list']
        else:
            row['stage'] = 'Candidate windows'
            row['reasons'] = ['No ranked window under the effective configuration']
            if key in retained and retained[key] == 0:
                row['reasons'].append('No retained epitope predictions')
            elif key in retained and not target_overlap[key]:
                row['reasons'].append('No retained target-overlapping epitope')
    return dict(sources=list(rows.values()),
                original_inputs=[asdict(p) for p in dataset.provenance],
                current_inputs=[asdict(r.input_provenance) for r in reports])


def coverage_view(coverage, source_evidence, ranked_candidates, products):
    """Attach actual product outcomes to recorded input/ranking observations."""
    coverage = coverage or {}
    rows = [dict(row) for row in coverage.get('sources', ())]
    all_products = [p for group in products.values() for p in group]
    mappings_complete = all(p.antigen_names and set(p.antigen_sources) == set(p.antigen_names)
                            and all(p.antigen_sources.values())
                            for p in all_products)
    input_labels = {p['input_id']: 'Input %d' % (i + 1)
                    for i, p in enumerate(coverage.get('original_inputs', ()))}
    current_ranking = {row['source_key']: row for row in ranked_candidates}
    for row in rows:
        row['origin'] = '; '.join(input_labels.get(identity, identity) for identity in row['inputs']) or 'Origin not recorded'
        row['selected'] = None
        current = current_ranking.get(row.get('source_key'))
        if current is not None:
            row['rank'] = current['rank']
            row['product_outcome'] = current['products']
            if products and mappings_complete:
                row['selected'] = bool(current['selected_products'])
        elif row['rank'] is not None:
            row['rank'] = None
            row['stage'] = 'Current report ranking'
            row['product_outcome'] = 'Current source mapping unavailable'
            row['reasons'] = row['reasons'] + ['Previously ranked source absent from the current report ranking']
        else:
            row['product_outcome'] = 'Not ranked'
            if products and mappings_complete:
                row['selected'] = False

    observations, producer_rows, expression_rows = [], [], []
    supplied = {row['family'] for row in rows}
    for evidence_id, evidence in source_evidence.items():
        if 'direct_input_coverage' in evidence:
            observations.append(dict(evidence['direct_input_coverage'], evidence_id=evidence_id))
            supplied.add('mutation')
        if 'isovar_result' in evidence:
            result = evidence['isovar_result']
            family = evidence['family']
            supplied.add(family)
            models = result.get('models', ())
            producer_rows.append(dict(identity=evidence_id, family=FAMILIES.get(family, family),
                sample=result.get('sample_id') or result.get('provenance', {}).get('sample_id'),
                candidates=len(evidence.get('decisions', ())), models=len(models) if 'models' in result else None,
                rna_reconstructions=len(result['rna_candidates']) if 'rna_candidates' in result else None,
                event=result.get('variant') or result.get('event_id'),
                state=result.get('assembly_state') or result.get('status'),
                model_outcomes=[dict(identity=m.get('model_id') or m.get('id') or str(i),
                    reference=m.get('reference', {}).get('orf_id') or m.get('reference', {}).get('kind'),
                    candidates=len(m['candidate_ids']) if 'candidate_ids' in m else None,
                    compatible_rna=m.get('compatible_rna_candidate_count'),
                    state=m.get('state') or m.get('status'), reason=m.get('reason'))
                    for i, m in enumerate(models, 1)],
                unresolved=[dict(identity=m.get('model_id') or m.get('id') or str(i),
                                 context=' / '.join(str(v) for v in (
                                     m.get('reference', {}).get('orf_id'),
                                     m.get('reference', {}).get('transcript_id')) if v),
                                 reason=m.get('reason'), status=m.get('state') or m.get('status'))
                            for i, m in enumerate(models, 1) if m.get('reason')]))
        if 'cta_expression_admission' in evidence:
            from .cta_expression import load_cta_expression_result
            admission = load_cta_expression_result(evidence['cta_expression_admission'])
            supplied.add('cta')
            expression_rows.append(dict(identity=evidence_id, contract=asdict(admission.input_contract),
                                        decisions=[asdict(d) for d in admission.decisions]))

    families = []
    for family in dict.fromkeys((*FAMILIES, *sorted(supplied))):
        sources = [row for row in rows if row['family'] == family]
        known_admission = all(r['admitted'] is not None for r in sources)
        families.append(dict(name=FAMILIES.get(family, family),
            status='Supplied' if family in supplied else 'No supplied evidence recorded',
            records=len(sources) if coverage else None,
            admitted=sum(r['admitted'] for r in sources) if known_admission and sources else None,
            ranked=sum(r['rank'] is not None for r in sources) if coverage else None,
            represented=(sum(r['selected'] for r in sources)
                         if products and mappings_complete and coverage
                         and all(r['selected'] is not None for r in sources) else None)))
    return dict(families=families, sources=rows, direct=observations,
                producers=producer_rows, expression=expression_rows,
                original_inputs=coverage.get('original_inputs', ()),
                current_inputs=coverage.get('current_inputs', ()),
                available=bool(coverage),
                input_count=len({p['input_id'] for p in coverage.get('original_inputs', ())}))


def direct_ranking_coverage(dataset, results, ranked):
    """Adapt already computed direct outcomes to the shared report inventory."""
    from .external_input import ExternalRankingResult, ExternalVariantEntry
    from .external_report import ExternalReport
    provenance, = dataset.provenance
    report = ExternalReport('vcf_bam', provenance.path, input_provenance=provenance, dataset=dataset)
    outcome = ExternalRankingResult(entries=tuple(ExternalVariantEntry(
        variant=result.variant, vaccine_peptides=tuple(results.variant_to_vaccine_peptides_dict.get(result.variant, ())))
        for result in results.isovar_results))
    return ranking_coverage(dataset, [report], [outcome], ranked)


def record_direct_coverage(dataset, results, loaded_variants, *, window_selection_performed=True,
                           filter_thresholds=None):
    """Retain direct VCF/BAM denominators and outcomes before native export.

    These are variant-record counts at the original pipeline stages, not
    antigen counts or a claim of complete peptide prediction coverage. Distinct
    observations of the same input remain distinct catalog records when their
    outcomes differ; pooled reports must not add them as independent inputs.
    """
    provenance, = dataset.provenance
    properties = results.variant_properties()
    for row, result in zip(properties, results.isovar_results):
        row['isovar_filter_results'] = {name: bool(passed) for name, passed in
                                      getattr(result, 'filter_values', {}).items()}
    if not window_selection_performed:
        properties = [dict(row, has_vaccine_peptide=None) for row in properties]
    processed = {result.variant for result in results.isovar_results}
    unavailable = [dict(contig=v.contig, start=v.start, ref=v.ref, alt=v.alt)
                   for v in loaded_variants if v not in processed]
    observation = dict(
        schema='vaxrank.direct_input_coverage.v1', input_id=provenance.input_id,
        producer_version=__version__, unit='parsed variant records',
        stage='original direct pipeline, before shared ranking and product assembly',
        window_selection_performed=window_selection_performed,
        isovar_filter_thresholds=filter_thresholds,
        counts=dict(loaded=len(loaded_variants), annotatable=len(properties),
                    coding_effect=sum(p['is_coding_nonsynonymous'] for p in properties),
                    coding_with_rna_protein=sum(p['is_coding_nonsynonymous'] and p['rna_support']
                                               for p in properties),
                    with_candidate_windows=(sum(p['has_vaccine_peptide'] for p in properties)
                                            if window_selection_performed else None)),
        variant_outcomes=properties, annotation_unavailable=unavailable)
    identity = 'direct-coverage:' + canonical_digest(observation)
    dataset.selection.setdefault('source_evidence', {})[identity] = {
        'direct_input_coverage': observation}
    return observation

