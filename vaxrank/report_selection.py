"""Report recorded selection policy and calculations without running selectors."""

import json
from types import SimpleNamespace

from .combined_score_dsl import combined_score_binding_names, combined_score_bindings
from .epitope_dsl import default_score_expr


UNRECORDED = 'Not recorded'


def display(value):
    if value is None:
        return 'Not set'
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    if isinstance(value, float):
        return '%.6g' % value
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True)
    return str(value)


def _rows(configuration, fields):
    none_labels = {'max_ranked_sources': 'No source limit',
                   'num_target_epitopes_to_keep': 'No additional cap',
                   'ranking_rules': 'Default order; not expanded in this record',
                   'window_selection': 'Disabled', 'default_methods': 'Resolved per input'}
    return [(label, (none_labels.get(key, 'Not set') if configuration[key] is None
                     else display(configuration[key])) if key in configuration else UNRECORDED)
            for key, label in fields]


def selection_policy_view(selection):
    """Use actual assembly settings and resolved policy, never invocation defaults."""
    record = selection.get('run_configuration') or {}
    epitopes = record.get('effective_epitopes') or {}
    vaccine = record.get('effective_vaccine_peptides') or {}
    policy = epitopes.get('selection_policy') or {}
    row_filter = policy.get('filter_by') if policy else epitopes.get('filter_expr')
    score = policy.get('score_by') or epitopes.get('score_expr')
    if score is None:
        required = ({'scoring_mode', 'percentile_rank_cutoff'}
                    if epitopes.get('scoring_mode') == 'percentile_rank' else
                    {'scoring_mode', 'binding_affinity_cutoff', 'logistic_epitope_score_midpoint',
                     'logistic_epitope_score_width'})
        if required <= epitopes.keys():
            score = default_score_expr(SimpleNamespace(**epitopes))
    sections = [dict(title='Epitope scoring and eligibility', rows=[
        ('Per-allele score expression', display(score) if score is not None else UNRECORDED),
        ('Selection policy criteria', display(policy.get('criteria', {})) if policy else 'No named policy criteria recorded'),
        ('Explicit row filter', display(row_filter) if row_filter is not None
         else 'None; see the recorded score expression and minimum score' if epitopes else UNRECORDED),
        ('Minimum score', display(policy.get('min_score')) if policy else
         display(epitopes['min_epitope_score']) if 'min_epitope_score' in epitopes else UNRECORDED),
    ] + (_rows(policy, [('unknown', 'Unknown criterion outcome'),
                       ('ranking_by', 'Policy ranking terms'), ('strata', 'Policy strata')]) if policy else [])
      + _rows(epitopes, [('default_methods', 'Default methods'),
                        ('comparison_reference', 'Reference comparator policy'),
                        ('allele_free_evidence', 'Credit for evidence without an allele'),
                        ('allele_selection_axis', 'Allele selection axis')]))]
    if epitopes.get('comparison_reference') == 'patient_background':
        sections[0]['rows'] += _rows(epitopes, [
            ('patient_background_baseline', 'Patient-background somatic reversion policy')])
    sections.append(dict(title='Candidate window settings', rows=_rows(vaccine, [
        ('combined_score_expr', 'Configured combined-score expression'),
        ('min_peptide_length', 'Minimum candidate length (aa)'),
        ('preferred_peptide_length', 'Preferred candidate length (aa)'),
        ('max_peptide_length', 'Maximum candidate length (aa)'),
        ('padding_around_mutation', 'Padding around the altered region (aa)'),
        ('score_fraction_of_best', 'Score fraction retained before window tie breaks'),
        ('max_vaccine_peptides_per_variant', 'Maximum windows per source'),
        ('num_target_epitopes_to_keep', 'Maximum target epitopes per window'),
        ('ranking_rules', 'Window ranking rules'),
        ('require_target_epitopes_in_variant', 'Require a target epitope'),
        ('window_selection', 'Shared window-selection policy')])))
    sections.append(dict(title='Design scope', rows=_rows(record.get('effective_design', {}), [
        ('max_ranked_sources', 'Maximum ranked sources'), ('vaccine_types', 'Requested modalities')])))
    assembly = selection.get('construct_assembly') or {}
    modalities = {}
    for modality in dict.fromkeys((*record.get('effective_constructs', {}), *assembly)):
        configuration = assembly.get(modality, {}).get('configuration')
        origin = 'Saved assembly configuration'
        if configuration is None:
            configuration = record.get('effective_constructs', {}).get(modality, {})
            origin = 'Resolved construct configuration'
        fields = [('antigen_content', 'Antigen content'), ('candidates_per_slot', 'Candidate windows per source'),
                  ('epitopes_per_antigen', 'Epitopes per antigen'),
                  ('min_antigen_length_aa', 'Minimum refined antigen length (aa)' if configuration.get('window_selection') is not None
                   else 'Minimum antigen length warning threshold (aa)'),
                  ('max_antigen_length_aa', 'Maximum antigen length (aa)'),
                  ('antigens_per_construct', 'Antigen slots per construct'),
                  ('max_constructs', 'Maximum constructs'), ('target_alleles', 'Alleles used for coverage selection')]
        fields += ([('window_selection', 'Peptide window-selection policy')]
                   if modality == 'peptide' else [('max_length_nt', 'Maximum mRNA length (nt)'),
                       ('optimize_linkers', 'Optimize linkers'), ('junction_rank_kind', 'Junction prediction endpoint')])
        modalities[modality] = dict(origin=origin, rows=_rows(configuration, fields))
    return dict(available=bool(record or assembly), name=record.get('name') or 'No named policy recorded',
                identity=record.get('effective_sha256') or UNRECORDED,
                version=record.get('vaxrank_version') or UNRECORDED,
                sections=sections, modalities=modalities)


def candidate_calculation(peptide):
    """Show the actual expression and bound values, including missing-count behavior."""
    values = combined_score_bindings(peptide)
    fragment = peptide.mutant_protein_fragment
    terms = []
    for name in sorted(combined_score_binding_names(peptide.combined_score_expr)):
        value, unit, evidence = values.get(name), 'score', ''
        if name.startswith(('n_rna_', 'n_alt_', 'n_ref_', 'n_overlapping_')):
            unit = (getattr(fragment, 'rna_evidence_subject', '') or 'count unit not recorded'
                    if name.startswith('n_rna_') else 'reads')
            if fragment is not None and getattr(fragment, name, None) is None:
                evidence = 'Source value unrecorded; the score expression receives zero'
        elif name.startswith('n_dna_'):
            unit = 'DNA observations'
        elif name == 'n_mutant_amino_acids':
            unit = 'aa'
        terms.append(dict(name=name, value=value, display=display(value), unit=unit, evidence=evidence))
    return dict(expression=peptide.combined_score_expr, terms=terms,
                result=display(peptide.combined_score), raw_result=peptide.combined_score,
                target_score=display(peptide.target_epitope_score),
                targets=[dict(sequence=e.sequence, offset=e.offset,
                    allele_scores=dict(e.per_allele_scores), score=e.epitope_score) for e in peptide.target_epitopes])


def score_calculations(ranked, rows):
    """Share identical calculations while retaining every linked report source."""
    groups, by_value = [], {}
    peptides = [peptides[0] for _, peptides in ranked if peptides]
    for peptide, row in zip(peptides, rows):
        calculation = candidate_calculation(peptide)
        key = json.dumps(calculation, sort_keys=True)
        if key not in by_value:
            by_value[key] = len(groups)
            groups.append(dict(calculation, sources=[]))
        groups[by_value[key]]['sources'].append(row['source_id'])
    return groups


def window_view(audit, sequence):
    """Summarize recorded alternatives, with no retrospective selection."""
    if not audit:
        return None
    alternatives = audit.get('alternatives', ())
    # New records have exact ranks. Historical records retain producer order;
    # without ordering keys we cannot call that order a tie-break decision.
    ordered = sorted(alternatives, key=lambda a: (
        a.get('selection_rank') is None, a.get('selection_rank') or 0))
    decision = audit.get('decision')
    others = [a for a in ordered if (a.get('input_order') != decision.get('input_order')
              if decision else a.get('sequence') != sequence)]
    # Show nearest ranked alternatives and one rejected window if available.
    preview = others[:3]
    rejected = next((a for a in others if not a.get('eligible')), None)
    if rejected is not None and rejected not in preview:
        preview.append(rejected)
    fields = ('sequence', 'eligible', 'ineligibility_reason', 'target_score', 'utility', 'combined_score',
              'ordering_key', 'selection_rank', 'selected', 'decision_reason', 'parent_window')
    policy = audit.get('policy', {})
    utility = audit.get('algorithm') != 'score_fraction_then_ranking'
    return dict(sequence=sequence, parent=audit.get('parent_window'),
        utility_selector=utility,
        ordering_description=(UNRECORDED if not audit.get('decision') else
                              'negative utility, preferred-length distance, ranking tuple, stable input order'
                              if utility else 'score-fraction eligibility, then ranking tuple and stable input order'),
        none_eligible=audit.get('selected') == [],
        decision=audit.get('decision'), policy=policy,
        metrics=_rows(audit, [('target_score', 'Unique sequence/allele target score' if utility else 'Target-occurrence score'),
            ('combined_score', 'Combined score'),
            ('non_cta_self_score', 'Recorded non-CTA self score'), ('at_risk_target_score', 'At-risk target score'),
            ('utility', 'Window utility'), ('preferred_length', 'Preferred length (aa)'),
            ('limit', 'Maximum windows retained'), ('self_provenance_complete', 'Self-reference provenance complete')]),
        serum_status=('Not evaluated by the score-fraction selector' if not utility else
                      'Skipped: serum weight is zero' if policy.get('serum_weight') == 0 else
                      'Unassessed models: ' + ', '.join(audit['unassessed_models'])
                      if audit.get('unassessed_models') else
                      'Recorded cleavage profiles evaluated' if audit.get('cleavage_profiles') else UNRECORDED),
        self_status=('No optional self-weight penalty; the ranking tuple can use self-epitope score'
                     if not utility else 'No self penalty applied; this is not a completed self assessment'
                     if policy.get('self_weight') == 0 else 'Self penalty active'
                     if policy.get('self_weight') else UNRECORDED),
        targets=audit.get('target_details', []),
        alternatives=[{key: a[key] for key in fields if key in a} for a in preview],
        alternative_count=len(alternatives))


def junction_view(product):
    """Explain the recorded construct-level comparison and each actual linker."""
    record = product.elements.get('junction_swap')
    if record is None:
        return dict(status=UNRECORDED, rows=[], linkers=[], decisions=[])
    enabled = record.get('enabled')
    status = ('Optimization evaluated' if enabled else 'Not assessed: ' + record.get('note', UNRECORDED))
    rows = _rows(record, [('note', 'Decision'), ('policy', 'Model selection policy'),
        ('burden_strong', 'Chosen strong burden'), ('burden_mild', 'Chosen mild burden'),
        ('default_burden_strong', 'Default strong burden'), ('default_burden_mild', 'Default mild burden')]) if enabled else []
    prediction = record.get('prediction') or {}
    if enabled:
        rows += _rows(prediction, [('rank_kind', 'Prediction endpoint'), ('rank_metric', 'Rank metric'),
            ('rank_strong', 'Strong threshold (%)'), ('rank_mild', 'Mild threshold (%)'),
            ('alleles', 'Evaluated alleles'), ('coverage', 'Prediction coverage'),
            ('scored_prediction_count', 'Scored predictions')])
    return dict(status=status, rows=rows,
                decisions=prediction.get('junction_decisions', []),
                linkers=product.elements.get('linkers_per_junction', []))


def attach_selection_decisions(data, selection, ranked, products):
    """Link compact decision summaries to existing report sources and products."""
    from .assembly_antigens import source_matches_candidate
    from .input_coverage import report_source_key
    from .native_serialization import from_native_json, to_native_json
    source_ids = {row['source_key']: row['source_id'] for row in data['ranked_candidates']}
    candidates = [(source_ids[report_source_key(source)], peptide)
                  for source, peptides in ranked if peptides for peptide in peptides]
    assembly_records = selection.get('construct_assembly', {})
    graphs = {modality: json.loads(record['source_graph']) for modality, record in assembly_records.items()
              if 'source_graph' in record}
    recorded_windows = {}
    for graph in graphs.values():
        for source in graph['inputs']:
            for candidate in source['candidates']:
                key = (json.dumps(candidate['antigen'], sort_keys=True), candidate['amino_acids'])
                audit = candidate.get('window_selection', {})
                serialized = json.dumps(audit, sort_keys=True)
                recorded_windows.setdefault(key, {})[serialized] = audit
    groups, by_value = [], {}

    def register(audit, sequence, source_ids=(), product=None):
        view = window_view(audit, sequence)
        if view is None:
            return None
        key = json.dumps(view, sort_keys=True)
        if key not in by_value:
            by_value[key] = len(groups)
            groups.append(dict(view, identifier='W%d' % (len(groups) + 1), sources=[], products=[]))
        group = groups[by_value[key]]
        group['sources'] = sorted(set(group['sources']) | set(source_ids), key=lambda source: int(source[1:]))
        if product and product not in group['products']:
            group['products'].append(product)
        return group['identifier']

    for row in data['ranked_candidates']:
        row['assembly_outcomes'] = []
        identifiers = []
        for source, peptide in candidates:
            if source != row['source_id']:
                continue
            key = (json.dumps(json.loads(to_native_json(peptide.antigen)), sort_keys=True), peptide.amino_acids)
            # A current reranking can produce new audit details. Historical
            # product reports use only decisions stored with those products.
            audits = ([from_native_json(json.dumps(audit), dict)
                       for audit in recorded_windows.get(key, {}).values()] if graphs else
                      [peptide.window_selection_audit])
            identifiers.extend(identifier for audit in audits
                               if (identifier := register(audit, peptide.amino_acids, [source])))
        row['window_decisions'] = list(dict.fromkeys(identifiers))
    for group in data['product_groups']:
        assembly = selection.get('construct_assembly', {}).get(group['modality'], {})
        recorded = graphs.get(group['modality'], {}).get('assembly_decisions')
        decisions = None
        if recorded is not None:
            decisions = dict(order=recorded['order'], rows=[])
            for item in recorded['antigens']:
                ids = sorted({identifier for source in item['sources'] for identifier, peptide in candidates
                    if source_matches_candidate(source, peptide.antigen, peptide.amino_acids,
                        {e.prediction_group_source for e in peptide.epitopes},
                        {e.input_provenance.input_id for e in peptide.epitopes if e.input_provenance})},
                    key=lambda source: int(source[1:]))
                decisions['rows'].append(dict(order=item['order'], sequence=item['sequence'],
                    sources=ids, mapping_count=len(item['sources']), products=item['products'], outcome=item['outcome']))
                for candidate in data['ranked_candidates']:
                    if candidate['source_id'] in ids:
                        candidate['assembly_outcomes'].append(dict(
                            modality=group['modality'], outcome=item['outcome'],
                            shared=len(item['sources']) > 1))
        group['assembly_decisions'] = decisions
        for audit in assembly.get('window_audit') or []:
            if audit.get('selected') != []:
                continue
            sources = [a.get('source') for a in audit.get('alternatives', ()) if a.get('source')]
            ids = [identifier for identifier, peptide in candidates if any(
                json.dumps(source, sort_keys=True) == json.dumps(
                    {key: getattr(peptide.antigen, key) for key in source}, sort_keys=True)
                for source in sources)]
            identifier = register(audit, None, ids)
            for candidate in data['ranked_candidates']:
                if candidate['source_id'] in ids:
                    if identifier not in candidate['window_decisions']:
                        candidate['window_decisions'].append(identifier)
                    candidate['assembly_outcomes'].append(dict(
                        modality=group['modality'], outcome='window_eligibility', shared=False))
        for row, product in zip(group['products'], products[group['modality']]):
            ids = [identifier for a in row['antigens'] for identifier in a['source_ids']]
            window = register(product.components.get('window_selection'),
                              product.sequence, ids, product.name)
            row['window_decisions'] = ([window] if window else sorted({identifier
                for candidate in data['ranked_candidates'] if candidate['source_id'] in ids
                for identifier in candidate['window_decisions']}))
            row['junction'] = junction_view(product) if group['modality'] == 'mrna' else None
    data['window_decisions'] = groups
