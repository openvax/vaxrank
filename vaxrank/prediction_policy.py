"""Describe effective prediction stages without loading or executing models."""

from dataclasses import asdict

import pandas as pd

from .processing_prediction import ProcessingPrediction


def record_candidate_models(args, models):
    """Retain Topiary's configured model facts separately from observed rows."""
    from topiary import TopiaryPredictor, describe_prediction_models
    predictor = (models if isinstance(models, TopiaryPredictor) else
                 TopiaryPredictor(models=models if isinstance(models, (list, tuple)) else [models]))
    inventory = describe_prediction_models(predictor)
    runs = getattr(args, '_candidate_model_runs', [])
    if inventory not in runs:
        runs.append(inventory)
    args._candidate_model_runs = runs


def retained_processing(datasets):
    """Merge validated occurrence annotations; conflicting evidence is explicit."""
    records = {}
    for dataset in datasets:
        for values in dataset.selection.get('processing_annotations', ()):
            record = ProcessingPrediction(**values)
            if record.key() in records and records[record.key()] != record:
                raise ValueError('Inputs contain conflicting processing annotations for one occurrence')
            records[record.key()] = record
    return records


def _models(predictions):
    records = []
    for prediction in predictions:
        record = dict(model=prediction.predictor_name or None,
                      version=prediction.predictor_version or None, endpoint=str(prediction.kind))
        if record not in records:
            records.append(record)
    return records


def effective_prediction_configuration(args, ranked, predictions, products, processing):
    """Record actual evidence and execution status, keeping stage semantics apart."""
    dataset = getattr(args, '_design_dataset', None)
    epitopes = (dataset.epitopes if dataset is not None else predictions) or tuple(
        epitope for _, peptides in ranked for peptide in peptides for epitope in peptide.epitopes)
    leaves = [prediction for epitope in epitopes for prediction in epitope.predictions_flat()]
    executed = getattr(args, '_candidate_prediction_mode', None)
    from .epitope_config import EpitopeConfig
    from .epitope_dsl import default_score_expr, resolve_default_methods, resolve_default_versions
    cfg = (dataset.config if dataset is not None and dataset.config is not None
           else EpitopeConfig(**getattr(args, '_effective_epitope_configuration', {})))
    policy = (dict(mode='selection_policy', definition=cfg.selection_policy)
              if cfg.selection_policy else dict(
                  mode='legacy DSL translation', filter=cfg.filter_expr,
                  score=cfg.score_expr or default_score_expr(cfg),
                  min_epitope_score=cfg.min_epitope_score, default_methods=cfg.default_methods))
    if not cfg.selection_policy:
        identities = pd.DataFrame([
            dict(kind=model['endpoint'], prediction_method_name=model['model'], predictor_version=model['version'])
            for model in _models(leaves)])
        policy['unqualified_methods'] = resolve_default_methods(cfg, identities)
        policy['unqualified_versions'] = [dict(kind=kind, model=model, version=version)
                                           for (kind, model), version in resolve_default_versions(cfg, identities).items()]
    candidate = dict(
        status=executed or 'recorded evidence; no new inference', models=_models(leaves),
        configured_model_runs=getattr(args, '_candidate_model_runs', getattr(
            args, '_recorded_prediction_configuration', {}).get('candidate', {}).get(
                'configured_model_runs', [])),
        alleles=sorted({p.allele for p in leaves if p.allele}),
        lengths=sorted({len(e.sequence) for e in epitopes}),
        context=('peptide_only' if getattr(args, 'external_peptide_only', False)
                 else 'source occurrences' if executed else 'recorded occurrence context'),
        cache=(dict(path=args.prediction_cache,
                    behavior='exact context/allele/kind lookup; version-checked fallback',
                    hit_counts='not recorded')
               if executed and getattr(args, 'prediction_cache', None) else None),
        coverage=('validated by Topiary for supplied occurrences' if executed in ('fresh', 'additive')
                  else 'recorded rows; completeness not inferred'),
        policy=policy,
        scoring=('selection_policy: ' + str(cfg.selection_policy.get('name', 'unnamed'))
                 if cfg.selection_policy else 'score = %s; filter = %s; minimum score = %s' % (
                     policy['score'], policy['filter'] or 'none', policy['min_epitope_score'])))
    if dataset is not None:
        candidate['additive_runs'] = dataset.result.extra.get('candidate_rescoring', {})
        for run in candidate['additive_runs'].values():
            for feature in run.get('features', {}).values():
                model = dict(model=feature['model'], version=feature.get('version'), endpoint=feature['kind'])
                if model not in candidate['models']:
                    candidate['models'].append(model)
    junction_records = [construct.elements.get('junction_swap', {})
                        for construct in products.get('mrna', ())]
    evidence = [record['prediction'] for record in junction_records if record.get('prediction')]
    junction = dict(
        status=('replayed evidence; no new inference' if evidence and getattr(args, '_junction_replayed', False)
                else 'executed' if evidence else 'skipped'),
        reason=None if evidence else 'no enabled junction queries in emitted products',
        context='peptide_only', cache=None, runs=[{
            key: value for key, value in record.items() if key != 'requests'
        } for record in evidence])
    processing_models = []
    for prediction in processing.values():
        model = dict(model=prediction.predictor_name, version=prediction.predictor_version,
                     endpoint='proteasome_cleavage')
        if model not in processing_models:
            processing_models.append(model)
    annotation = dict(
        status=('executed' if processing and getattr(args, 'processing_aware_annotation', False)
                else 'replayed evidence; no new inference' if processing else 'skipped'),
        reason=(None if processing else 'disabled' if not getattr(args, 'processing_aware_annotation', False)
                else 'enabled but no annotations returned: model unavailable or no assessable occurrences; see run log'),
        models=processing_models, alleles=[],
        lengths=sorted({len(p.peptide_sequence) for p in processing.values()}),
        context='exact source sequence and peptide offset',
        coverage=dict(annotated_occurrences=len(processing),
                      assessed_c_boundaries=sum(p.c_term_cleavage_prob is not None
                                                for p in processing.values())),
        scoring='annotation only: sqrt(c_term * (1 - max_internal)); not a calibrated probability',
        configuration=(dict(human_only=getattr(args, 'pepsickle_human_only', False),
                            threshold=getattr(args, 'pepsickle_threshold', 0.5))
                       if getattr(args, 'processing_aware_annotation', False) else getattr(
                           args, '_recorded_prediction_configuration', {}).get(
                               'processing', {}).get('configuration')),
        cache='retained native annotations' if processing else None)
    configuration = dict(candidate=candidate, junction=junction, processing=annotation)
    if dataset is not None:
        dataset.selection['processing_annotations'] = [asdict(p) for p in processing.values()]
    return configuration


def prediction_summary(configuration):
    """Compact readable rows shared by the text and rendered reports."""
    rows = []
    for name, stage in configuration.items():
        models = stage.get('models', [])
        if name == 'junction':
            models = [dict(model=m.get('predictor_name'), version=m.get('predictor_version'),
                           endpoint=m.get('kind'))
                      for run in stage.get('runs', []) for m in run.get('models', [])]
        labels = list(dict.fromkeys('%s [%s]: %s' % (
            model['model'] or 'unknown model', model['version'] or 'unknown version', model['endpoint'])
            for model in models))
        details = ['Models: ' + ('; '.join(labels) or 'none recorded'),
                   'Context: ' + stage['context']]
        inventories = (stage.get('configured_model_runs', []) if name == 'candidate' else
                       [run.get('configured_models', []) for run in stage.get('runs', [])])
        for inventory in inventories:
            for model in inventory:
                details.append('Configured %s: %s [%s]; alleles: %s; lengths: %s; flanks: %s' % (
                    model['model_class'], model['prediction_method_name'] or 'unknown method',
                    model['predictor_version'] or 'unknown version',
                    ', '.join(model['alleles']) if model['alleles'] is not None else 'unknown',
                    model['default_peptide_lengths'], model['protein_flank_lengths']))
        if name == 'junction':
            for run in stage.get('runs', []):
                details.append('Rank: %s; strong ≤ %s%%; mild ≤ %s%%; alleles: %s; lengths: %s; coverage: %s' % (
                    run.get('rank_kind', run.get('rank_kind_policy', 'unknown')),
                    run.get('rank_strong', 'unknown'), run.get('rank_mild', 'unknown'),
                    ', '.join(run.get('alleles', [])), ', '.join(map(str, run.get('peptide_lengths', []))),
                    run.get('coverage', 'not recorded')))
        else:
            coverage = stage['coverage']
            if isinstance(coverage, dict):
                coverage = '%d annotated occurrences; %d assessed C-boundaries' % (
                    coverage['annotated_occurrences'], coverage['assessed_c_boundaries'])
            details += ['Recorded alleles: ' + (', '.join(stage['alleles']) or 'none'),
                        'Lengths: ' + (', '.join(map(str, stage['lengths'])) or 'none'),
                        'Coverage: ' + coverage, 'Scoring: ' + stage['scoring']]
            if name == 'candidate' and stage['policy']['mode'] != 'selection_policy':
                policy = stage['policy']
                if policy.get('unqualified_methods'):
                    details.append('Default models for unqualified DSL references: ' + '; '.join(
                        '%s: %s' % item for item in policy['unqualified_methods'].items()))
                if policy.get('unqualified_versions'):
                    details.append('Default versions for unqualified DSL references: ' + '; '.join(
                        '%s/%s: %s' % (item['kind'], item['model'], item['version'])
                        for item in policy['unqualified_versions']))
        cache = stage.get('cache')
        if isinstance(cache, dict):
            cache = '%s (%s; hit counts %s)' % (cache['path'], cache['behavior'], cache['hit_counts'])
        details.append('Cache: ' + (cache or 'none used'))
        if stage.get('reason'):
            details.append(stage['reason'])
        if name == 'processing' and stage.get('configuration'):
            details.append('Human-only model: %s; cleavage threshold: %s' % (
                stage['configuration']['human_only'], stage['configuration']['threshold']))
        rows.append(dict(stage=name, status=stage['status'], details=details))
    return rows
