"""Combine external reports and select historical or fresh prediction evidence.

The candidate set is the reported peptide occurrences, not all windows in the
source proteins. Source identities and biological admission stay with the
readers; mhctools supplies new predictions and Topiary supplies scoring.
"""
from collections import defaultdict
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
from functools import partial

import pandas as pd

from .epitope_dsl import attach_per_allele_scores
from .epitope_io import (
    annotate_credited_alleles, attach_source_annotations, neoepitope_core_row,
    normalize_hla_allele, read_lens_report, read_pvacseq_report, save_predictions,
)
from .external_report import ExternalRecord
from .epitope_dataset import EpitopeDataset, read_epitope_report
from .input_scope import (
    InputProvenance, SHARED_FIELDS, combine_scopes, normalize_alleles, provenance_columns,
    read_input_manifest, report_declarations, resolve_input_genome, scope_from_mapping,
    validate_input_scopes,
)


def _exacto_identity(path, options):
    contents = {'table': hashlib.sha256(Path(path).read_bytes()).hexdigest()}
    contents.update({name: hashlib.sha256(Path(options[name]).read_bytes()).hexdigest()
                     for name in ('primary_structures', 'transcript_read_support') if name in options})
    return 'exacto:' + hashlib.sha256(json.dumps(contents, sort_keys=True).encode()).hexdigest()


def read_exacto_report(path, *, scope_declarations=None, **_):
    """Use Topiary's native reader, retaining unknown admission and MHC scope."""
    from topiary import read_exacto
    from .external_report import ExternalReport
    declarations = dict(scope_declarations or {})
    options = declarations.pop('exacto', {})
    scope = scope_from_mapping(declarations, str(path))
    if not scope.sample_id:
        raise ValueError('Native Exacto inputs require an explicit sample_id in --input-manifest')
    result = read_exacto(path, sample_name=scope.sample_id, reference_name=scope.reference_assembly,
                        library_id=scope.library_id, **options)
    # These coordinates/context are producer-supplied only with the companion
    # primary structure. Preserve unknown source coordinates in the table.
    result.df['peptide_offset'] = result.df.peptide_start
    result.df['source_sequence'] = result.df.sequence.where(result.df.peptide_start.notna(), None)
    label = _exacto_identity(path, options) + ':' + hashlib.sha256(
        json.dumps(asdict(scope), sort_keys=True).encode()).hexdigest()
    dataset = EpitopeDataset.from_topiary(result, label=label)
    return ExternalReport('exacto', str(path), report_df=dataset.report_frame(),
                          epitopes=dataset.epitopes, dataset=dataset)


READERS = {
    "lens": read_lens_report, "pvacseq": read_pvacseq_report,
    "topiary": partial(read_epitope_report, source_format="topiary"),
    "epitopes": partial(read_epitope_report, source_format="epitopes"),
    "exacto": read_exacto_report,
}


def external_inputs(args):
    """Return ordered, validated (format, path) pairs from CLI arguments."""
    return [(fmt, path) for fmt, path, _ in external_input_specs(args)]


def external_input_specs(args):
    """Resolve a manifest or the existing single-report flags."""
    inputs = [(fmt, getattr(args, 'input_' + fmt)) for fmt in READERS
              if getattr(args, 'input_' + fmt, None)]
    for value in getattr(args, 'external_input', None) or ():
        fmt, sep, path = value.partition('=')
        if not sep or fmt not in READERS or not path:
            raise ValueError("--external-input requires " + ", ".join(fmt + "=PATH" for fmt in READERS))
        inputs.append((fmt, path))
    manifest = getattr(args, 'input_manifest', None)
    if manifest:
        if inputs:
            raise ValueError("--input-manifest cannot be combined with --input-lens, "
                             "--input-pvacseq or --external-input; list every file in the manifest")
        return read_input_manifest(manifest)
    return [(fmt, path, {}) for fmt, path in inputs]


def _namespace_report(report, provenance):
    if report.dataset is not None:
        dataset = report.dataset
        inherited = tuple(replace(p, scope=combine_scopes(
            p.scope, provenance.scope, report.path, SHARED_FIELDS))
            for p in dataset.provenance)
        candidates = tuple(replace(e, input_provenance=e.input_provenance or provenance)
                           for e in dataset.epitopes)
        dataset = replace(dataset, provenance=inherited or (provenance,), epitopes=candidates)
        return replace(report, dataset=dataset, epitopes=candidates, input_provenance=provenance)
    input_id = provenance.input_id
    scope_columns = provenance_columns(provenance)
    keys = {r.key.identifier: replace(r.key, input_id=input_id)
            for r in report.records}
    ids = {old: key.identifier for old, key in keys.items()}

    def row_copy(row):
        row = dict(row)
        if row.get('prediction_id') in ids:
            row['prediction_id'] = ids[row['prediction_id']]
        row['input_source'] = input_id
        row['input_path'] = report.path
        row['input_format'] = report.source_format
        row.update(scope_columns)
        return row

    frame = report.report_df.copy()
    if not frame.empty:
        frame['Prediction identity'] = frame['Prediction identity'].map(ids)
    frame['Input source'] = input_id
    frame['Input file'] = report.path
    frame['Input format'] = report.source_format
    scoring = report.report_df.attrs['topiary_df'].copy()
    if not scoring.empty:
        scoring['prediction_id'] = scoring['prediction_id'].map(ids)
    scoring['input_source'] = input_id
    scoring['input_path'] = report.path
    scoring['input_format'] = report.source_format
    for name, value in scope_columns.items():
        frame[name] = value
        scoring[name] = value
    frame.attrs = {'topiary_df': scoring}
    # Topiary owns evidence normalization and units. Keep its canonical and
    # source-specific evidence side by side, even when observations disagree.
    from topiary import EVIDENCE_COLUMNS
    from topiary.evidence import source_columns
    evidence_columns = set(EVIDENCE_COLUMNS) | set(source_columns(scoring))
    evidence_columns.update(c for c in scoring if c.startswith(('n_rna_', 'n_dna_')))
    evidence = defaultdict(list)
    for record in report.records:
        row = {name: None if pd.isna(value) else value
               for name, value in record.row.items() if name in evidence_columns}
        if row and row not in evidence[record.key.identifier]:
            evidence[record.key.identifier].append(row)
    return replace(
        report, report_df=frame, input_provenance=provenance,
        epitopes=tuple(replace(e, prediction_id=ids[e.prediction_id],
                              input_provenance=provenance,
                              input_evidence=tuple(evidence[e.prediction_id]))
                       for e in report.epitopes),
        records=tuple(ExternalRecord(keys[r.key.identifier], row_copy(r.row))
                      for r in report.records),
        rows=tuple(row_copy(r) for r in report.rows))


class InputTablePredictor:
    """Replay historical predictions for their exact source occurrences.

    Unlike a live-model cache, this does not infer a predictor version, complete
    a genotype, pool conflicting files, or use one source's scores for another.
    A miss raises. The same reports can instead be sent to fresh prediction.
    """
    def __init__(self, epitopes):
        self._entries = {}
        for epitope in epitopes:
            key = epitope.prediction_group_key
            if not epitope.prediction_id or key in self._entries:
                raise ValueError("Input predictions require unique source occurrence identities")
            self._entries[key] = epitope

    def predict_candidates(self, epitopes):
        result = []
        for epitope in epitopes:
            original = self._entries[epitope.prediction_group_key]
            if (epitope.sequence, epitope.source_sequence, epitope.offset,
                    epitope.n_flank, epitope.c_flank) != (
                    original.sequence, original.source_sequence, original.offset,
                    original.n_flank, original.c_flank):
                raise ValueError("Input prediction occurrence/context does not match")
            result.append(replace(epitope, predictions=original.predictions,
                                  comparators=original.comparators,
                                  patient_alleles=original.patient_alleles,
                                  per_allele_scores={}, allele_attributions=()))
        return result


def _context(peptide):
    # Readers retain the real source window even when their Prediction leaves
    # contain no flanks. Use that window, never a guessed reference extension.
    if peptide.source_sequence:
        start, end = peptide.offset, peptide.offset + len(peptide.sequence)
        if peptide.source_sequence[start:end] != peptide.sequence:
            raise ValueError("Peptide does not match its source context")
        return (peptide.sequence, peptide.source_sequence[:start],
                peptide.source_sequence[end:])
    return peptide.sequence, peptide.n_flank, peptide.c_flank


def rescore_candidates(epitopes, models, alleles, *, use_flanks=True):
    """Predict exact supplied windows with Topiary's public occurrence API."""
    return _predict_selected_occurrences(epitopes, models, alleles, use_flanks=use_flanks)[0]


def _predict_selected_occurrences(epitopes, models, alleles, *, use_flanks=True):
    from dataclasses import fields
    from mhctools.pred import Prediction
    from topiary import predict_peptide_occurrences, is_stated
    from topiary.ranking import prediction_mhc_scope
    from .candidate_epitope import stated_or_blank

    if not models:
        raise ValueError("Fresh rescoring requires predictors")
    alleles = tuple(sorted({normalize_hla_allele(a) for a in alleles}))
    requests, owners = [], []
    for index, epitope in enumerate(epitopes):
        peptides = [(None, epitope), *epitope.comparators.items()]
        for comparator, peptide in peptides:
            if not peptide.sequence:
                continue
            sequence, n, c = _context(peptide)
            request = dict(prediction_id=str(len(requests)), peptide=sequence,
                           peptide_offset=peptide.offset, n_flank=n, c_flank=c)
            if epitope.input_provenance is not None:
                request['sample_name'] = (epitope.input_provenance.scope.sample_id
                                          or epitope.input_provenance.scope.patient_id
                                          or epitope.input_provenance.input_id)
            requests.append(request)
            owners.append((index, comparator))
    leaves = defaultdict(list)
    predicted_frames = []
    prediction_fields = {field.name for field in fields(Prediction)}
    for model in models:
        support = model.kind_support()
        dependent = any(spec.get('mhc_dependence') != 'none' for spec in support.values())
        configured = tuple(sorted({normalize_hla_allele(a) for a in getattr(model, 'alleles', ())}))
        if dependent and (not alleles or configured != alleles):
            raise ValueError("Predictor allele coverage must match the explicit patient genotype")
        predicted = predict_peptide_occurrences(requests, model, use_flanks=use_flanks)
        for row in predicted.to_dict('records'):
            request_index = int(row['prediction_id'])
            index, comparator = owners[request_index]
            values = {key: value for key, value in row.items() if key in prediction_fields}
            values['predictor_name'] = row['prediction_method_name']
            values['offset'] = int(row['peptide_offset'])
            for key in ('value', 'score', 'percentile_rank', 'measurement_context', 'peptide_input', 'cache_key'):
                if key in values and not is_stated(values[key]):
                    values[key] = None
            values['allele'] = stated_or_blank(values.get('allele'))
            if isinstance(values.get('measurement_context'), dict) and not values['measurement_context'].get('unit'):
                # Topiary aliases dimensionless primary scores to value.
                # mhctools.value is a unit-bearing measurement; retain the
                # canonical table value without fabricating a native unit.
                values['value'] = None
            leaves[(index, comparator)].append(Prediction.from_dict(values))
        primary = predicted.loc[[owners[int(identity)][1] is None for identity in predicted.prediction_id]].copy()
        def match_key(row):
            return (owners[int(row['prediction_id'])][0], row['kind'], row['prediction_method_name'],
                    row.get('predictor_version'), prediction_mhc_scope(row['allele'],
                        dependence=row['prediction_mhc_dependence'], allele_set=row['allele_set']))
        wildtypes = {match_key(row): row for row in predicted.to_dict('records')
                     if owners[int(row['prediction_id'])][1] == 'wt'}
        matches = [wildtypes.get(match_key(row), {}) for row in primary.to_dict('records')]
        for field in ('value', 'score', 'affinity', 'percentile_rank', 'prediction_method_name',
                      'predictor_version', 'allele', 'allele_set', 'measurement_context'):
            if wildtypes:
                primary['wt_' + field] = [row.get(field) for row in matches]
        primary['prediction_id'] = [epitopes[owners[int(identity)][0]].prediction_group_source
                                    for identity in primary.prediction_id]
        predicted_frames.append(primary)
    fresh = []
    for index, epitope in enumerate(epitopes):
        fresh.append(replace(
            epitope, predictions=tuple(leaves[(index, None)]),
            comparators={name: replace(peptide, predictions=tuple(leaves[(index, name)]))
                         for name, peptide in epitope.comparators.items()},
            patient_alleles=alleles, per_allele_scores={}, allele_attributions=()))
    # Native leaves preserve comparators; the public prediction frame retains
    # declared MHC dependence, configured genotype and supplied flank context.
    scoring = pd.concat(predicted_frames, ignore_index=True) if predicted_frames else pd.DataFrame()
    return fresh, scoring


def _fresh_report_frame(report, epitopes):
    """Display fresh values separately from each unchanged original row."""
    originals = defaultdict(dict)
    for row in report.report_df.to_dict('records'):
        originals[row['Prediction identity']][row['Allele']] = row
    rows = []
    for epitope in epitopes:
        original_by_allele = originals[epitope.prediction_id]
        metadata = next(iter(original_by_allele.values()))
        for allele in epitope.patient_alleles:
            row = {key: metadata[key] for key in (
                'Prediction identity', 'Peptide offset', 'Source sequence name',
                'Input source', 'Input file', 'Input format', 'Gene name',
                'Genomic variant') if key in metadata}
            row.update(neoepitope_core_row(
                allele, epitope.sequence, None,
                epitope.wt.sequence if epitope.wt else '', None,
                metadata.get('Gene name'), metadata.get('Genomic variant')))
            if report.input_provenance is not None:
                row.update(provenance_columns(report.input_provenance))
            row['Prediction evidence'] = 'fresh'
            # Every historical field is retained under an explicit namespace;
            # no old value is presented as a new-model result or a new allele.
            for key, value in original_by_allele.get(allele, {}).items():
                row['Input ' + key] = value
            for prefix, peptide in (('', epitope), ('WT ', epitope.wt)):
                if peptide is None:
                    continue
                for p in peptide.predictions_flat():
                    if p.allele and p.allele != allele:
                        continue
                    label = '%s%s[%s] %s' % (
                        prefix, p.predictor_name, p.predictor_version or 'version unknown', p.kind)
                    row[label + ' value'] = p.value
                    row[label + ' score'] = p.score
                    row[label + ' percentile rank'] = p.percentile_rank
            rows.append(row)
    return pd.DataFrame(rows)


def _annotate_sequence_matches(frame, epitopes):
    """Link exact sequence matches without merging their source evidence.

    These hashes identify amino-acid strings, not ORFs or biological events.
    A shared short context cannot establish full-length ORF equivalence.
    """
    contexts = {e.prediction_group_key: e.source_sequence or e.sequence for e in epitopes}
    peptides = {e.prediction_group_key: e.sequence for e in epitopes}
    if frame.empty:
        return frame
    frame = frame.copy()
    identity = pd.Series(list(frame[[
        'Prediction identity', 'Mutant peptide sequence', 'Peptide offset'
    ]].itertuples(index=False, name=None)), index=frame.index)
    frame['Reported sequence context'] = identity.map(contexts)
    for label, sequences in (('Context', contexts), ('Peptide', peptides)):
        frame[label + ' sequence SHA256'] = identity.map({
            key: hashlib.sha256(sequence.encode('ascii')).hexdigest()
            for key, sequence in sequences.items()})
    frame['Context extent'] = [
        'epitope only' if contexts[key] == peptides[key] else 'reported window'
        for key in identity]
    return frame


def read_reports(inputs, *, scopes=None, manifest_path=None, genome=None, alleles=()):
    """Load and validate source evidence before any models or scoring run."""
    inputs = list(inputs)
    if scopes is None:
        scopes = [{} for _ in inputs]
    if len(scopes) != len(inputs):
        raise ValueError("Each external input must have its own scope declaration")
    reports, seen = [], set()
    for (fmt, path), declarations in zip(inputs, scopes):
        if fmt not in READERS:
            raise ValueError("Unsupported external format: %s" % fmt)
        content = Path(path).read_bytes()
        content_sha256 = hashlib.sha256(content).hexdigest()
        content_id = (_exacto_identity(path, declarations.get('exacto', {}))
                      if fmt == 'exacto' else fmt + ':' + content_sha256)
        if content_id in seen:
            raise ValueError("The same external input was supplied more than once: %s" % path)
        seen.add(content_id)
        reader_kwargs = {'score_epitopes': False}
        if fmt in ('topiary', 'epitopes', 'exacto'):
            reader_kwargs['scope_declarations'] = declarations
        scope_declarations = {key: value for key, value in declarations.items() if key != 'exacto'}
        report = READERS[fmt](path, **reader_kwargs)
        producer_scope, producer_declarations = report_declarations(
            path, fmt, frame=report.dataset.result.df if report.dataset is not None else None)
        if report.dataset is not None and report.dataset.provenance:
            declared = scope_from_mapping(scope_declarations, str(path))
            stored_scope = validate_input_scopes([
                replace(p, scope=combine_scopes(p.scope, declared, str(path), SHARED_FIELDS))
                for p in report.dataset.provenance], genome=genome)
            producer_scope = combine_scopes(producer_scope, stored_scope, str(path), SHARED_FIELDS)
        scope = combine_scopes(
            producer_scope, scope_from_mapping(scope_declarations, str(path)), str(path))
        # An identical file assigned to a different patient/sample must not
        # acquire the same occurrence IDs in separately saved native outputs.
        scope_id = hashlib.sha256(json.dumps(asdict(scope), sort_keys=True).encode()).hexdigest()
        input_id = content_id + ':' + scope_id
        observed = sorted({a for e in report.epitopes for a in e.patient_alleles})
        provenance = InputProvenance(
            input_id=input_id, source_format=fmt, path=str(path), scope=scope,
            content_sha256=content_sha256,
            manifest_path=str(Path(manifest_path).resolve()) if manifest_path else None,
            manifest_declarations=dict(declarations), report_declarations=producer_declarations,
            observed_mhc_alleles=normalize_alleles(observed, str(path)) if observed else ())
        reports.append(_namespace_report(report, provenance))
    validate_input_scopes([r.input_provenance for r in reports],
                          genome=genome, prediction_alleles=alleles)
    return reports


def report_config(reports, epitope_config=None):
    """Resolve an explicit policy or require agreement between saved policies."""
    if epitope_config is None:
        from .epitope_config import EpitopeConfig
        configs = [r.dataset.config for r in reports
                   if r.dataset is not None and r.dataset.config is not None]
        if configs and any(cfg != configs[0] for cfg in configs):
            raise ValueError("Saved inputs use different scoring policies; supply an explicit epitope config")
        epitope_config = configs[0] if configs else EpitopeConfig()
    return epitope_config


def _add_predictions(reports, models, prefix, use_flanks):
    """Project public additive features through retained consumer identities."""
    from copy import deepcopy
    from topiary import TopiaryResult, combine_sources, reconcile_evidence
    from topiary import rescore_candidates as add_features
    frames, metadata = [], {}
    runs = {}
    for index, report in enumerate(reports):
        source = (report.dataset.result if report.dataset is not None
                  else TopiaryResult(report.report_df.attrs['topiary_df']))
        if 'source_observation_id' not in source.df:
            scope = report.input_provenance.scope
            source = combine_sources({report.input_provenance.input_id: source}, sample_name=(
                scope.sample_id or scope.patient_id or report.input_provenance.input_id))
            frame = source.long_df.copy()
            frame['vaxrank_prediction_id'] = frame.prediction_id
        else:
            frame = source.long_df.copy()
        frame['_vaxrank_report_index'] = index
        frames.append(frame)
        metadata[report.input_provenance.input_id] = asdict(source.metadata)
        for name, run in source.extra.get('candidate_rescoring', {}).items():
            if name in runs and runs[name] != run:
                raise ValueError(f'Inputs disagree on saved prediction run {name!r}')
            runs[name] = run
    combined = reconcile_evidence(TopiaryResult(pd.concat(frames, ignore_index=True), extra={
        'input_results': metadata, 'candidate_rescoring': deepcopy(runs)}))
    added = add_features(combined, models, prefix=prefix, use_flanks=use_flanks)
    updated = []
    for index, report in enumerate(reports):
        frame = added.long_df.loc[added.long_df['_vaxrank_report_index'].eq(index)].copy()
        frame = frame.drop(columns=['_vaxrank_report_index']).reset_index(drop=True)
        source_metadata = deepcopy(report.dataset.result.metadata) if report.dataset is not None else None
        result = TopiaryResult(frame, metadata=source_metadata, form='long')
        for key in ('candidate_rescoring', 'evidence_reconciliation'):
            result.extra[key] = deepcopy(added.extra[key])
        dataset = (replace(report.dataset, result=result) if report.dataset is not None
                   else EpitopeDataset.from_predictions(report.epitopes, result=result))
        scoring = dataset.scoring_frame()
        display = report.report_df.copy()
        features = [column for column in frame if column.startswith(prefix + '__')]
        if features and not display.empty:
            values = scoring[['prediction_id', 'allele', *features]].drop_duplicates()
            values = values.rename(columns={'prediction_id': 'Prediction identity', 'allele': 'Allele'})
            display = display.merge(values, on=['Prediction identity', 'Allele'], how='left', validate='many_to_one')
        display.attrs = {'topiary_df': scoring}
        # Source report adapters keep their recorded construction evidence.
        updated.append(replace(report, report_df=display,
                               dataset=dataset if report.dataset is not None else None))
    added = TopiaryResult(added.long_df.drop(columns=['_vaxrank_report_index']), metadata=added.metadata)
    return updated, added


def prepare_reports(inputs, epitope_config=None, *, mode='input', models=(),
                    alleles=(), input_predictions_path=None, scopes=None,
                    manifest_path=None, genome=None, model_factory=None,
                    prepared_reports=(), use_flanks=True, prediction_prefix=None):
    """Read once and score with fresh common models or source input values."""
    from topiary import is_stated
    if mode not in ('input', 'fresh', 'additive'):
        raise ValueError("Prediction mode must be input, fresh or additive")
    if prediction_prefix is not None and mode != 'additive':
        raise ValueError('A prediction prefix requires additive prediction mode')
    reports = list(prepared_reports) + read_reports(
        inputs, scopes=scopes, manifest_path=manifest_path, genome=genome, alleles=alleles)
    validate_input_scopes([r.input_provenance for r in reports],
                          genome=genome, prediction_alleles=alleles)
    epitope_config = report_config(reports, epitope_config)
    policy_evaluations = []
    original = [e for report in reports for e in report.epitopes]
    identities = [e.prediction_group_key for e in original]
    if len(set(identities)) != len(identities):
        raise ValueError("Inputs contain overlapping candidate identities; do not load an input alongside its native export")
    if not original and not any(r.dataset is not None and not r.dataset.result.df.empty for r in reports):
        raise ValueError("External inputs contain no candidate peptides to score")
    if input_predictions_path:
        Path(input_predictions_path).parent.mkdir(parents=True, exist_ok=True)
        save_predictions(original, input_predictions_path)
    additive_result = None
    if mode == 'additive':
        if not prediction_prefix:
            raise ValueError('Additive predictions require an explicit prediction prefix')
        if model_factory is not None:
            models = model_factory()
        reports, additive_result = _add_predictions(reports, models, prediction_prefix, use_flanks)
    if mode in ('input', 'additive'):
        if models and mode == 'input':
            raise ValueError("Input-table prediction mode does not run fresh models")
        frames_by_report = [r.dataset.scoring_frames() if r.dataset is not None else
                            [r.report_df.attrs['topiary_df']] for r in reports]
        input_frames = [f for frames in frames_by_report for f in frames]
        scoring = pd.concat(input_frames, ignore_index=True) if input_frames else pd.DataFrame()
        scored = []
        for report, frames in zip(reports, frames_by_report):
            for source_frame in frames:
                if source_frame.empty:
                    continue
                if 'kind' in source_frame and not source_frame.kind.map(is_stated).any():
                    continue  # Unmeasured sequence/RNA evidence remains report-only.
                identities = set(source_frame.prediction_id)
                candidates = [e for e in report.epitopes if e.prediction_group_source in identities]
                scored.extend(attach_per_allele_scores(
                    candidates, epitope_config, topiary_df=source_frame,
                    policy_evaluations=policy_evaluations))
    else:
        if any(r.dataset is not None for r in reports):
            raise ValueError(
                "Add predictions to generalized inputs with Topiary's additive rescore_candidates "
                "API before loading; --external-predictions fresh replaces legacy report predictions")
        if model_factory is not None:
            models = model_factory()
        epitopes, fresh_frame = _predict_selected_occurrences(original, models, alleles, use_flanks=use_flanks)
        # Historical metric columns must never look like fresh model values.
        # Keep the common biological evidence vocabulary directly accessible,
        # and retain every original annotation in an explicit input namespace.
        from topiary import EVIDENCE_COLUMNS
        biological = set(EVIDENCE_COLUMNS) | {
            'prediction_id', 'variant', 'gene', 'gene_id', 'transcript',
            'transcript_id', 'antigen_source', 'input_source', 'input_format', 'input_path',
        }
        biological.update(provenance_columns(reports[0].input_provenance))
        records = [replace(record, row={
            **{'input_' + key: value for key, value in record.row.items()},
            **{key: value for key, value in record.row.items() if key in biological},
        }) for report in reports for record in report.records]
        scoring = attach_source_annotations(fresh_frame, records)
        scored = attach_per_allele_scores(epitopes, epitope_config, topiary_df=scoring,
                                          policy_evaluations=policy_evaluations)
    by_id = {e.prediction_group_key: e for e in original if not e.predictions_flat()}
    by_id.update((e.prediction_group_key, e) for e in scored)
    scored = [by_id[e.prediction_group_key] for e in original]
    updated = []
    for report in reports:
        candidates = tuple(by_id[e.prediction_group_key] for e in report.epitopes)
        frame = (_fresh_report_frame(report, candidates) if mode == 'fresh'
                 else report.report_df.copy())
        frame = _annotate_sequence_matches(frame, candidates)
        frame['Prediction evidence'] = mode
        frame = annotate_credited_alleles(frame, candidates)
        frame.attrs = {}
        dataset = report.dataset
        if dataset is not None:
            selection = dict(dataset.selection)
            selection.pop('representatives', None)
            ids = {e.prediction_group_source for e in candidates}
            dataset = replace(dataset, epitopes=candidates, config=epitope_config, selection=selection,
                              policy_evaluations=[evaluation for evaluation in policy_evaluations
                                  if set(evaluation.evidence.df.prediction_id) <= ids])
        updated.append(replace(report, report_df=frame, epitopes=candidates, dataset=dataset))
    frames = [r.report_df for r in updated if not r.report_df.empty]
    frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    frame.attrs['topiary_df'] = scoring
    if mode in ('input', 'additive'):
        frame.attrs['input_score_frames'] = input_frames
    from topiary import TopiaryResult
    # One native dataset carries the scoring frame and the resolved policy.
    # Keep the richer original Topiary result when the input supplies one.
    if mode == 'additive':
        result = additive_result
    elif len(updated) == 1 and updated[0].dataset is not None and mode == 'input':
        result = updated[0].dataset.result
    elif mode == 'input':
        original_frames = []
        for report in reports:
            original_frame = (report.dataset.result.long_df if report.dataset is not None
                              else report.report_df.attrs['topiary_df']).copy()
            original_frame.attrs = {}
            original_frames.append(original_frame)
        result = TopiaryResult(pd.concat(original_frames, ignore_index=True), extra={
            'input_results': {r.input_provenance.input_id: asdict(r.dataset.result.metadata)
                              for r in updated if r.dataset is not None}})
    else:
        result = TopiaryResult(scoring)
    dataset = EpitopeDataset.from_predictions(
        scored, result=result, config=epitope_config, policy_evaluations=policy_evaluations)
    dataset.provenance = tuple(p for r in updated for p in (
        r.dataset.provenance if r.dataset is not None else (r.input_provenance,)))
    dataset.antigens = {key: antigen for r in updated if r.dataset is not None
                        for key, antigen in r.dataset.antigens.items()}
    dataset.mutation_fragments = {key: fragment for r in updated if r.dataset is not None
                                  for key, fragment in r.dataset.mutation_fragments.items()}
    from .native_references import NativeReferences
    dataset.native_references = NativeReferences.combine(
        r.dataset.native_references for r in updated if r.dataset is not None)
    from .native_serialization import to_native_json
    def native_row(row):
        frame = pd.DataFrame([row]).astype(object)
        return frame.where(frame.notna(), None).to_dict('records')[0]

    for report in updated:
        if report.dataset is not None:
            dataset.construction_reports.extend(report.dataset.construction_reports)
            dataset.direct_sources.extend(report.dataset.direct_sources)
        elif report.source_format in ('lens', 'pvacseq'):
            dataset.construction_reports.append({
                'source_format': report.source_format, 'path': report.path,
                'records': [to_native_json(replace(r, row=native_row(r.row)))
                            for r in report.records],
                'rows': to_native_json([native_row(row) for row in report.rows]),
                'genome': to_native_json(resolve_input_genome(genome)) if genome is not None else None,
            })
    dataset.selection = (dict(updated[0].dataset.selection)
                         if len(updated) == 1 and updated[0].dataset is not None else {})
    if mode in ('input', 'additive'):
        dataset.selection['score_groups'] = [list(f.prediction_id.unique()) for f in input_frames]
    frame.attrs['epitope_dataset'] = dataset
    frame.attrs['saved_evidence_input'] = any(r.dataset is not None for r in reports)
    return updated, frame, scored


def load_unified_external(args, epitope_config, options, genome):
    """Connect unified predictions to the existing vaccine/report dispatch."""
    from .epitope_dsl import epitopes_for_ranking
    from .external_input import (
        ExternalInputSummary, lens_ranking_result, pvacseq_ranking_result,
        patient_info_from_external,
    )
    from .ranking import rank_constructs
    from .epitope_dataset import dataset_ranking_result
    genome = resolve_input_genome(genome)
    mode = getattr(args, 'external_predictions', 'input')
    if mode == 'input' and getattr(args, 'external_peptide_only', False):
        raise ValueError('--external-peptide-only requires fresh or additive prediction mode')
    model_factory, alleles = None, []
    direct = bool(getattr(args, 'vcf', None) or getattr(args, 'bam', None))
    cta = bool(getattr(args, 'input_cta_expression', None))
    if direct and mode == 'fresh':
        raise ValueError('Direct/table composition requires --external-predictions input or additive')
    if mode in ('fresh', 'additive'):
        from mhctools.cli import predictors_from_args, mhc_alleles_from_args
        if not getattr(args, 'mhc_predictor', None):
            raise ValueError("Fresh/additive external predictions require --mhc-predictor")
        alleles = sorted({normalize_hla_allele(a) for a in mhc_alleles_from_args(args)})
        def model_factory():
            models = predictors_from_args(args)
            cache_path = getattr(args, 'prediction_cache', None)
            if cache_path:
                from topiary import CachedPredictor
                if not getattr(args, 'external_peptide_only', False):
                    raise ValueError('External prediction caches require --external-peptide-only; '
                                     'contextual exact-cache queries are tracked in Topiary #468')
                if len(models) != 1:
                    raise ValueError('--prediction-cache requires one configured external predictor')
                models = [CachedPredictor.from_topiary_output(cache_path, fallback=models[0])]
            return models
    elif not direct and not cta and getattr(args, 'mhc_predictor', None):
        raise ValueError("Use --external-predictions fresh or additive to run --mhc-predictor")
    elif not direct and not cta and (getattr(args, 'mhc_alleles', None)
                         or getattr(args, 'mhc_alleles_file', None)):
        raise ValueError(
            "Input prediction mode uses the alleles recorded in each report. "
            "Remove --mhc-alleles/--mhc-alleles-file, or use "
            "--external-predictions fresh with --mhc-predictor to predict "
            "for an explicit HLA set.")
    specs = external_input_specs(args)
    loaded_reports = read_reports(
        [(fmt, path) for fmt, path, _ in specs],
        scopes=[scope for _, _, scope in specs],
        manifest_path=getattr(args, 'input_manifest', None), genome=genome)
    if cta:
        from .cta_input import prepare_cta_report
        from .external_input import ExternalConstructOptions
        from .vaccine_config import VaccineConfig
        from .window_selection import WindowSelection
        loaded_reports = [prepare_cta_report(args, genome, epitope_config)]
        vaccine = options.vaccine_config or VaccineConfig()
        if vaccine.window_selection is None:
            import msgspec
            vaccine = msgspec.structs.replace(vaccine, window_selection=WindowSelection())
        options = ExternalConstructOptions.from_configs(
            vaccine_config=vaccine, manufacturability_config=options.manufacturability_config)
    if getattr(args, 'index_native_references', False):
        native = [report.dataset for report in loaded_reports
                  if report.source_format == 'epitopes' and report.dataset is not None]
        if not native:
            raise ValueError('--index-native-references requires a native epitope input')
        for dataset in native:
            dataset.index_native_references()
    from .config.loader import saved_construct_configuration, load_vaxrank_config
    saved_configs = [saved_construct_configuration(r.dataset.selection.get('run_configuration'))
                     for r in loaded_reports if r.dataset is not None]
    saved_configs = [cfg for cfg in saved_configs if cfg is not None]
    if saved_configs:
        # Compare effective merged settings: an explicit current YAML can
        # reconcile different source defaults, but input order cannot choose.
        resolved = [load_vaxrank_config(args, base_config=cfg) for cfg in saved_configs]
        if any(cfg != resolved[0] for cfg in resolved):
            raise ValueError('Saved inputs use different construct configurations; '
                             'supply explicit vaccine configuration overrides')
        args._saved_construct_config = saved_configs[0]
        from .cli.vaccine_config_args import vaccine_config_from_args
        from .external_input import ExternalConstructOptions
        cli_args = getattr(args, '_external_vaccine_args', args)
        vaccine = vaccine_config_from_args(cli_args, merged_config=resolved[0])
        options = ExternalConstructOptions.from_configs(
            vaccine_config=vaccine, manufacturability_config=options.manufacturability_config)
        args.vaccine_peptide_length = vaccine.preferred_peptide_length
        args.num_epitopes_per_vaccine_peptide = vaccine.num_target_epitopes_to_keep
        args.max_vaccine_peptides_per_variant = vaccine.max_vaccine_peptides_per_variant
        args.included_antigen_sources = list(vaccine.included_antigen_sources)
    epitope_config = report_config(loaded_reports, epitope_config)
    if direct:
        from .direct_input import prepare_direct_report
        loaded_reports.insert(0, prepare_direct_report(args, loaded_reports, epitope_config))
    reports, frame, predictions = prepare_reports(
        [], epitope_config, mode=mode, prepared_reports=loaded_reports,
        model_factory=model_factory, alleles=alleles,
        use_flanks=not getattr(args, 'external_peptide_only', False),
        prediction_prefix=getattr(args, 'external_prediction_prefix', None),
        manifest_path=getattr(args, 'input_manifest', None), genome=genome,
        input_predictions_path=getattr(args, 'output_input_predictions', None))
    if cta:
        frame.attrs['cta_predictor'] = loaded_reports[0].report_df.attrs['cta_predictor']
    epitope_config = frame.attrs['epitope_dataset'].config
    dataset = frame.attrs['epitope_dataset']
    if dataset.selection.get('cta_expression_admission'):
        args._cta_assembly = dataset.selection.setdefault('cta_assembly', {})
        args._saved_cta_construct_config = {
            modality: saved['configuration']
            for modality, saved in dataset.selection.get('cta_assembly', {}).items()}
    from .selection_policy import write_run_policy
    duplicate_policy = getattr(args, 'duplicate_candidates', None)
    if duplicate_policy is not None and 'candidate_id' not in dataset.result.df:
        raise ValueError("--duplicate-candidates requires generalized Topiary observations")
    selected = dataset.select_representatives(duplicate_policy)
    run_policy = write_run_policy(args, dataset.policy_evaluations, epitope_config, options.vaccine_config)
    if run_policy is not None:
        dataset.selection['run_configuration'] = run_policy
    generalized_ids = (set(dataset.scoring_frame().loc[
        dataset.result.long_df.source_observation_id.notna(), 'prediction_id'].dropna())
                       if 'source_observation_id' in dataset.result.df else set())
    if selected is not None and not frame.empty:
        frame['Selected observation'] = [
            (identity, allele) in selected if identity in generalized_ids else None
            for identity, allele in
            frame[['Prediction identity', 'Allele']].itertuples(index=False, name=None)]
    rankers = {'lens': lens_ranking_result, 'pvacseq': pvacseq_ranking_result,
               'topiary': dataset_ranking_result, 'epitopes': dataset_ranking_result,
               'vcf_bam': dataset_ranking_result, 'exacto': dataset_ranking_result,
               'cta_expression': dataset_ranking_result}
    outcomes = []
    for report in reports:
        ranker = rankers[report.source_format]
        if report.dataset is not None:
            for key in ('duplicates', 'representatives'):
                if key in dataset.selection:
                    report.dataset.selection[key] = dataset.selection[key]
        ranking_candidates = epitopes_for_ranking(report.epitopes, epitope_config)
        if selected is not None:
            ranking_candidates = [replace(e, per_allele_scores={
                a: score for a, score in e.per_allele_scores.items()
                if e.prediction_group_source not in generalized_ids or (e.prediction_group_source, a) in selected})
                for e in ranking_candidates]
        outcomes.append(ranker(report, ranking_candidates,
                               genome=genome, options=options))
    # A discovery in two tables is not two variants or twice the RNA support.
    # Keep one best source-derived construct per variant; never blend evidence
    # or predictions between contexts. All observations remain in the report.
    observations = defaultdict(list)
    for outcome, report in zip(outcomes, reports):
        for entry in outcome.entries:
            source = entry.ranking_source
            scope = report.input_provenance.scope
            observations[(scope.patient_id, scope.reference_assembly,
                          type(source).__name__, str(source))].append(entry)
    ranked, dna_vaf = [], {}
    for entries in observations.values():
        candidates = [(e.ranking_source, list(e.ranking_peptides)) for e in entries
                      if e.ranking_peptides]
        if options.vaccine_config is not None and options.vaccine_config.window_selection is not None:
            from .window_selection import optimize_antigen_windows
            candidates = [item for candidate in candidates for item in (
                [candidate] if all(p.window_selection_audit for p in candidate[1]) else
                optimize_antigen_windows([candidate], options.vaccine_config))]
        if candidates:
            ranked.append(rank_constructs(candidates)[0])
        values = {e.dna_vaf for e in entries if e.dna_vaf is not None}
        if len(values) == 1 and entries[0].variant is not None:
            dna_vaf[entries[0].variant] = next(iter(values))
    ranked = rank_constructs(ranked)
    passing_path = getattr(args, 'output_passing_variants_csv', None)
    if direct and passing_path:
        from varcode import Variant, StructuralVariant
        from .native_serialization import from_native_json
        selected_variants = {source for source, _ in ranked if isinstance(source, Variant)}
        rows = [{**source['properties'], 'has_vaccine_peptide':
                 from_native_json(source['variant'], (Variant, StructuralVariant)) in selected_variants}
                for source in dataset.direct_sources]
        Path(passing_path).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(passing_path, index=False)
    summary = ExternalInputSummary(
        num_somatic_variants=len(observations),
        num_coding_effect_variants=sum(any(e.resolved_protein_context for e in es)
                                       for es in observations.values()),
        num_variants_with_rna_support=sum(any(e.has_rna_support for e in es)
                                         for es in observations.values()))
    patient = patient_info_from_external(
        ranked, '', getattr(args, 'output_patient_id', '') or '', summary,
        predictions=predictions)
    patient.inputs = [(('VCF/BAM' if r.source_format == 'vcf_bam' else r.source_format + ' report'),
                       r.path) for r in reports]
    patient.input_provenance = [r.input_provenance for r in reports]
    if dataset.selection.get('cta_expression_admission'):
        from .cta_expression import CTAExpressionResult
        from .native_serialization import from_native_json
        admission = from_native_json(dataset.selection['cta_expression_admission'], CTAExpressionResult)
        patient.cta_expression_summary = dict(
            input_features=len(admission.decisions), admitted_targets=len(admission.admitted_antigens),
            selected_targets=len(ranked), measurement_level=admission.input_contract.measurement_level,
            expression_unit=admission.input_contract.expression_unit)
        patient.num_somatic_variants = patient.num_coding_effect_variants = 0
        patient.num_variants_with_rna_support = patient.num_variants_with_vaccine_peptides = 0
    # Preparation already validated every input before scoring. Pooled inputs
    # all declare the same patient/genotype; single-input unknowns stay unknown.
    scope = patient.input_provenance[0].scope
    if scope.patient_id:
        patient.patient_id = patient.patient_id or scope.patient_id
    if scope.mhc_alleles is not None:
        patient.mhc_alleles = list(scope.mhc_alleles)
    elif mode in ('fresh', 'additive'):
        patient.mhc_alleles = alleles
    return ranked, frame, predictions, patient, dna_vaf
