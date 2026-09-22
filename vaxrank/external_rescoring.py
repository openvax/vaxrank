"""Combine external reports and select historical or fresh prediction evidence.

The candidate set is the reported peptide occurrences, not all windows in the
source proteins. Source identities and biological admission stay with the
readers; mhctools supplies new predictions and Topiary supplies scoring.
"""
from collections import defaultdict
from dataclasses import replace
import hashlib
from pathlib import Path

import pandas as pd

from .epitope_dsl import attach_per_allele_scores, epitopes_to_topiary_df
from .epitope_io import (
    annotate_credited_alleles, attach_source_annotations, neoepitope_core_row,
    normalize_hla_allele, read_lens_report, read_pvacseq_report, save_predictions,
)
from .external_report import ExternalRecord

READERS = {"lens": read_lens_report, "pvacseq": read_pvacseq_report}


def external_inputs(args):
    """Return ordered, validated (format, path) pairs from CLI arguments."""
    inputs = [(fmt, getattr(args, 'input_' + fmt)) for fmt in READERS
              if getattr(args, 'input_' + fmt, None)]
    for value in getattr(args, 'external_input', None) or ():
        fmt, sep, path = value.partition('=')
        if not sep or fmt not in READERS or not path:
            raise ValueError("--external-input requires lens=PATH or pvacseq=PATH")
        inputs.append((fmt, path))
    return inputs


def _namespace_report(report, input_id):
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
    frame.attrs = {'topiary_df': scoring}
    return replace(
        report, report_df=frame,
        epitopes=tuple(replace(e, prediction_id=ids[e.prediction_id])
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
            key = epitope.prediction_id
            if not key or key in self._entries:
                raise ValueError("Input predictions require unique source identities")
            self._entries[key] = epitope

    def predict_candidates(self, epitopes):
        result = []
        for epitope in epitopes:
            original = self._entries[epitope.prediction_id]
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
    return peptide.sequence, peptide.n_flank or '', peptide.c_flank or ''


def rescore_candidates(epitopes, models, alleles):
    """Re-predict candidate and known comparator sequences using mhctools.

    Each batch has at most one occurrence of any peptide so backends that
    return peptide-grouped results cannot join different flanks accidentally.
    Identical inference contexts across source files are evaluated once.
    """
    if not models or not alleles:
        raise ValueError("Fresh rescoring requires predictors and explicit patient alleles")
    alleles = sorted({normalize_hla_allele(a) for a in alleles})
    contexts = dict.fromkeys(_context(e) for e in epitopes)
    for e in epitopes:
        contexts.update(dict.fromkeys(_context(p) for p in e.comparators.values()
                                      if p.sequence))
    predictions = {key: [] for key in contexts}
    for model in models:
        if not callable(getattr(model, 'predict_with_flanks', None)):
            raise ValueError("%s has no contextual mhctools prediction API" % type(model).__name__)
        support = model.kind_support() if callable(getattr(model, 'kind_support', None)) else {}
        if any(meta.get('mhc_dependence') == 'haplotype' for meta in support.values()):
            raise ValueError("Fresh external rescoring requires per-allele models; "
                             "haplotype-scoped transport is tracked in Topiary #367")
        # Trimming is part of the model's declared inference contract.
        use_flanks = getattr(model, 'uses_flanking_sequences', False)
        default = getattr(model, 'flank_length', 15)
        n_length = getattr(model, 'n_flank_length', None)
        c_length = getattr(model, 'c_flank_length', None)
        n_length = default if n_length is None else n_length
        c_length = default if c_length is None else c_length
        queries = defaultdict(list)
        for key in contexts:
            sequence, n, c = key
            inference = (sequence, n[-n_length:] if use_flanks and n_length else '',
                         c[:c_length] if use_flanks and c_length else '')
            queries[inference].append(key)
        queues = defaultdict(list)
        for query in queries:
            queues[query[0]].append(query)
        while queues:
            batch = [queue.pop(0) for queue in queues.values()]
            queues = {seq: queue for seq, queue in queues.items() if queue}
            by_sequence = {q[0]: q for q in batch}
            results = model.predict_with_flanks(
                [q[0] for q in batch], [q[1] for q in batch], [q[2] for q in batch])
            observed = set()
            for result in results:
                if result.peptide not in by_sequence or result.peptide in observed:
                    raise ValueError("Predictor returned an unexpected or duplicate peptide")
                observed.add(result.peptide)
                query = by_sequence[result.peptide]
                if not result.preds:
                    raise ValueError("Predictor returned no values for %s" % result.peptide)
                for pred in result.preds:
                    if pred.peptide != result.peptide:
                        raise ValueError("Predictor returned mismatched peptide identity")
                    if pred.allele and pred.allele not in alleles:
                        raise ValueError("Predictor returned an allele outside the patient genotype")
                for kind, meta in support.items():
                    if meta.get('mhc_dependence') == 'single_allele':
                        observed_alleles = {p.allele for p in result.preds if p.kind == kind}
                        if observed_alleles != set(getattr(model, 'alleles', ())):
                            raise ValueError("Predictor omitted requested allele coverage for %s" % kind)
                for key in queries[query]:
                    predictions[key].extend(result.preds)
            missing = set(by_sequence) - observed
            if missing:
                raise ValueError("Predictor omitted requested peptides: %s" % sorted(missing))
    result = []
    for epitope in epitopes:
        comparators = {
            name: replace(p, predictions=tuple(predictions[_context(p)]) if p.sequence else ())
            for name, p in epitope.comparators.items()}
        result.append(replace(
            epitope, predictions=tuple(predictions[_context(epitope)]),
            comparators=comparators, patient_alleles=tuple(alleles),
            per_allele_scores={}, allele_attributions=()))
    return result


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
    contexts = {e.prediction_id: e.source_sequence or e.sequence for e in epitopes}
    peptides = {e.prediction_id: e.sequence for e in epitopes}
    if frame.empty:
        return frame
    frame = frame.copy()
    identity = frame['Prediction identity']
    frame['Reported sequence context'] = identity.map(contexts)
    for label, sequences in (('Context', contexts), ('Peptide', peptides)):
        frame[label + ' sequence SHA256'] = identity.map({
            key: hashlib.sha256(sequence.encode('ascii')).hexdigest()
            for key, sequence in sequences.items()})
    frame['Context extent'] = [
        'epitope only' if contexts[key] == peptides[key] else 'reported window'
        for key in identity]
    return frame


def prepare_reports(inputs, epitope_config=None, *, mode='input', models=(),
                    alleles=(), input_predictions_path=None):
    """Read once and score with fresh common models or source input values."""
    if mode not in ('input', 'fresh'):
        raise ValueError("Prediction mode must be input or fresh")
    reports, seen = [], set()
    for fmt, path in inputs:
        if fmt not in READERS:
            raise ValueError("Unsupported external format: %s" % fmt)
        input_id = fmt + ':' + hashlib.sha256(Path(path).read_bytes()).hexdigest()
        if input_id in seen:
            raise ValueError("The same external input was supplied more than once: %s" % path)
        seen.add(input_id)
        report = READERS[fmt](path, score_epitopes=False)
        reports.append(_namespace_report(report, input_id))
    original = [e for report in reports for e in report.epitopes]
    if not original:
        raise ValueError("External inputs contain no candidate peptides to score")
    if input_predictions_path:
        Path(input_predictions_path).parent.mkdir(parents=True, exist_ok=True)
        save_predictions(original, input_predictions_path)
    if mode == 'input':
        if models:
            raise ValueError("Input-table prediction mode does not run fresh models")
        epitopes = InputTablePredictor(original).predict_candidates(original)
        input_frames = [r.report_df.attrs['topiary_df'] for r in reports]
        scoring = pd.concat(input_frames, ignore_index=True) if input_frames else pd.DataFrame()
        cached = {e.prediction_id: e for e in epitopes}
        scored = [e for report, source_frame in zip(reports, input_frames)
                  for e in attach_per_allele_scores(
                      [cached[e.prediction_id] for e in report.epitopes],
                      epitope_config, topiary_df=source_frame)]
    else:
        epitopes = rescore_candidates(original, models, alleles)
        # Historical metric columns must never look like fresh model values.
        # Keep the common biological evidence vocabulary directly accessible,
        # and retain every original annotation in an explicit input namespace.
        from topiary import EVIDENCE_COLUMNS
        biological = set(EVIDENCE_COLUMNS) | {
            'prediction_id', 'variant', 'gene', 'gene_id', 'transcript',
            'transcript_id', 'antigen_source', 'input_source', 'input_format', 'input_path',
        }
        records = [replace(record, row={
            **{'input_' + key: value for key, value in record.row.items()},
            **{key: value for key, value in record.row.items() if key in biological},
        }) for report in reports for record in report.records]
        scoring = attach_source_annotations(epitopes_to_topiary_df(epitopes), records)
        scored = attach_per_allele_scores(epitopes, epitope_config, topiary_df=scoring)
    by_id = {e.prediction_id: e for e in scored}
    updated = []
    for report in reports:
        candidates = tuple(by_id[e.prediction_id] for e in report.epitopes)
        frame = (_fresh_report_frame(report, candidates) if mode == 'fresh'
                 else report.report_df.copy())
        frame = _annotate_sequence_matches(frame, candidates)
        frame['Prediction evidence'] = mode
        frame = annotate_credited_alleles(frame, candidates)
        frame.attrs = {}
        updated.append(replace(report, report_df=frame, epitopes=candidates))
    frames = [r.report_df for r in updated if not r.report_df.empty]
    frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    frame.attrs['topiary_df'] = scoring
    if mode == 'input':
        frame.attrs['input_score_frames'] = input_frames
    return updated, frame, scored


def load_unified_external(args, epitope_config, options, genome):
    """Connect unified predictions to the existing vaccine/report dispatch."""
    from .epitope_dsl import epitopes_for_ranking
    from .external_input import (
        ExternalInputSummary, lens_ranking_result, pvacseq_ranking_result,
        patient_info_from_external, ranked_sorted_by_target_score,
    )
    mode = getattr(args, 'external_predictions', 'input')
    models, alleles = [], []
    if mode == 'fresh':
        from mhctools.cli import predictors_from_args, mhc_alleles_from_args
        if not getattr(args, 'mhc_predictor', None):
            raise ValueError("Fresh external predictions require --mhc-predictor")
        alleles = sorted({normalize_hla_allele(a) for a in mhc_alleles_from_args(args)})
        models = predictors_from_args(args)
    elif getattr(args, 'mhc_predictor', None):
        raise ValueError("Use --external-predictions fresh to run --mhc-predictor")
    reports, frame, predictions = prepare_reports(
        external_inputs(args), epitope_config, mode=mode, models=models,
        alleles=alleles,
        input_predictions_path=getattr(args, 'output_input_predictions', None))
    outcomes = []
    for report in reports:
        ranker = lens_ranking_result if report.source_format == 'lens' else pvacseq_ranking_result
        outcomes.append(ranker(report, epitopes_for_ranking(report.epitopes, epitope_config),
                               genome=genome, options=options))
    # A discovery in two tables is not two variants or twice the RNA support.
    # Keep one best source-derived construct per variant; never blend evidence
    # or predictions between contexts. All observations remain in the report.
    observations = defaultdict(list)
    for outcome in outcomes:
        for entry in outcome.entries:
            source = entry.ranking_source
            observations[(type(source).__name__, str(source))].append(entry)
    ranked, dna_vaf = [], {}
    for entries in observations.values():
        candidates = [(e.ranking_source, [e.vaccine_peptide]) for e in entries
                      if e.vaccine_peptide is not None]
        if candidates:
            ranked.append(ranked_sorted_by_target_score(candidates)[0])
        values = {e.dna_vaf for e in entries if e.dna_vaf is not None}
        if len(values) == 1 and entries[0].variant is not None:
            dna_vaf[entries[0].variant] = next(iter(values))
    ranked = ranked_sorted_by_target_score(ranked)
    summary = ExternalInputSummary(
        num_somatic_variants=len(observations),
        num_coding_effect_variants=sum(any(e.resolved_protein_context for e in es)
                                       for es in observations.values()),
        num_variants_with_rna_support=sum(any(e.has_rna_support for e in es)
                                         for es in observations.values()))
    patient = patient_info_from_external(
        ranked, '', getattr(args, 'output_patient_id', '') or '', summary,
        predictions=predictions)
    patient.inputs = [(r.source_format + ' report', r.path) for r in reports]
    if mode == 'fresh':
        patient.mhc_alleles = alleles
    return ranked, frame, predictions, patient, dna_vaf
