"""Read installed Hitlist observations without acquiring or rebuilding data."""

import hashlib
from importlib.metadata import version
import json

import pandas as pd
from mhcgnomes import parse

from .allele_validation import parse_prediction_allele


def capture_local_hitlist_evidence(sequences, alleles):
    """Capture exact MS queries, or why the local data could not be queried.

    Only current indexes are queried. Unavailable data are not negative
    observations. Captured rows retain Hitlist's source/assay context; this
    does not infer CTA specificity or donor-resolved tissue exclusions.
    """
    from hitlist.dataset import dataset_path, dataset_snapshot
    from hitlist.observations import (
        load_ms_observations, observations_cache_is_current, observations_path,
    )

    peptides = sorted(set(sequences))
    species = sorted({parse(parse_prediction_allele(allele),
                            infer_class2_pairing=False).species.name for allele in alleles if allele})
    evidence = dict(schema_version=1, hitlist_version=version('hitlist'),
                    query=dict(peptides=peptides, mhc_species=species),
                    status='not_installed', species_observations={})
    if not peptides or not species:
        evidence['status'] = 'no_candidates'
        return evidence
    with dataset_snapshot() as directory:
        evidence['data_directory'] = str(directory)
        path = observations_path()
        if not path.exists():
            return evidence
        current = observations_cache_is_current()
        if current is not True:
            evidence['status'] = 'stale' if current is False else 'unverifiable'
            return evidence
        metadata_path = dataset_path('observations_meta.json')
        metadata = metadata_path.read_bytes()
        before = path.stat()
        for name in species:
            frame = load_ms_observations(peptide=peptides, species=name)
            evidence['species_observations'][name] = dict(
                columns=list(frame.columns), rows=json.loads(frame.to_json(orient='records')))
        after = path.stat()
        if (metadata != metadata_path.read_bytes()
                or (before.st_ino, before.st_size, before.st_mtime_ns)
                != (after.st_ino, after.st_size, after.st_mtime_ns)
                or observations_cache_is_current() is not True):
            raise ValueError('Hitlist data changed during evidence capture; retry with a stable dataset')
        evidence.update(status='captured', index_metadata=json.loads(metadata),
                        index_metadata_sha256=hashlib.sha256(metadata).hexdigest())
    return evidence


def hitlist_annotations(sequences, bundle, alleles=None, *, local_evidence=None):
    """Exact public observations; an unqueried sequence has unknown support."""
    columns = ('hitlist_ms_observations', 'hitlist_cta_specific', 'hitlist_avoid_sequence',
               'hitlist_tissue_review_required', 'hitlist_essential_tissue_donors')
    summaries = ({row['peptide']: row for row in bundle['artifacts']['peptides.parquet']['rows']}
                 if bundle else {})
    source_fields = ('n_ms_observations', 'cta_specific', 'avoid_sequence',
                     'tissue_review_required', 'n_essential_tissue_donors')
    from mhcgnomes import ParseError, parse
    from .allele_validation import parse_prediction_allele
    observations = {}
    if bundle:
        for row in bundle['artifacts'].get('presentation.parquet', {}).get('rows', []):
            observations.setdefault(row['peptide'], []).append(row)
    local_observations = {}
    queried_peptides, queried_species, species_by_allele = set(), set(), {}
    if local_evidence and local_evidence['status'] == 'captured':
        queried_peptides = set(local_evidence['query']['peptides'])
        queried_species = set(local_evidence['query']['mhc_species'])
        for species, table in local_evidence['species_observations'].items():
            for record in table['rows']:
                local_observations.setdefault((species, record['peptide']), []).append(record)
    rows = []
    for sequence, allele in zip(sequences, alleles if alleles is not None else [''] * len(sequences)):
        allele = parse_prediction_allele(allele) if allele else ''
        summary = summaries.get(sequence)
        measured = observations.get(sequence, ())
        if not bundle and local_evidence and local_evidence['status'] == 'captured' and allele:
            if allele not in species_by_allele:
                species_by_allele[allele] = parse(allele, infer_class2_pairing=False).species.name
            species = species_by_allele[allele]
            if sequence in queried_peptides and species in queried_species:
                measured = local_observations.get((species, sequence), ())
                summary = dict(n_ms_observations=len(measured))
        row = dict(zip(columns, (summary.get(k) if summary is not None else None
                                 for k in source_fields)))
        row['hitlist_query_status'] = 'captured' if summary is not None else 'unqueried'
        counts = dict(hitlist_ms_monoallelic=0, hitlist_ms_experimental_restriction=0,
                      hitlist_ms_predicted_restriction=0, hitlist_ms_donor_hla_candidate=0,
                      hitlist_ms_coarse_or_untyped=0, hitlist_ms_unknown_restriction_evidence=0)
        for observation in measured:
            try:
                restriction = parse_prediction_allele(observation.get('mhc_restriction', ''))
            except (ParseError, TypeError, ValueError):
                restriction = None  # Coarse/untyped restriction remains in captured raw evidence.
            evidence = observation.get('restriction_evidence')
            counts['hitlist_ms_coarse_or_untyped'] += int(restriction is None)
            counts['hitlist_ms_unknown_restriction_evidence'] += int(evidence not in {'monoallelic', 'experimental', 'predicted'})
            if restriction == allele:
                key = {'monoallelic': 'hitlist_ms_monoallelic',
                       'experimental': 'hitlist_ms_experimental_restriction',
                       'predicted': 'hitlist_ms_predicted_restriction'}.get(evidence)
                if key:
                    counts[key] += 1
            if observation.get('mhc_allele_provenance') in {'sample_allele_match', 'sample_locus_match'}:
                candidates = observation.get('mhc_allele_set') or []
                if isinstance(candidates, str):
                    # Portable bundles encode a JSON list; the local index uses
                    # Hitlist's semicolon-delimited candidate set.
                    candidates = json.loads(candidates) if bundle else candidates.split(';')
                counts['hitlist_ms_donor_hla_candidate'] += int(allele in candidates)
        row.update({k: v if summary is not None else None for k, v in counts.items()})
        rows.append(row)
    return pd.DataFrame(rows)



def annotate_hitlist_frame(frame, evidence):
    """Annotate exact rows in place without changing their index or identities."""
    annotations = hitlist_annotations(
        frame.peptide.tolist(), None, frame.allele.fillna('').tolist(),
        local_evidence=evidence)
    for name in annotations:
        frame[name] = annotations[name].to_numpy()


class LocalHitlistLookup:
    """Capture local queries before scoring and retain them for native export.

    Each query keeps its own peptide/species scope, including unavailable data.
    Previously annotated rows (including unknown results) are frozen on replay.
    """

    def __init__(self):
        self.queries = []

    def __call__(self, frame):
        missing = (frame.hitlist_query_status.isna() if 'hitlist_query_status' in frame
                   else pd.Series(True, index=frame.index))
        if not missing.any():
            return
        rows = frame.loc[missing].copy()
        evidence = capture_local_hitlist_evidence(rows.peptide, rows.allele.fillna(''))
        if self.queries:
            identity = ('status', 'hitlist_version', 'data_directory', 'index_metadata_sha256')
            if any(evidence.get(key) != self.queries[0].get(key) for key in identity):
                raise ValueError('Hitlist data changed between queries; retry with a stable dataset')
        annotate_hitlist_frame(rows, evidence)
        for name in rows:
            if name.startswith('hitlist_'):
                # Assign positionally: Topiary frames need not have a RangeIndex.
                if name not in frame:
                    frame[name] = rows[name].reindex(frame.index)
                else:
                    frame[name] = frame[name].astype(object)
                    frame.loc[missing, name] = rows[name].to_numpy()
        self.queries.append(evidence)

    def annotate_captured(self, frame):
        """Restore annotations on exported predictions without querying Hitlist."""
        species = frame.allele.fillna('').map(
            lambda allele: parse(parse_prediction_allele(allele),
                                 infer_class2_pairing=False).species.name if allele else '')
        for evidence in self.queries:
            query = evidence['query']
            selected = frame.peptide.isin(query['peptides']) & (species.isin(query['mhc_species']) | species.eq(''))
            rows = frame.loc[selected].copy()
            if rows.empty:
                continue
            annotate_hitlist_frame(rows, evidence)
            for name in rows:
                if name.startswith('hitlist_'):
                    if name not in frame:
                        frame[name] = rows[name].reindex(frame.index)
                    else:
                        frame[name] = frame[name].astype(object)
                        frame.loc[selected, name] = rows[name].to_numpy()


def annotate_hitlist_reports(reports):
    """Enrich input rows after occurrence adaptation and preserve source capture."""
    from dataclasses import replace
    from topiary import TopiaryResult

    updated, evidence = [], {}
    for report in reports:
        dataset = report.dataset
        frame = (dataset.result.long_df if dataset is not None
                 else report.report_df.attrs['topiary_df']).copy()
        sources = pd.Series(report.input_provenance.input_id, index=frame.index)
        if dataset is not None:
            by_candidate = {e.prediction_group_source: e.input_provenance.input_id
                            for e in report.epitopes if e.input_provenance is not None}
            sources = dataset.scoring_frame().prediction_id.map(by_candidate).fillna(sources)
        captured = annotate_hitlist_sources(frame, sources)
        if not captured:
            updated.append(report)
            continue
        for identity, queries in captured.items():
            evidence.setdefault(identity, {'hitlist_local_queries': []})['hitlist_local_queries'].extend(queries)
        if dataset is not None:
            dataset = replace(dataset, result=TopiaryResult(frame, metadata=dataset.result.metadata))
            display = dataset.report_frame()
            display.attrs.update(report.report_df.attrs)
            display.attrs['topiary_df'] = dataset.scoring_frame()
        else:
            display = report.report_df.copy()
            columns = [name for name in frame if name.startswith('hitlist_')]
            values = frame[['prediction_id', 'peptide', 'allele', *columns]].drop_duplicates()
            values = values.rename(columns={'prediction_id': 'Prediction identity',
                'peptide': 'Mutant peptide sequence', 'allele': 'Allele'})
            display = display.drop(columns=columns, errors='ignore').merge(values,
                on=['Prediction identity', 'Mutant peptide sequence', 'Allele'],
                how='left', validate='many_to_one')
            display.attrs['topiary_df'] = frame
        updated.append(replace(report, dataset=dataset, report_df=display))
    return updated, evidence


def annotate_hitlist_sources(frame, sources):
    """Capture missing annotations separately for each retained source identity."""
    captured = {}
    for identity in sources.unique():
        selected = sources.eq(identity)
        rows = frame.loc[selected].copy()
        lookup = LocalHitlistLookup()
        lookup(rows)
        if not lookup.queries:
            continue
        for name in rows:
            if name.startswith('hitlist_'):
                if name not in frame:
                    frame[name] = rows[name].reindex(frame.index)
                else:
                    frame[name] = frame[name].astype(object)
                    frame.loc[selected, name] = rows[name].to_numpy()
        captured[identity] = lookup.queries
    return captured
