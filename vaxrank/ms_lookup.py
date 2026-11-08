"""Capture source-verified Hitlist evidence for the MS scoring profile."""

from functools import lru_cache
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import pandas as pd

from .ms_evidence import _texts


RECORDS_COLUMN = 'hitlist_ms_records_json'


def normalize_observation(row):
    """Normalize a Hitlist positive-MS row without inventing source identity."""
    from .cells import text
    reasons = []
    if row.get('has_ptm') is not False:
        reasons.append('modified_or_unknown_chemical_form')
    if row.get('is_non_peptide_ligand') is True:
        reasons.append('non_peptide_ligand')
    if row.get('is_potential_contaminant') is True:
        reasons.append('potential_contaminant')
    references = []
    pmid = text(row.get('pmid'))
    if pmid:
        try:
            pmid = str(int(float(pmid)))
        except ValueError:
            raise ValueError(f'Invalid Hitlist PMID: {pmid!r}') from None
        references.append('pmid:' + pmid)
    if text(row.get('reference_iri')):
        references.append(text(row['reference_iri']))
    identity = text(row.get('provenance_id')) or text(row.get('assay_iri'))
    if not identity:
        # An unresolved source row can support presence; its content digest is
        # an audit locator, never a donor, specimen or independent observation.
        identity = 'unresolved:' + hashlib.sha256(
            json.dumps(row, sort_keys=True, allow_nan=False).encode()).hexdigest()
    typing = (row.get('mhc_allele_set') if row.get('mhc_allele_provenance')
              in ('sample_allele_match', 'sample_locus_match') else None)
    mhc_class = text(row.get('mhc_class')) or None
    restriction = text(row.get('mhc_restriction'))
    basis = row.get('restriction_evidence')
    # Unresolved class supports sequence transfer only.
    if mhc_class not in ('I', 'II'):
        typing, restriction, basis = None, None, 'unknown'
    result = dict(observation_id=identity, peptide=row['peptide'],
        eligible=not reasons, exclusion_reasons=reasons,
        mhc_species=text(row.get('mhc_species')), mhc_class=mhc_class,
        restriction=restriction, restriction_basis=basis,
        monoallelic_verified=(row.get('is_monoallelic') is True and basis == 'monoallelic'),
        typed_alleles=list(_texts(typing)), references=references,
        source_prediction_qualified=False, quality=1., context_weight=1.,
        source_dataset=text(row.get('source')), source_assay_id=text(row.get('assay_iri')),
        chemical_form_status='unmodified' if row.get('has_ptm') is False else 'modified_or_unknown')
    # Only explicit reviewed axes qualify. The flat-reader contract normally
    # has none; labels and descriptive count_samples fallbacks are not IDs.
    engineered = row.get('is_engineered_mhc') is True
    axis = 'experimental_origin' if engineered else 'donor'
    units = _texts(row.get(axis + '_ids'))
    if row.get(axis + '_status') == 'resolved' and len(units) == 1:
        result.update(unit_id=units[0],
                      unit_kind='engineered_origin' if engineered else 'native_donor')
    return result


@lru_cache(maxsize=4)
def verified_package(directory, manifest_digest):
    """Verify a complete producer package once per manifest identity."""
    from hitlist.mhc_evidence import verify_mhc_evidence_package
    path = Path(directory)
    verified = verify_mhc_evidence_package(path)
    if hashlib.sha256((path / 'manifest.json').read_bytes()).hexdigest() != manifest_digest:
        raise ValueError('Hitlist package changed during verification')
    return verified


def capture_package(directory, sequences, alleles):
    """Query the original positive partition of a verified MHC package."""
    import pyarrow.dataset as ds
    from mhcgnomes import parse
    from .allele_validation import parse_prediction_allele

    directory = Path(directory).resolve()
    manifest_bytes = (directory / 'manifest.json').read_bytes()
    digest = hashlib.sha256(manifest_bytes).hexdigest()
    manifest = verified_package(str(directory), digest)
    path = directory / 'native' / 'observations.parquet'
    before = path.stat()
    species = sorted({parse(parse_prediction_allele(a), infer_class2_pairing=False).species.name
                      for a in alleles if a})
    peptides = sorted(set(sequences))
    frame = ds.dataset(path, format='parquet').to_table(filter=(
        ds.field('peptide').isin(peptides) & ds.field('mhc_species').isin(species))).to_pandas()
    # Current package verification establishes the positive-ligand partition's
    # modality/outcome contract. Retain every field and the package identity.
    after = path.stat()
    if (before.st_ino, before.st_mtime_ns, before.st_size) != (
            after.st_ino, after.st_mtime_ns, after.st_size) or manifest_bytes != (
                directory / 'manifest.json').read_bytes():
        raise ValueError('Hitlist package changed during capture')
    return dict(schema_version=1, hitlist_version=version('hitlist'), status='captured',
        query=dict(peptides=peptides, mhc_species=species), data_directory=str(directory),
        index_metadata_sha256=digest, index_metadata=manifest,
        source_contract='verified_hitlist_mhc_package', species_observations={
            species_name: dict(columns=list(frame), rows=json.loads(
                frame.loc[frame.mhc_species.eq(species_name)].to_json(orient='records')))
            for species_name in species})


class MSLookup:
    """One profile lookup with portable source captures and normalized rows."""

    def __init__(self, config):
        from .ms_evidence import validate_config
        self.config = validate_config(config)
        self.queries = []

    def __call__(self, frame):
        from .hitlist_evidence import capture_local_hitlist_evidence, annotate_hitlist_frame
        missing = (frame[RECORDS_COLUMN].isna() if RECORDS_COLUMN in frame
                   else pd.Series(True, index=frame.index))
        if not missing.any():
            return
        rows = frame.loc[missing].copy()
        if self.config.get('lookup', 'auto') == 'off':
            evidence = dict(status='disabled', query=dict(peptides=sorted(set(rows.peptide)),
                mhc_species=[]), species_observations={})
        elif self.config.get('package'):
            evidence = capture_package(self.config['package'], rows.peptide, rows.allele.fillna(''))
        else:
            evidence = capture_local_hitlist_evidence(rows.peptide, rows.allele.fillna(''))
        if self.queries:
            previous = self.queries[0]
            identity = ('status', 'hitlist_version', 'data_directory', 'index_metadata_sha256')
            if any(evidence.get(key) != previous.get(key) for key in identity):
                raise ValueError('Hitlist data changed between MS queries')
        if self.config.get('lookup') == 'required' and evidence['status'] != 'captured':
            raise ValueError('Required MS evidence is ' + evidence['status'])
        if evidence['status'] != 'disabled':
            annotate_hitlist_frame(rows, evidence)
        self._annotate_records(rows, evidence)
        for name in rows:
            if name.startswith('hitlist_'):
                if name not in frame:
                    frame[name] = rows[name].reindex(frame.index)
                else:
                    frame[name] = frame[name].astype(object)
                    frame.loc[missing, name] = rows[name].to_numpy()
        self.queries.append(evidence)

    @staticmethod
    def _annotate_records(frame, evidence):
        records = {}
        for table in evidence.get('species_observations', {}).values():
            for row in table['rows']:
                records.setdefault(row['peptide'], []).append(normalize_observation(row))
        frame[RECORDS_COLUMN] = [json.dumps(records.get(p, []), sort_keys=True, allow_nan=False)
                                 for p in frame.peptide]
        frame['hitlist_ms_availability'] = evidence['status']

    def annotate_captured(self, frame):
        for evidence in self.queries:
            selected = frame.peptide.isin(evidence['query']['peptides'])
            rows = frame.loc[selected].copy()
            self._annotate_records(rows, evidence)
            for name in (RECORDS_COLUMN, 'hitlist_ms_availability'):
                if name not in frame:
                    frame[name] = rows[name].reindex(frame.index)
                else:
                    frame[name] = frame[name].astype(object)
                    frame.loc[selected, name] = rows[name].to_numpy()


def complete_candidate_frame(frame, sequences, predictor, *, alleles=None):
    """Enumerate declared native windows even when a model emits no result.

    Missing rows carry no prediction kind, method or value. They are candidate
    identities, not synthetic predictions. Existing predictor rows are retained.
    """
    from .allele_validation import parse_prediction_allele
    models = getattr(predictor, 'models', (predictor,))
    requests = set()
    for model in models:
        model_alleles = getattr(model, 'alleles', ()) or alleles or ()
        lengths = getattr(model, 'default_peptide_lengths', ())
        if model_alleles and not lengths:
            raise ValueError('MS candidate enumeration requires declared model peptide lengths')
        for allele in model_alleles:
            for length in lengths:
                if type(length) is not int or length < 1:
                    raise ValueError('MS candidate peptide lengths must be positive integers')
                requests.add((parse_prediction_allele(allele), length))
    if not requests:
        raise ValueError('MS native prediction requires declared patient alleles and peptide lengths')
    keys = ['source_sequence_name', 'peptide_offset', 'peptide', 'allele']
    existing = set(frame[keys].itertuples(index=False, name=None)) if not frame.empty else set()
    missing = []
    for name, sequence in sequences.items():
        for allele, length in sorted(requests):
            for offset in range(len(sequence) - length + 1):
                peptide = sequence[offset:offset + length]
                if (name, offset, peptide, allele) not in existing:
                    missing.append(dict(source_sequence_name=name, peptide_offset=offset,
                        peptide=peptide, peptide_length=length, allele=allele,
                        kind=None, prediction_method_name=None, predictor_version=None,
                        value=None, score=None, percentile_rank=None))
    result = pd.concat([frame, pd.DataFrame(missing)], ignore_index=True) if missing else frame.copy()
    result.attrs = dict(frame.attrs)
    return result
