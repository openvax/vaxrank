"""Read installed Hitlist observations without acquiring or rebuilding data."""

import hashlib
from importlib.metadata import version
import json

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
