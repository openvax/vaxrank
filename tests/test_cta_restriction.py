"""Synthetic RNA bounds test a declared research policy, never a real CTA panel."""

from dataclasses import asdict, replace
import json
from types import SimpleNamespace

import pytest

from vaxrank.cta_restriction import (
    CTARestrictionEvidence, CTARestrictionObservation, CTARestrictionOccurrence,
    CTARestrictionPolicy, assess_restriction,
)


def policy():
    return CTARestrictionPolicy('synthetic-test', '1', 'TPM', 'source_native',
        1, 5, ('heart', 'kidney'), ('testis',), 2, 2, 'Synthetic fixture only')


def occurrence(gene='1', transcript='NM_1.1'):
    return CTARestrictionOccurrence(gene, transcript, 'NP_1.1', 'a' * 64)


def snapshot(taxon=9615, facts=()):
    return SimpleNamespace(taxon_id=taxon, snapshot_id='openvax:cta-evidence:v2:sha256:' + 'b' * 64,
        candidates=lambda: [dict(id='candidate', gene_ids=['1'])],
        as_dict=lambda: dict(facts=list(facts), evidence=[dict(id='native', candidate_ids=['candidate'])]))


def evidence(taxon=9615, targets=None):
    rows = []
    for target in targets or (occurrence(),):
        for tissue in ('heart', 'kidney', 'testis'):
            for donor in ('donor1', 'donor2'):
                key = f'{target.gene_id}:{target.transcript_id}:{target.entry_id}:{tissue}:{donor}'
                rows.append(CTARestrictionObservation(key, (target,), tissue,
                    'reproductive' if tissue == 'testis' else 'somatic', 'transcript',
                    'TPM', 'source_native', 6 if tissue == 'testis' else 0,
                    8 if tissue == 'testis' else 0.5, 'independent',
                    'Synthetic allocated transcript bounds', f'{tissue}:{donor}',
                    (donor,), False, 'synthetic-study', '1', 'c' * 64, key,
                    ('synthetic-study',)))
    return CTARestrictionEvidence(taxon, 'openvax:reference:v1:sha256:' + 'd' * 64,
        'e' * 64, snapshot(taxon).snapshot_id, tuple(rows))


def codes(result):
    return {f.code for f in result.findings}


@pytest.mark.parametrize('taxon', [9606, 9615, 10090])
def test_equal_rna_evidence_has_equal_policy_result_for_every_species(taxon):
    result = assess_restriction(policy(), evidence(taxon), (occurrence(),), snapshot(taxon))
    assert result.passed


def test_missing_rna_panel_and_donors_are_unknown_not_zero():
    missing = assess_restriction(policy(), None, (occurrence(),), snapshot())
    assert codes(missing) == {'missing_restriction_evidence', 'insufficient_independent_donors',
                             'insufficient_independent_samples'}
    original = evidence()
    rows = tuple(r for r in original.observations if r.tissue != 'kidney')
    result = assess_restriction(policy(), replace(original, observations=rows), (occurrence(),), snapshot())
    assert not result.passed
    assert {f.tissue for f in result.findings} == {'kidney'}


@pytest.mark.parametrize('change', [dict(pooled=True, donor_ids=('a', 'b')),
    dict(donor_ids=()), dict(allocation='unknown'), dict(feature_level='gene'),
    dict(exclusion_reason='quality exclusion'), dict(rna_unit='counts')])
def test_unresolved_or_pooled_measurements_cannot_supply_independent_coverage(change):
    original = evidence()
    changed = replace(original, observations=tuple(replace(r, **change) for r in original.observations))
    result = assess_restriction(policy(), changed, (occurrence(),), snapshot())
    assert not result.passed
    assert 'insufficient_independent_donors' in codes(result)


def test_repeated_assays_and_sources_do_not_inflate_donor_or_sample_counts():
    original = evidence()
    rows = []
    for row in original.observations:
        if row.donor_ids == ('donor1',):
            rows.extend((row, replace(row, evidence_id=row.evidence_id + ':repeat',
                source_id='same-data-other-paper', lineage_ids=('synthetic-study', 'secondary-paper'))))
    result = assess_restriction(policy(), replace(original, observations=tuple(rows)),
                                (occurrence(),), snapshot())
    assert codes(result) == {'insufficient_independent_donors', 'insufficient_independent_samples'}


def test_all_source_heart_counterexample_survives_exclusion_and_good_target_panel():
    other = occurrence('2', 'NM_2.1')
    original = evidence(targets=(occurrence(), other))
    rows = tuple(replace(r, lower_bound=20, upper_bound=22, exclusion_reason='caller prefers another source')
        if r.occurrences == (other,) and r.tissue == 'heart' else r for r in original.observations)
    result = assess_restriction(policy(), replace(original, observations=rows),
                                (occurrence(), other), snapshot())
    assert not result.passed
    assert any(f.code == 'somatic_rna_above_bound' and f.subject.startswith('2:') for f in result.findings)


@pytest.mark.parametrize('bounds', [(None, None), (0, None), (0, 10)])
def test_uncertain_somatic_bounds_survive_other_complete_panel_rows(bounds):
    original = evidence()
    unknown = replace(original.observations[0], evidence_id='uncertain',
                      lower_bound=bounds[0], upper_bound=bounds[1])
    result = assess_restriction(policy(), replace(original, observations=original.observations + (unknown,)),
                                (occurrence(),), snapshot())
    assert 'unresolved_somatic_rna' in codes(result)


def test_native_somatic_facts_survive_supplemental_negative_panel():
    fact = dict(id='native:qualitative', evidence_id='native', tissue='heart', tissue_role='somatic',
                observation=dict(status='qualitative'))
    result = assess_restriction(policy(), evidence(), (occurrence(),), snapshot(facts=(fact,)))
    assert codes(result) == {'somatic_fact_requires_review'}
    assert result.findings[0].evidence_id == 'native'


def test_reference_bundle_json_round_trip_and_closed_fields(tmp_path):
    original = evidence()
    path = tmp_path / 'evidence.json'
    path.write_text(json.dumps(asdict(original)))
    assert CTARestrictionEvidence.read(path) == original
    data = asdict(original)
    data['treat_unknown_as_zero'] = True
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='Unknown'):
        CTARestrictionEvidence.read(path)


def test_conflicting_sample_identity_and_panel_roles_are_rejected():
    original = evidence()
    row = replace(original.observations[0], evidence_id='conflict', donor_ids=('another',))
    with pytest.raises(ValueError, match='Conflicting CTA sample'):
        replace(original, observations=original.observations + (row,))
    rows = tuple(replace(r, tissue_role='reproductive') if r.tissue == 'heart' else r
                 for r in original.observations)
    with pytest.raises(ValueError, match='tissue role'):
        assess_restriction(policy(), replace(original, observations=rows), (occurrence(),), snapshot())
