"""Carrier orientation from new exports and conservative historical replay."""

from copy import deepcopy
import json

import pytest

from vaxrank import AntigenDataset, antigens_from_isovar
from vaxrank.isovar_antigens import _variant_frame
from .test_isovar_inputs import attestation, result


def paired_transcripts(left, right, kind):
    producer = result('ncorf')
    first_model, = producer['models']
    first_candidate, = producer['protein_candidates']
    first_model['reference']['resolved_transcript']['strand'] = left
    first_model['reference']['kind'] = kind
    first_model['reference']['model'] = None
    if kind == 'annotated':
        first_model['reference']['orf_type'] = None
    other_model, other_candidate = deepcopy(first_model), deepcopy(first_candidate)
    other_model['model_id'] += ':other-transcript'
    resolved = other_model['reference']['resolved_transcript']
    resolved.update(id=resolved['id'] + ':other', strand=right)
    other_model['reference']['source_transcript_id'] = resolved['id']
    other_candidate.update(candidate_id=other_candidate['candidate_id'] + ':other-transcript',
                           model_id=other_model['model_id'], transcript_id=resolved['id'])
    other_model['candidate_ids'] = [other_candidate['candidate_id']]
    producer['models'].append(other_model)
    producer['protein_candidates'].append(other_candidate)
    return producer


@pytest.mark.parametrize('kind', ['annotated', 'extra_orf'])
@pytest.mark.parametrize('left,right,ambiguous', [
    ('+', '+', False), ('-', '-', False), ('+', '-', True),
    ('+', None, True), (None, None, True),
])
def test_carrier_strands_determine_cross_transcript_frame_ambiguity(tmp_path, kind, left, right, ambiguous):
    producer = paired_transcripts(left, right, kind)
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation())
    assert len(dataset.antigens) == 2
    assert all(a.tumor_specificity.admits_construct == (not ambiguous)
               for a in dataset.antigens.values())
    assert {json.loads(dict(a.source_metadata)['isovar_candidate'])['reference_strand']
            for a in dataset.antigens.values()} == {left, right}
    evidence, = dataset.source_evidence.values()
    assert evidence['isovar_result'] == producer
    path = tmp_path / 'antigens.json'
    dataset.save(path)
    assert AntigenDataset.load(path) == dataset


@pytest.mark.parametrize('resolved,model,expected', [
    ({'strand': '+'}, {'strand': '-'}, '+'),
    ({'strand': '-'}, {'strand': '+'}, '-'),
    ({'strand': None}, {'strand': '+'}, None),
    ({}, {'strand': '+'}, '+'),
    ({}, {'strand': '-'}, '-'),
    ({}, None, None),
])
def test_carrier_is_authoritative_and_legacy_fallback_requires_an_absent_field(resolved, model, expected):
    row = dict(nucleotide_sequence='ATGAAATAA', coding_interval=[0, 9], variant_cdna_interval=[4, 5])
    reference = dict(resolved_transcript=resolved, model=model, translation_table=1)
    assert _variant_frame(row, reference) == (expected, 1, 1)


@pytest.mark.parametrize('invalid', ['', 'forward', 1, True, [], {}])
def test_invalid_carrier_strand_is_not_hidden_by_a_valid_source_model(invalid):
    row = dict(nucleotide_sequence='ATGAAATAA', coding_interval=[0, 9], variant_cdna_interval=[4, 5])
    reference = dict(resolved_transcript={'strand': invalid}, model={'strand': '+'}, translation_table=1)
    with pytest.raises(ValueError, match='strand must be'):
        _variant_frame(row, reference)
