"""Native reference batching must equal exhaustive per-antigen matching."""

from dataclasses import replace
from types import SimpleNamespace

import pandas as pd
import pytest
from topiary import TopiaryPredictor

from vaxrank.antigen_input import AntigenDataset, _batched_self_matches, predict_antigens
from vaxrank.input_scope import InputScope
from vaxrank.reference_proteome import self_reference_matches
from vaxrank.vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation, VaccineAntigen,
)

from .test_external_rescoring import ALLELES, ContextPredictor
from .test_reference_proteome import create_mock_transcript


def reference_inputs():
    antigen = VaccineAntigen(
        kind='CTA', amino_acids='ACDEFGHIKL', gene_id='GCTA',
        source_identifier='source-a', species='Homo sapiens',
        self_reference_excluded_gene_ids=('GCTA',),
        targetable_mask=TargetableMask((AminoAcidInterval(0, 10),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture', evidence_source=__name__))
    records = {
        'cta-a': antigen,
        'cta-b': replace(antigen, source_identifier='source-b'),
        'cta-unexcluded': replace(antigen, self_reference_excluded_gene_ids=()),
        'viral': replace(antigen, kind='viral', self_reference_excluded_gene_ids=()),
        'held-viral': replace(antigen, kind='viral', self_reference_excluded_gene_ids=(),
                              tumor_specificity=replace(antigen.tumor_specificity, status='held_out')),
    }
    inputs = AntigenDataset(records, InputScope(
        patient_id='test', sample_id='sample', reference_assembly='GRCh38',
        annotation='ensembl:93', mhc_alleles=tuple(ALLELES)))
    transcripts = [
        create_mock_transcript('CTA', 'ACDEFGHIKL', gene_id='GCTA'),
        create_mock_transcript('SELF1', 'MACDEFGHIKV', gene_id='GSELF'),
        create_mock_transcript('SELF2', 'MACDEFGHIKV', gene_id='GSELF'),
    ]
    genome = SimpleNamespace(transcripts=lambda: transcripts, release=93,
                             reference_name='GRCh38',
                             species=SimpleNamespace(latin_name='Homo sapiens'))
    return inputs, genome


def test_batched_native_prediction_equals_exhaustive_sources_and_measurements(monkeypatch):
    inputs, genome = reference_inputs()
    calls = []

    def tracked(peptides, antigen, reference):
        calls.append((antigen.kind, antigen.self_reference_excluded_gene_ids, tuple(peptides)))
        return self_reference_matches(peptides, antigen, reference)

    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches', tracked)
    batched = inputs.predict(TopiaryPredictor(models=[ContextPredictor()]), genome=genome)
    assert len(calls) == 3  # CTA excluded, CTA whole host, viral whole host.
    assert all(len(peptides) == 2 for _, _, peptides in calls)
    exhaustive_calls = []

    def exhaustive(peptides, antigen, reference):
        exhaustive_calls.append(antigen.source_identifier)
        return self_reference_matches(peptides, antigen, reference)

    exhaustive_result = predict_antigens(
        inputs.antigens, TopiaryPredictor(models=[ContextPredictor()]), genome=genome,
        provenance=inputs.prediction_provenance(genome=genome), self_matcher=exhaustive)
    assert len(exhaustive_calls) == 7  # Custom matchers still see each source/policy.
    assert batched.epitopes == exhaustive_result.epitopes
    assert batched.antigens == exhaustive_result.antigens
    assert batched.selection == exhaustive_result.selection
    pd.testing.assert_frame_equal(batched.result.df, exhaustive_result.result.df)
    shared = next(e for e in batched.epitopes
                  if e.source_name == 'cta-a' and e.sequence == 'ACDEFGHIK')
    assert [s.transcript_id for s in shared.self_reference_match.sources] == ['SELF1', 'SELF2']
    assert all(s.gene_id == 'GSELF' for s in shared.self_reference_match.sources)
    cta_only = next(e for e in batched.epitopes
                   if e.source_name == 'cta-a' and e.sequence == 'CDEFGHIKL')
    assert cta_only.occurs_in_reference and not cta_only.occurs_in_non_CTA_reference
    assert len({e.prediction_id for e in batched.epitopes}) == 10


@pytest.mark.parametrize('genome_available', [False, True])
def test_repeated_windows_share_queries_without_losing_source_occurrences(monkeypatch, genome_available):
    inputs, genome = reference_inputs()
    genome = genome if genome_available else None
    antigen = inputs.antigens['viral']
    antigens = {f'window-{i}': replace(antigen, source_identifier=f'window-{i}') for i in range(300)}
    frame = pd.DataFrame([dict(source_sequence_name=name, peptide=p)
                          for name in antigens for p in ['ACDEFGHIK', 'CDEFGHIKL', 'ACDEFGHIK']])
    frame.source_sequence_name = pd.Categorical(frame.source_sequence_name,
                                                categories=[*antigens, 'unused-category'])
    calls = []

    def tracked(peptides, source, reference):
        calls.append(tuple(peptides))
        return self_reference_matches(peptides, source, reference)

    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches', tracked)
    matches, host = _batched_self_matches(frame, antigens, genome)
    assert calls == [('ACDEFGHIK', 'CDEFGHIKL')]
    assert set(matches) == set(host) == set(antigens)
    expected = self_reference_matches(calls[0], antigen, genome)
    assert all(matches[name] == host[name] == expected for name in antigens)
    assert all(m.source_provenance_complete == genome_available for m in expected.values())
