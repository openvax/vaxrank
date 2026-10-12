"""Native source evidence drives shared ranking and explicit source strata."""

from dataclasses import replace
import json

import pandas as pd
import pytest

from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.vaccine_antigen import AminoAcidInterval, TargetableMask

from .test_core_logic_config import _make_epitope
from .test_epitope_dataset import forbidden
from . import test_report_annotation_requirements as fixtures


@pytest.fixture
def mixed_native(tmp_path):
    mutation, _ = fixtures.mutation_native.__wrapped__(tmp_path)
    other = EpitopeDataset.load(fixtures.native_antigens.__wrapped__(tmp_path))
    candidates = list(mutation.epitopes)
    antigens = dict(mutation.antigens)
    for antigen, affinity in zip(other.antigens.values(), [50., 200.]):
        key = antigen.source_identifier
        candidates.append(replace(_make_epitope(
            antigen.amino_acids, ic50=affinity, overlaps_mutation=False), prediction_id=key))
        antigens[key] = antigen
    # Strongest binding alone must not admit a held-out viral candidate.
    held = replace(next(a for a in other.antigens.values() if a.kind == 'viral'),
        amino_acids='SLYNTVATL', source_identifier='held-virus',
        targetable_mask=TargetableMask((AminoAcidInterval(0, 9),)))
    held = replace(held, tumor_specificity=replace(held.tumor_specificity, status='held_out'))
    candidates.append(replace(_make_epitope(held.amino_acids, ic50=1.), prediction_id='held-virus'))
    antigens['held-virus'] = held
    dataset = EpitopeDataset.from_predictions(candidates, antigens=antigens,
        mutation_fragments=mutation.mutation_fragments,
        config=EpitopeConfig(score_expr="1 / affinity['test'].value"))
    path = tmp_path / 'mixed.tsv'
    dataset.save(path)
    return path


@pytest.mark.parametrize('kind,expected', [(None, ['CTA', 'mutation', 'viral']),
                                          ('CTA', ['CTA']), ('viral', ['viral'])])
def test_mixed_native_design_and_source_strata(mixed_native, tmp_path, monkeypatch, kind, expected):
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    output = tmp_path / 'design'
    saved = tmp_path / 'selected.tsv'
    options = ['--input-epitopes', str(mixed_native), '--output-dir', str(output),
               '--output-epitopes', str(saved), '--vaccine-type', 'peptide',
               '--pdf-backend', 'weasyprint',
               '--index-native-references', '--vaccine-peptide-length', '10',
               '--no-processing-aware-annotation',
               '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score']
    if kind:
        options += ['--config-text', f'epitopes.filter_expr=vaxrank_antigen_kind == "{kind}"']
    run_cli(options)
    manifest = json.loads((output / 'manifest.json').read_text())
    sequences = [entry['sequence'] for entry in manifest]
    by_kind = {a.kind: a.amino_acids for a in EpitopeDataset.load(mixed_native).antigens.values()
               if a.tumor_specificity.admits_construct}
    assert len(sequences) == len(expected)
    assert all(sequence in by_kind[k] for k, sequence in zip(expected, sequences))
    assert 'SLYNTVATL' not in sequences
    table = pd.read_csv(output / 'ranked_vaccine_peptides.csv')
    held = table.loc[table['vaxrank_admission_status'] == 'held_out']
    assert len(held) == 1 and 'held out' in held.iloc[0]['Construction limitation']
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(saved), '--output-dir', str(replay),
             '--index-native-references', '--pdf-backend', 'weasyprint',
             '--no-processing-aware-annotation'])
    assert (replay / 'vaccine.fasta').read_bytes() == (output / 'vaccine.fasta').read_bytes()
    assert (replay / 'manifest.json').read_bytes() == (output / 'manifest.json').read_bytes()


def test_native_source_fields_ignore_unverified_producer_labels(tmp_path):
    dataset = EpitopeDataset.load(fixtures.native_antigens.__wrapped__(tmp_path))
    dataset.result.df['vaxrank_antigen_kind'] = 'mutation'
    assert set(dataset.scoring_frame().vaxrank_antigen_kind) == {'CTA', 'viral'}
    dataset.antigens.clear()
    assert dataset.scoring_frame().vaxrank_antigen_kind.isna().all()
    assert set(dataset.result.df.vaxrank_antigen_kind) == {'mutation'}
