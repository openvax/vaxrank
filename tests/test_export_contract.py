"""Table flags and replay bundles have the same meaning across input routes."""

import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pandas as pd
import pytest

from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.peptide import PeptideConstructConfig
from vaxrank.report import make_csv_report, new_columns
from vaxrank.selection_policy import record_construct_configuration

from .test_assembly_antigens import candidate
from .test_epitope_dataset import forbidden

DATA = Path(__file__).parent / 'data/epitope_fixtures'


def test_candidate_table_preserves_sources_and_empty_mutation_fields(tmp_path):
    ranked = [candidate(0), candidate(1, prefix='GGG'), candidate(2, 'GILGFVFTL')]
    csv, xlsx = tmp_path / 'candidates.csv', tmp_path / 'candidates.xlsx'
    make_csv_report(ranked, csv_report_path=csv, excel_report_path=xlsx)
    frame = pd.read_csv(csv)
    pd.testing.assert_frame_equal(frame, pd.read_excel(xlsx, sheet_name='All'), check_dtype=False)
    assert list(frame) == list(new_columns())
    assert list(frame.source_rank) == [1, 2, 3] and set(frame.window_rank) == {1}
    assert frame.chr.isna().all() and frame.expression_score.isna().all()
    assert frame.mutation_start.isna().all()
    assert [json.loads(x) for x in frame.gene_ids] == [['gene-0'], ['gene-1'], ['gene-2']]
    assert list(frame.amino_acids) == [peptides[0].amino_acids for _, peptides in ranked]
    assert frame.source_identifier.nunique() == 3
    assert len(pd.ExcelFile(xlsx).sheet_names) == 4


@pytest.mark.parametrize('flag,filename', [('--input-lens', 'lens_example.tsv'),
                                           ('--input-pvacseq', 'pvacseq_example.tsv')])
def test_cli_rank_tables_have_parity_with_native_reload(tmp_path, monkeypatch, flag, filename):
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    native = tmp_path / 'native.tsv'
    frames = []
    for label, inputs in [('source', [flag, str(DATA / filename)]),
                          ('native', ['--input-epitopes', str(native)])]:
        csv, xlsx, predictions = [tmp_path / (label + suffix) for suffix in ('.csv', '.xlsx', '-predictions.csv')]
        run_cli(inputs + ['--output-csv', str(csv), '--output-xlsx-report', str(xlsx),
                         '--output-neoepitope-report', str(predictions), '--no-processing-aware-annotation']
                + (['--output-epitopes', str(native)] if label == 'source' else []))
        frame = pd.read_csv(csv)
        pd.testing.assert_frame_equal(frame, pd.read_excel(xlsx, sheet_name='All'), check_dtype=False)
        assert list(frame) == list(new_columns()) and not frame.empty
        evidence = pd.read_csv(predictions)
        assert 'Prediction identity' in evidence and 'vaxrank_score' in evidence
        assert 'source_kind' not in evidence
        frames.append(frame)
    pd.testing.assert_frame_equal(*frames)


def test_successful_empty_native_writes_requested_table_schemas(tmp_path, monkeypatch):
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    native = tmp_path / 'empty.tsv'
    EpitopeDataset.from_predictions(()).save(native)
    csv, xlsx, predictions = [tmp_path / n for n in ('rank.csv', 'rank.xlsx', 'predictions.csv')]
    run_cli(['--input-epitopes', str(native), '--output-csv', str(csv),
             '--output-xlsx-report', str(xlsx), '--output-neoepitope-report', str(predictions),
             '--no-processing-aware-annotation'])
    assert pd.read_csv(csv).empty and pd.read_excel(xlsx).empty
    assert list(pd.read_csv(csv)) == list(pd.read_excel(xlsx)) == list(new_columns())
    assert pd.read_csv(predictions).empty and 'vaxrank_score' in pd.read_csv(predictions)


def test_native_bundle_replays_without_optional_policy_exports(tmp_path, monkeypatch):
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.make_pdf_report', lambda **kwargs: None)
    first, second, policy = [tmp_path / p for p in ('first', 'second', 'policy')]
    run_cli(['--input-lens', str(DATA / 'lens_example.tsv'), '--config', 'builtin:lens-v1',
             '--output-dir', str(first), '--vaccine-type', 'peptide', '--output-policy-dir', str(policy),
             '--no-processing-aware-annotation'])
    native = first / 'candidate_predictions.tsv'
    saved = EpitopeDataset.load(native)
    assert saved.policy_evaluations and saved.provenance
    assert (policy / 'selection_policy.json').is_file()
    assert list((policy / 'policy_evidence').glob('*-decisions.tsv'))
    assert not (first / 'input_provenance.json').exists()
    assert not (first / 'selection_policy.json').exists()
    assert not (first / 'ranked_vaccine_peptides.json').exists()
    assert not Path(str(native) + '.policy').exists()
    shutil.rmtree(policy)
    run_cli(['--input-epitopes', str(native), '--output-dir', str(second),
             '--no-processing-aware-annotation'])
    assert (first / 'vaccine.fasta').read_bytes() == (second / 'vaccine.fasta').read_bytes()
    restored = EpitopeDataset.load(second / 'candidate_predictions.tsv')
    assert restored.provenance == saved.provenance and restored.policy_evaluations
    assert restored.selection['construct_assembly'] == saved.selection['construct_assembly']


def test_construct_configuration_never_reads_a_stale_output_file(tmp_path):
    path = tmp_path / 'selection_policy.json'
    content = '{"effective_constructs": {"unrelated": true}}'
    path.write_text(content)
    record_construct_configuration(SimpleNamespace(output_policy_dir=str(tmp_path)), 'peptide',
                                   PeptideConstructConfig())
    assert path.read_text() == content


def test_conflicting_output_paths_fail_before_loading_input(tmp_path):
    path = tmp_path / 'same.csv'
    with pytest.raises(ValueError, match='different output files'):
        run_cli(['--input-epitopes', str(tmp_path / 'missing.tsv'),
                 '--output-csv', str(path), '--output-neoepitope-report', str(path)])
