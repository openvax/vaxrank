"""Whole antigen sequences obey the same hard packing policy in both modalities."""

import pytest

from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs

from .test_assembly_antigens import candidate


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
@pytest.mark.parametrize('content', ['mutation_spanning', 'minimal_epitope'])
@pytest.mark.parametrize('per_construct', [1, 2])
def test_oversized_sequence_does_not_consume_construct_or_type_budget(modality, content, per_construct):
    ranked = [candidate(0, 'A' * 200), candidate(1), candidate(2, 'GILGFVFTL')]
    settings = dict(antigen_content=content, antigens_per_construct=per_construct,
                    max_constructs=1, min_antigen_length_aa=9, max_antigen_length_aa=25,
                    antigen_allocation={'max_per_kind': {'mutation': per_construct}})
    audit = {}
    if modality == 'peptide':
        products = assemble_peptide_constructs(ranked, PeptideConstructConfig(**settings), selection_audit=audit)
    else:
        products = assemble_mrna_constructs(ranked,
            RNAConstructConfig(**settings, max_length_nt=600, optimize_linkers=False,
                               signal_peptide=None, include_mitd=False, utr_3p='HBB'),
            selection_audit=audit, host_taxon_id=9606)
        assert all(len(p.no_polya_nt) <= 600 for p in products)
    assert len(products) == 1
    assert [s['source_identifier'] for a in products[0].antigens for s in a['sources']] == [
        'candidate-%d' % i for i in range(1, per_construct + 1)]
    excluded = audit['antigens'][0]
    assert excluded['outcome'] == 'construct_length_limit'
    assert excluded['length_constraint']['length'] > excluded['length_constraint']['limit']
    assert audit['allocation']['emitted_counts'] == {'mutation': per_construct}


@pytest.mark.parametrize('signal', [None, 'HLA_B'])
@pytest.mark.parametrize('mitd', [False, True])
@pytest.mark.parametrize('starts_with_m', [False, True])
def test_mrna_exact_boundary_includes_every_element(signal, mitd, starts_with_m):
    ranked = [candidate(0, prefix='MA' if starts_with_m else 'AA')]
    settings = dict(signal_peptide=signal, include_mitd=mitd, optimize_linkers=False,
                    max_antigen_length_aa=25, antigens_per_construct=1)
    [baseline] = assemble_mrna_constructs(ranked, RNAConstructConfig(**settings), host_taxon_id=9606)
    length = len(baseline.no_polya_nt)
    [exact] = assemble_mrna_constructs(ranked, RNAConstructConfig(**settings, max_length_nt=length),
                                       host_taxon_id=9606)
    assert exact.no_polya_nt == baseline.no_polya_nt
    audit = {}
    assert assemble_mrna_constructs(ranked, RNAConstructConfig(**settings, max_length_nt=length - 1),
                                    selection_audit=audit, host_taxon_id=9606) == []
    assert audit['antigens'][0]['length_constraint'] == dict(length=length, limit=length - 1, unit='nt')


def test_length_exclusions_report_and_replay_offline(tmp_path, monkeypatch):
    from vaxrank.cli.entry_point import run_cli
    from vaxrank.epitope_dataset import EpitopeDataset
    from vaxrank.native_serialization import from_native_json
    from .test_epitope_dataset import forbidden
    from .test_source_windows_allocation import mixed_source_table
    source = tmp_path / 'mixed.tsv'
    mixed_source_table(source)
    out = tmp_path / 'design'
    run_cli(['--input-topiary', str(source), '--output-dir', str(out),
        '--output-ascii-report', str(out / 'review.txt'),
        '--output-html-report', str(out / 'review.html'),
        '--no-processing-aware-annotation', '--mrna-no-optimize-linkers',
        '--mrna-max-length-nt', '1', '--config-text', 'epitopes.score_expr=1 / affinity.value'])
    native = out / 'candidate_predictions.tsv'
    saved = EpitopeDataset.load(native)
    record = saved.selection['construct_assembly']['mrna']
    assert from_native_json(record['products'], list) == []
    graph = from_native_json(record['source_graph'], dict)
    assert {a['outcome'] for a in graph['assembly_decisions']['antigens']} == {'construct_length_limit'}
    for suffix in ('txt', 'html'):
        text = (out / ('review.' + suffix)).read_text()
        assert 'construct length limit' in text and 'limit 1 nt' in text
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    source.unlink()
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(native), '--output-dir', str(replay),
        '--output-ascii-report', str(replay / 'review.txt'), '--no-processing-aware-annotation'])
    assert 'limit 1 nt' in (replay / 'review.txt').read_text()
    replayed = EpitopeDataset.load(replay / 'candidate_predictions.tsv')
    assert replayed.selection['construct_assembly']['mrna'] == record
