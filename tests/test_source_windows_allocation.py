"""Distinct biological windows and allocation budgets must survive both writers."""

from dataclasses import replace

import pytest

from vaxrank.antigen_allocation import AntigenAllocation
from vaxrank.cli.arg_parser import external_input_arg_parser
from vaxrank.cli.entry_point import _peptide_construct_options
from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
from vaxrank.config.loader import load_vaxrank_config, extract_vaccine_config_kwargs
from vaxrank.epitope_logic import slice_epitopes
from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.vaccine_config import VaccineConfig
from vaxrank.window_selection import select_source_windows, optimize_antigen_windows, target_content

from .test_assembly_antigens import candidate
from .test_window_selection import vaccine


def windows():
    full = vaccine('AASIINFEKLAA' + 'AAGILGFVFTLA' + 'AANLVPMVATVA' + 'AAGLLGATCMFA' + 'AASIINFEKLAA',
                   [(2, 10, 10., False), (14, 23, 3., False), (26, 35, 2., False),
                    (38, 47, 1., False), (50, 58, 10., False)])
    full = replace(full, ranking_rules=('target_epitope_score',))
    spans = [(0, 12), (1, 13), (12, 24), (24, 36), (36, 48), (48, 60)]
    return full, [replace(full, antigen=full.antigen.sliced(start, end),
        epitopes=slice_epitopes(full.epitopes, start, end),
        window_selection_audit={'parent_window': dict(sequence=full.amino_acids, start=start, end=end)})
        for start, end in spans]


def test_three_regions_reject_overlap_and_repeated_content_but_allow_touching():
    full, choices = windows()
    selected = select_source_windows(choices)
    assert [p.window_selection_audit['parent_window']['start'] for p in selected] == [0, 12, 24]
    # The second region is 30% of the first score: the near-best cutoff resets.
    assert [p.target_epitope_score for p in selected] == [10., 3., 2.]
    reasons = [a['decision_reason'] for a in selected[0].window_selection_audit['alternatives']]
    assert reasons == ['selected_source_window', 'overlaps_selected_window', 'selected_source_window',
                       'selected_source_window', 'source_window_limit', 'no_additional_target_content']
    alternatives = selected[0].window_selection_audit['manufacturing_alternatives']
    # The shifted first window would overlap the second; it is not a valid alternative.
    assert alternatives == []
    _, choices = windows()
    first = select_source_windows(choices, limit=1)[0]
    assert first.window_selection_audit['manufacturing_alternatives'][0]['parent_window']['start'] == 1
    assert target_content(first, positioned=True) == {('SIINFEKL', 'HLA-A*02:01', 2)}


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_all_biological_windows_reach_assembly_with_one_alternative_per_slot(modality):
    full, choices = windows()
    selected = select_source_windows(choices)
    options = dict(min_antigen_length_aa=8, max_antigen_length_aa=12, max_constructs=3,
                   antigens_per_construct=1, candidates_per_slot=1)
    if modality == 'peptide':
        products = assemble_peptide_constructs([(full.antigen, selected)], PeptideConstructConfig(**options))
    else:
        products = assemble_mrna_constructs([(full.antigen, selected)],
            RNAConstructConfig(**options, optimize_linkers=False), host_taxon_id=9606)
    assert len(products) == 3
    assert [p.antigens[0]['aa'] for p in products] == [p.amino_acids for p in selected]
    assert all('window' in p.antigen_names[0] for p in products)


def test_native_full_antigen_enumeration_and_content_preserving_refinement():
    full, _ = windows()
    config = VaccineConfig(preferred_peptide_length=12, min_peptide_length=12, max_peptide_length=12)
    selected = optimize_antigen_windows([(full.antigen, [full])], config)[0][1]
    assert len(selected) == 3
    products = assemble_peptide_constructs([(full.antigen, selected)], PeptideConstructConfig(
        min_antigen_length_aa=8, max_antigen_length_aa=12, window_selection={'self_weight': .5}))
    assert len(products) == 3
    for product, original in zip(products, selected):
        assert all(epitope.sequence in product.sequence for epitope in original.target_epitopes)
        audit = product.components['window_selection']
        parent = audit['parent_window']
        assert full.amino_acids[parent['start']:parent['end']] == product.sequence
        assert audit['source_window']['index'] == original.window_selection_audit['source_window']['index']


def test_secondary_windows_compete_on_their_own_scores():
    full, choices = windows()
    selected = select_source_windows(choices)
    other = vaccine('AALLFGYPVYVAA', [(2, 11, 5., False)])
    other = replace(other, antigen=replace(other.antigen, source_identifier='another-source'))
    products = assemble_peptide_constructs([(full.antigen, selected), (other.antigen, [other])],
        PeptideConstructConfig(max_constructs=2, min_antigen_length_aa=8, max_antigen_length_aa=12))
    assert [p.sequence for p in products] == [selected[0].amino_acids, other.amino_acids]


def test_modality_trimming_preserves_selected_mutation_epitopes():
    from vaxrank.vaccine_antigen import TargetableMask, AminoAcidInterval
    antigen, peptides = candidate(0, 'SIINFEKLA', prefix='AAAA')
    antigen = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(6, 7),)))
    peptide = replace(peptides[0], antigen=antigen,
        window_selection_audit={'parent_window': dict(sequence=antigen.amino_acids, start=0,
                                                      end=len(antigen.amino_acids))})
    selected = select_source_windows([peptide])
    products = assemble_peptide_constructs([(antigen, selected)],
        PeptideConstructConfig(min_antigen_length_aa=9, max_antigen_length_aa=9))
    assert products[0].sequence == 'SIINFEKLA'


def test_native_rna_scoring_selects_measured_windows_before_ranking_parent():
    from types import SimpleNamespace
    from isovar import RNAFeature, measure_rna_feature
    from isovar.rna_path_linkage import OriginalPathMatch
    from vaxrank.epitope_dataset import EpitopeDataset, dataset_ranking_result
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import ExternalConstructOptions
    from vaxrank.rna_support import RNAFeatureEvidence
    full, _ = windows()
    records = tuple(RNAFeatureEvidence(measure_rna_feature(
        RNAFeature('window-%d' % start, 'synthetic-path', ((start * 3, (start + 12) * 3),), 'coding_interval'),
        [OriginalPathMatch(('RG', str(i), 0), 'placement', 0, len(full.amino_acids) * 3) for i in range(4)],
        evidence_scope=('synthetic-sample', 'synthetic-library')), 'window', (start, start + 12))
        for start in range(len(full.amino_acids) - 12 + 1))
    antigen = replace(full.antigen, rna_support=records)
    epitopes = [replace(e, prediction_id='synthetic-%d' % i) for i, e in enumerate(full.epitopes)]
    dataset = EpitopeDataset.from_predictions(epitopes, config=EpitopeConfig(),
        antigens={e.prediction_group_source: antigen for e in epitopes})
    options = ExternalConstructOptions.from_configs(VaccineConfig(
        preferred_peptide_length=12, min_peptide_length=12, max_peptide_length=12,
        combined_score_expr='sqrt(rna_fragments) * target_epitope_score'))
    result = dataset_ranking_result(SimpleNamespace(dataset=dataset, source_format='epitopes'), epitopes,
                                    options=options)
    selected = result.ranked[0][1]
    assert len(selected) == 3
    assert all(p.combined_score == 2 * p.target_epitope_score for p in selected)


def typed_candidate(index, kind, sequence):
    antigen, peptides = candidate(index, sequence)
    antigen = replace(antigen, kind=kind)
    return antigen, [replace(peptides[0], antigen=antigen)]


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_balanced_after_dedup_caps_goals_and_ranked_alternative(modality):
    ranked = [typed_candidate(0, 'CTA', 'SIINFEKLA'), typed_candidate(1, 'mutation', 'SIINFEKLA'),
              typed_candidate(2, 'CTA', 'GILGFVFTL'), typed_candidate(3, 'CTA', 'NLVPMVATV'),
              typed_candidate(4, 'mutation', 'GLLGATCMF'), typed_candidate(5, 'viral', 'LLFGYPVYV')]
    def build(policy):
        audit = {}
        kwargs = dict(max_constructs=3, antigens_per_construct=1, min_antigen_length_aa=9,
                      max_antigen_length_aa=13, antigen_allocation=policy)
        if modality == 'peptide':
            products = assemble_peptide_constructs(ranked, PeptideConstructConfig(**kwargs), selection_audit=audit)
        else:
            products = assemble_mrna_constructs(ranked, RNAConstructConfig(**kwargs, optimize_linkers=False),
                                               selection_audit=audit, host_taxon_id=9606)
        return products, audit
    balanced, audit = build(None)
    assert audit['allocation']['emitted_counts'] == {'CTA': 1, 'mutation': 1, 'viral': 1}
    assert len(balanced[0].antigens[0]['sources']) == 2
    assert [s['kind'] for s in balanced[0].antigens[0]['sources']] == ['CTA', 'mutation']
    _, plain = build({'mode': 'ranked'})
    assert plain['allocation']['emitted_counts'] == {'CTA': 3}
    _, limited = build({'mode': 'ranked', 'max_per_kind': {'CTA': 1},
                        'min_per_kind': {'viral': 2, 'ERV': 1}})
    assert limited['allocation']['emitted_counts'] == {'CTA': 1, 'mutation': 1, 'viral': 1}
    assert {d['kind'] for d in limited['allocation']['unmet_goals']} == {'viral', 'ERV'}
    assert sum(a['outcome'] == 'antigen_type_cap' for a in limited['antigens']) == 2


def test_cli_yaml_api_equivalence_and_frozen_default(tmp_path):
    parser = external_input_arg_parser()
    args = parser.parse_args(['--max-windows-per-antigen-source', '4', '--antigen-allocation', 'ranked',
                              '--min-antigens-per-type', 'mutation=2', '--max-antigens-per-type', 'CTA=3'])
    cli = vaccine_config_from_args(args)
    path = tmp_path / 'allocation.yaml'
    path.write_text('''vaccine_peptides:
  max_windows_per_antigen_source: 4
  antigen_allocation:
    mode: ranked
    min_per_kind: {mutation: 2}
    max_per_kind: {CTA: 3}
''')
    yaml_args = parser.parse_args(['--config', str(path)])
    assert cli == vaccine_config_from_args(yaml_args)
    assert _peptide_construct_options(args).antigen_allocation == _peptide_construct_options(yaml_args).antigen_allocation
    assert VaccineConfig().max_windows_per_antigen_source == 3
    legacy = VaccineConfig(**extract_vaccine_config_kwargs(load_vaxrank_config(config_path='builtin:openvax-v1')))
    assert legacy.max_windows_per_antigen_source is None
    assert legacy.antigen_allocation == {'mode': 'ranked'}
    with pytest.raises(ValueError, match='positive integer'):
        VaccineConfig(max_windows_per_antigen_source=0)
    with pytest.raises(ValueError, match='minimum exceeds maximum'):
        AntigenAllocation(min_per_kind={'CTA': 4}, max_per_kind={'CTA': 3})
    with pytest.raises(ValueError, match='Unknown antigen type'):
        AntigenAllocation(max_per_kind={'not_a_type': 1})


def test_oversubscribed_goals_report_capacity_and_zero_cap_excludes_kind():
    ranked = [typed_candidate(0, 'CTA', 'SIINFEKLA'), typed_candidate(1, 'CTA', 'GILGFVFTL'),
              typed_candidate(2, 'mutation', 'NLVPMVATV'), typed_candidate(3, 'mutation', 'GLLGATCMF')]
    audit = {}
    assemble_peptide_constructs(ranked, PeptideConstructConfig(max_constructs=1,
        antigen_allocation={'min_per_kind': {'CTA': 2, 'mutation': 2}}), selection_audit=audit)
    goals = audit['allocation']['unmet_goals']
    assert len(goals) == 2
    assert all(goal['reason'] == 'construct_capacity' and goal['available'] == 2 for goal in goals)
    audit = {}
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig(max_constructs=4,
        antigen_allocation={'max_per_kind': {'CTA': 0}}), selection_audit=audit)
    assert len(products) == 2
    assert audit['allocation']['emitted_counts'] == {'mutation': 2}


def mixed_source_table(path):
    """Synthetic, public test evidence; no patient or external model data."""
    import pandas as pd
    from topiary import TopiaryResult
    from vaxrank.epitope_dataset import ANTIGEN_COLUMN
    from vaxrank.native_serialization import to_native_json
    full, _ = windows()
    rows = []
    for kind, scale in [('CTA', 1.), ('viral', 2.), ('ERV', 3.)]:
        # Distinct flanks make different full products, even for shared ligands.
        sequence = (full.amino_acids if kind == 'CTA' else full.amino_acids.replace('L', 'Y')
                    if kind == 'viral' else full.amino_acids.replace('L', 'W'))
        antigen = replace(full.antigen, kind=kind, amino_acids=sequence,
                          gene_name='synthetic_' + kind, gene_id='',
                          source_identifier='synthetic-' + kind)
        for start, end, value in [(2, 10, 20.), (14, 23, 50.), (26, 35, 100.)]:
            rows.append(dict(peptide=sequence[start:end], allele='HLA-A*02:01',
                kind='pMHC_affinity', value=scale * value, prediction_method_name='synthetic',
                predictor_version='1', source_sequence=sequence, peptide_offset=start,
                prediction_id='%s-%d' % (kind, start), sample_name='synthetic-patient', host_species='human',
                **{ANTIGEN_COLUMN: to_native_json(antigen)}))
    TopiaryResult(pd.DataFrame(rows)).to_tsv(path)


def test_mixed_cli_report_and_native_replay(tmp_path, monkeypatch):
    from vaxrank.cli.entry_point import run_cli
    from vaxrank.epitope_dataset import EpitopeDataset
    from vaxrank.native_serialization import from_native_json
    from .test_epitope_dataset import forbidden
    source = tmp_path / 'mixed.tsv'
    mixed_source_table(source)
    out = tmp_path / 'design'
    run_cli(['--input-topiary', str(source), '--output-dir', str(out),
        '--output-html-report', str(out / 'review.html'), '--output-ascii-report', str(out / 'review.txt'),
        '--no-processing-aware-annotation', '--vaccine-peptide-length', '12',
        '--peptide-max-constructs', '5', '--mrna-max-constructs', '1',
        '--mrna-antigens-per-construct', '5', '--mrna-no-optimize-linkers',
        '--config-text', 'epitopes.score_expr=1 / affinity.value',
        '--min-antigens-per-type', 'mutation=1'])
    saved = EpitopeDataset.load(out / 'candidate_predictions.tsv')
    for modality, record in saved.selection['construct_assembly'].items():
        graph = from_native_json(record['source_graph'], dict)
        assert graph['assembly_decisions']['allocation']['emitted_counts'] == {'CTA': 2, 'viral': 2, 'ERV': 1}
        assert graph['assembly_decisions']['allocation']['unmet_goals'][0]['kind'] == 'mutation'
        assert sum(len(p.antigen_names) for p in from_native_json(record['products'], list)) == 5
    assert 'Emitted antigen sequences by type:' in (out / 'review.html').read_text()
    assert 'Unmet mutation goal' in (out / 'review.txt').read_text()
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    run_cli(['--input-epitopes', str(out / 'candidate_predictions.tsv'), '--output-dir', str(tmp_path / 'replay'),
             '--no-processing-aware-annotation'])
    for first in out.rglob('*.fasta'):
        assert first.read_bytes() == (tmp_path / 'replay' / first.relative_to(out)).read_bytes()
