from dataclasses import replace

import pytest
from mhctools.pred import Prediction

from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.core_logic import SOURCE_AGNOSTIC_RANKING_RULES
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.vaccine_peptide import VaccinePeptide
from vaxrank.window_selection import WindowSelection, window_metrics, select_windows
from .test_vaccine_antigen import cta_antigen


def vaccine(sequence, intervals):
    antigen = cta_antigen(sequence)
    epitopes = []
    for start, end, score, self_match in intervals:
        peptide = sequence[start:end]
        epitopes.append(CandidateEpitope(
            sequence=peptide, source_sequence=sequence, offset=start,
            overlaps_targetable=not self_match, occurs_in_reference=self_match,
            occurs_in_non_CTA_reference=self_match,
            per_allele_scores={'HLA-A*02:01': score},
            predictions=(Prediction(kind='pMHC_affinity', peptide=peptide,
                                    allele='HLA-A*02:01', value=100., score=.8),)))
    return VaccinePeptide(antigen=antigen, epitopes=epitopes,
                          combined_score_expr='window_epitope_score',
                          ranking_rules=SOURCE_AGNOSTIC_RANKING_RULES)


def test_serum_penalty_only_internal_target_cuts():
    vp = vaccine('AAGPNQSIINFEKL', [(2, 11, 1., False), (4, 13, .8, False)])
    policy = WindowSelection(serum_models=('fap-endo-gp',), serum_weight=.25)
    metric = window_metrics(vp, policy)
    assert metric['cuts'] == [4]
    assert metric['target_score'] == 1.8
    assert metric['at_risk_target_score'] == 1.
    assert metric['utility'] == 1.55
    assert metric['cleavage_profiles'][0]['model']['evidence'] == 'motif_rule'


def test_intact_duplicate_copy_preserves_epitope_content():
    vp = vaccine('ASIINFEKLASIINFEKL', [(0, 9, 1., False), (9, 18, 1., False)])
    policy = WindowSelection(serum_models=('dpp4-qpisa',), serum_weight=.5,
                             dpp4_depletion_threshold=0.)
    metric = window_metrics(vp, policy)
    assert metric['target_score'] == 1.
    assert metric['at_risk_target_score'] == 0.


def test_non_cta_self_burden_changes_choice_without_losing_target():
    contaminated = vaccine('AAAAASIINFEKLAAAAA', [(5, 13, 1., False), (0, 8, 2., True)])
    clean = vaccine('SIINFEKLAAAAA', [(0, 8, 1., False)])
    weak = vaccine('SIINFEKLAAA', [(0, 8, .5, False)])
    chosen = select_windows([contaminated, weak, clean], WindowSelection(self_weight=.5),
                             preferred_length=17, limit=1)
    assert chosen == [clean]
    assert chosen[0].window_selection_audit['alternatives'][1]['eligible'] is False
    assert chosen[0].combined_score == 1.


def test_peptide_assembly_trims_self_ligand_and_retains_audit():
    vp = vaccine('AAAAASIINFEKLAAAAA', [(5, 13, 1., False), (0, 8, 2., True)])
    options = PeptideConstructConfig(min_antigen_length_aa=10, max_antigen_length_aa=17,
                                     window_selection={'self_weight': .5})
    constructs = assemble_peptide_constructs([(vp.antigen, [vp])], options)
    assert len(constructs) == 1
    assert len(constructs[0].sequence) < len(vp.amino_acids)
    assert 'SIINFEKL' in constructs[0].sequence
    assert constructs[0].components['window_selection']['non_cta_self_score'] == 0.
    assert vp.window_selection_audit == {}  # shared/mRNA input stays unchanged


def test_terminal_chemistry_and_unknown_coverage_remain_visible():
    vp = vaccine('ASIINFEKL', [(0, 9, 1., False)])
    policy = WindowSelection(serum_models=('dpp4-qpisa',), serum_weight=.5)
    metric = window_metrics(vp, policy, n_term='acetylated')
    assert metric['unassessed_models'] == ['dpp4-qpisa']
    assert not metric['cuts']
    assert not metric['self_provenance_complete']
    with pytest.raises(ValueError, match='one SLP'):
        replace(PeptideConstructConfig(), antigens_per_construct=2,
                window_selection={'serum_weight': .25})


def test_cli_reads_peptide_policy_and_writes_rejected_window_decisions(tmp_path):
    import json
    from vaxrank.cli import make_vaxrank_arg_parser
    from vaxrank.cli.entry_point import _emit_peptide_constructs
    config = tmp_path / 'settings.yaml'
    config.write_text('''peptide:
  window_selection:
    self_weight: 0.5
vaccine_peptides:
  min_antigen_length_aa: 10
  max_antigen_length_aa: 17
''')
    args = make_vaxrank_arg_parser().parse_args([
        '--vcf', 'unused', '--bam', 'unused', '--mhc-predictor', 'random',
        '--mhc-alleles', 'HLA-A*02:01', '--config', str(config)])
    vp = vaccine('AAAAASIINFEKLAAAAA', [(5, 13, 1., False), (0, 8, 2., True)])
    output = tmp_path / 'constructs'
    _emit_peptide_constructs(args, [(vp.antigen, [vp])], str(output))
    manifest = json.loads((output / 'manifest.json').read_text())
    assert manifest[0]['components']['window_selection']['non_cta_self_score'] == 0.
    audit = json.loads((output / 'window_selection.json').read_text())
    assert audit[0]['selected'] == [manifest[0]['sequence']]
    assert any(not row['eligible'] for row in audit[0]['alternatives'])
