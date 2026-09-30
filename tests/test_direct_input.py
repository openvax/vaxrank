"""Direct prediction is one scoped producer in the shared evidence workflow."""

import pytest
from types import SimpleNamespace
from pyensembl import EnsemblRelease
from varcode import Variant, VariantCollection

from vaxrank.cli.arg_parser import parse_vaxrank_args
from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.external_input import load_external_ranked
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_config import VaccineConfig
from .input_scope_helpers import write_input_manifest
from .test_core_logic_config import _make_epitope
from .test_epitope_dataset import evidence, forbidden


@pytest.fixture
def direct(tmp_path, monkeypatch):
    import varcode.cli
    import mhctools.cli
    import vaxrank.cli.entry_point as entry
    variant = Variant('1', 100, 'A', 'T', ensembl=EnsemblRelease(110))
    fragment = MutantProteinFragment(
        variant=variant, gene_name='DIRECT', amino_acids='AASIINFEKLKK',
        mutant_amino_acid_start_offset=5, mutant_amino_acid_end_offset=6,
        supporting_reference_transcripts=[], n_overlapping_reads=30,
        n_alt_reads=12, n_ref_reads=18, n_alt_reads_supporting_protein_sequence=12,
        n_overlapping_fragments=20, n_alt_fragments=8, n_ref_fragments=12,
        n_alt_fragments_supporting_protein_sequence=8,
        rna_evidence_method='isovar')
    candidates = [_make_epitope('SIINFEKL', 25., source_sequence=fragment.amino_acids, offset=2)]
    vcf, bam = tmp_path / 'direct.vcf', tmp_path / 'direct.bam'
    vcf.write_text('synthetic VCF')
    bam.write_bytes(b'synthetic BAM')
    calls = []
    def predict(args, *, epitope_dataset, epitope_config_override):
        calls.append(epitope_config_override)
        epitope_dataset.add_mutation(fragment, candidates)
        return SimpleNamespace(
            isovar_results=[SimpleNamespace(variant=variant)],
            variant_properties=lambda **kwargs: [{'is_coding_nonsynonymous': True,
                                                  'rna_support': True, 'dna_vaf': None}])
    monkeypatch.setattr(entry, 'run_vaxrank_from_parsed_args', predict)
    monkeypatch.setattr(varcode.cli, 'variant_collection_from_args',
                        lambda args: VariantCollection([variant]))
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    monkeypatch.setattr(entry, 'annotate_predictions_with_processing', forbidden)
    return vcf, bam, fragment, calls


def inputs(tmp_path, direct, **scope):
    vcf, bam, _, _ = direct
    table = tmp_path / 'table.tsv'
    evidence().to_tsv(table)
    manifest = write_input_manifest(tmp_path / 'inputs.json', [('topiary', table)], **scope)
    return ['--vcf', str(vcf), '--bam', str(bam), '--mhc-predictor', 'random',
            '--mhc-alleles', 'HLA-A*02:01', '--ensembl-release', '110',
            '--input-manifest', manifest, '--duplicate-candidates', 'best']


def test_mixed_cli_and_native_reload_keep_original_values_and_direct_context(direct, tmp_path):
    native, reloaded = tmp_path / 'saved.tsv', tmp_path / 'again.tsv'
    config = ['--config-value', 'vaccine_peptides.min_length=8',
              '--config-value', 'vaccine_peptides.preferred_length=12']
    run_cli(inputs(tmp_path, direct) + config + [
        '--output-epitopes', str(native), '--output-json-file', str(tmp_path / 'run.json'),
        '--output-passing-variants-csv', str(tmp_path / 'variants.csv')])
    import json
    assert json.loads((tmp_path / 'run.json').read_text())
    assert 'has_vaccine_peptide' in (tmp_path / 'variants.csv').read_text()
    dataset = EpitopeDataset.load(native)
    assert len(direct[3]) == 1
    assert len(dataset.epitopes) == 3
    assert sorted(dataset.result.df.value) == [25., 50., 500.]
    assert len(dataset.mutation_fragments) == 1
    fragment = next(iter(dataset.mutation_fragments.values()))
    assert fragment.n_rna_alt == 8
    assert fragment.rna_evidence_method == 'isovar'
    assert fragment.variant == direct[2].variant
    assert next(p for p in dataset.provenance if p.source_format == 'vcf_bam').scope.annotation == 'ensembl:110'
    run_cli(['--input-epitopes', str(native), '--output-epitopes', str(reloaded)] + config)
    again = EpitopeDataset.load(reloaded)
    assert [e.per_allele_scores for e in again.epitopes] == [e.per_allele_scores for e in dataset.epitopes]
    assert sorted(again.result.df.value) == [25., 50., 500.]
    assert len(direct[3]) == 1
    args = parse_vaxrank_args(['--input-epitopes', str(reloaded)])
    vc = VaccineConfig(preferred_peptide_length=12, min_peptide_length=8)
    ranked, _, _, _, _ = load_external_ranked(args, vaccine_config=vc)
    assert len(ranked) == 1
    assert ranked[0][1][0].mutant_protein_fragment.n_rna_alt == 8
    assert ranked[0][1][0].amino_acids == fragment.amino_acids


@pytest.mark.parametrize('scope, message', [
    ({'patient_id': None}, 'patient_id'),
    ({'reference_assembly': 'GRCh37'}, 'reference_assembly'),
    ({'annotation': 'ensembl:109'}, 'annotation'),
    ({'mhc_alleles': ['HLA-B*07:02']}, 'alleles outside'),
])
def test_scope_conflicts_precede_direct_prediction(direct, tmp_path, scope, message):
    with pytest.raises(ValueError, match=message):
        run_cli(inputs(tmp_path, direct, **scope) + ['--output-epitopes', str(tmp_path / 'no.tsv')])
    assert direct[3] == []
    assert not (tmp_path / 'no.tsv').exists()


def test_direct_prediction_flags_do_not_rescore_tables(direct, tmp_path):
    with pytest.raises(ValueError, match='external-predictions input'):
        run_cli(inputs(tmp_path, direct) + ['--external-predictions', 'fresh',
                                          '--output-epitopes', str(tmp_path / 'no.tsv')])
    assert direct[3] == []


@pytest.mark.parametrize('fmt', ['lens', 'pvacseq'])
def test_report_constructs_survive_common_native_export(tmp_path, fmt):
    from vaxrank.external_rescoring import prepare_reports
    from .test_external_rescoring import DATA
    source = DATA / (fmt + '_example.tsv')
    _, frame, _ = prepare_reports([(fmt, source)])
    native = tmp_path / 'saved.tsv'
    frame.attrs['epitope_dataset'].save(native)
    options = VaccineConfig(combined_score_expr='target_epitope_score')
    before = load_external_ranked(parse_vaxrank_args(['--external-input', f'{fmt}={source}']),
                                  vaccine_config=options)
    after = load_external_ranked(parse_vaxrank_args(['--input-epitopes', str(native)]),
                                 vaccine_config=options)
    assert [(str(v), p[0].amino_acids, p[0].combined_score) for v, p in before[0]] == [
        (str(v), p[0].amino_acids, p[0].combined_score) for v, p in after[0]]
    assert before[0]
    assert before[3].num_somatic_variants == after[3].num_somatic_variants


def test_shared_dsl_changes_duplicate_selection_and_constructs(direct, tmp_path):
    from vaxrank.epitope_dataset import ANTIGEN_COLUMN
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.native_serialization import to_native_json
    from vaxrank.vaccine_antigen import (
        AminoAcidInterval, TargetableMask, TumorSpecificityAttestation, VaccineAntigen,
    )
    from vaxrank.peptide import assemble_peptide_constructs, PeptideConstructConfig
    from vaxrank.mrna import assemble_mrna_constructs, RNAConstructConfig
    argv = inputs(tmp_path, direct)
    antigens = [VaccineAntigen(
        kind='viral', amino_acids=sequence,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture', evidence_source=__name__),
        source_identifier='table-' + sequence) for sequence in ('SIINFEKL', 'GILGFVFTL')]
    evidence(**{ANTIGEN_COLUMN: [to_native_json(a) for a in antigens]}).to_tsv(tmp_path / 'table.tsv')
    options = VaccineConfig(preferred_peptide_length=12, min_peptide_length=8,
                            combined_score_expr='target_epitope_score')
    selections = []
    for index, expression in enumerate(('1 / affinity.value', 'affinity.value')):
        args = parse_vaxrank_args(argv)
        # Real run_cli resolves the release before loading inputs.
        args.genome = EnsemblRelease(110)
        ranked, frame, _, patient, _ = load_external_ranked(
            args, epitope_config=EpitopeConfig(score_expr=expression, min_epitope_score=0),
            vaccine_config=options)
        native = tmp_path / f'policy-{index}.tsv'
        frame.attrs['epitope_dataset'].save(native)
        reloaded = load_external_ranked(parse_vaxrank_args(['--input-epitopes', str(native)]),
                                       vaccine_config=options)
        assert len(ranked) == len(reloaded[0]) == 2
        assert patient.num_somatic_variants == reloaded[3].num_somatic_variants
        # Direct + normalized discovery of SIINFEKL yields one credited
        # occurrence, never two doses or summed RNA support.
        selected = frame.loc[frame['Selected observation']]
        assert selected['Mutant peptide sequence'].tolist().count('SIINFEKL') == 1
        for assembly, config in (
            (assemble_peptide_constructs, PeptideConstructConfig(
                mode='minimal_epitope', min_antigen_length_aa=8)),
            (assemble_mrna_constructs, RNAConstructConfig(
                antigen_content='minimal_epitope', optimize_linkers=False)),
        ):
            original = assembly(ranked, options=config)
            again = assembly(reloaded[0], options=config)
            assert original and [c.sequence for c in original] == [c.sequence for c in again]
        selections.append([str(source) for source, _ in ranked])
    assert selections[0] != selections[1]
    assert len(direct[3]) == 2


def test_direct_and_both_report_adapters_keep_scope_and_counts(direct, tmp_path):
    from .test_external_rescoring import INPUTS
    argv = inputs(tmp_path, direct)
    write_input_manifest(tmp_path / 'inputs.json', INPUTS,
                         direct={'sample_id': 'direct-tumor', 'library_id': 'rna-1'})
    native = tmp_path / 'all.tsv'
    run_cli(argv + ['--output-epitopes', str(native)])
    dataset = EpitopeDataset.load(native)
    assert {p.source_format for p in dataset.provenance} == {'vcf_bam', 'lens', 'pvacseq'}
    direct_scope = next(p.scope for p in dataset.provenance if p.source_format == 'vcf_bam')
    assert direct_scope.sample_id == 'direct-tumor'
    assert direct_scope.library_id == 'rna-1'
    assert len(dataset.construction_reports) == 2
    assert len(dataset.direct_sources) == 1
    assert next(iter(dataset.mutation_fragments.values())).n_rna_alt == 8
    _, audit, _, _, _ = load_external_ranked(
        parse_vaxrank_args(['--input-epitopes', str(native)]))
    legacy = audit.loc[audit.input_format.isin(['lens', 'pvacseq'])]
    assert not legacy.empty
    assert legacy['Selected observation'].isna().all()


def test_mixed_coverage_uses_declared_genotype():
    from vaxrank.cli.entry_point import resolve_target_alleles
    args = parse_vaxrank_args(['--vcf', 'test.vcf', '--bam', 'test.bam',
                              '--mhc-predictor', 'random', '--mhc-alleles', 'HLA-A*02:01'])
    args._declared_mhc_alleles = ['HLA-A*02:01', 'HLA-B*07:02']
    assert resolve_target_alleles(args) == args._declared_mhc_alleles


def test_native_direct_fragment_preserves_structural_variant(direct):
    from dataclasses import replace
    from varcode import StructuralVariant
    from vaxrank.native_serialization import from_native_json, to_native_json
    variant = StructuralVariant('1', 100, 'DEL', end=200, genome=EnsemblRelease(110))
    fragment = replace(direct[2], variant=variant)
    loaded = from_native_json(to_native_json(fragment), MutantProteinFragment)
    assert isinstance(loaded.variant, StructuralVariant)
    assert loaded.variant == variant
    assert loaded.amino_acids == fragment.amino_acids
    assert loaded.n_rna_alt == 8


def test_b16_mixed_cli_assembly_and_reload(tmp_path, mouse_genome, monkeypatch):
    """Exercise real VCF/BAM reconstruction and random-model prediction."""
    import numpy as np
    from .testing_helpers import data_path
    from vaxrank.peptide import assemble_peptide_constructs, PeptideConstructConfig
    from vaxrank.mrna import assemble_mrna_constructs, RNAConstructConfig
    table = tmp_path / 'mouse-table.tsv'
    evidence(allele=['H2-Kb', 'H2-Kb']).to_tsv(table)
    manifest = write_input_manifest(
        tmp_path / 'mouse-inputs.json', [('topiary', table)],
        reference_assembly=mouse_genome.reference_name, mhc_alleles=['H2-Kb', 'H2-Db'])
    native = tmp_path / 'mixed.tsv'
    common = ['--vaccine-peptide-length', '15', '--no-processing-aware-annotation',
              '--config-text', 'epitopes.score_expr=1 / affinity.value',
              '--min-epitope-score', '0']
    run_cli(['--vcf', data_path('b16.f10/b16.vcf'),
             '--bam', data_path('b16.f10/b16.combined.bam'),
             '--mhc-predictor', 'random', '--mhc-alleles', 'H2-Kb,H2-Db',
             '--mhc-epitope-lengths', '8', '--padding-around-mutation', '0',
             '--input-manifest', manifest, '--output-epitopes', str(native)] + common)
    import vaxrank.cli.entry_point as entry
    import mhctools.cli
    monkeypatch.setattr(entry, 'predictors_from_args', forbidden)
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    reloaded = tmp_path / 'again.tsv'
    run_cli(['--input-epitopes', str(native), '--output-epitopes', str(reloaded)] + common)
    before, after = EpitopeDataset.load(native), EpitopeDataset.load(reloaded)
    assert len(before.epitopes) > 2
    assert len(before.direct_sources) == 5
    assert before.mutation_fragments
    for first, second in zip(before.epitopes, after.epitopes):
        assert first.prediction_id == second.prediction_id
        assert first.per_allele_scores.keys() == second.per_allele_scores.keys()
        # Topiary #441 tracks decimal CSV/TSV parser rounding (native object
        # scores themselves are exact). Check the recomputed numeric result.
        np.testing.assert_allclose(list(first.per_allele_scores.values()),
                                   list(second.per_allele_scores.values()), rtol=1e-14, atol=0)
    options = VaccineConfig(preferred_peptide_length=15, min_peptide_length=15,
                            max_peptide_length=15)
    sequences = []
    for path in (native, reloaded):
        ranked, _, _, patient, _ = load_external_ranked(
            parse_vaxrank_args(['--input-epitopes', str(path)]), vaccine_config=options)
        assert patient.num_somatic_variants == 5
        peptide = assemble_peptide_constructs(ranked, options=PeptideConstructConfig())
        mrna = assemble_mrna_constructs(ranked, options=RNAConstructConfig(optimize_linkers=False))
        assert peptide and mrna
        sequences.append(([c.sequence for c in peptide], [c.sequence for c in mrna]))
    assert sequences[0] == sequences[1]


def test_direct_configured_alternative_windows_survive_reload(direct, tmp_path, monkeypatch):
    from dataclasses import replace
    import vaxrank.cli.entry_point as entry
    fragment = replace(direct[2], amino_acids='AAAAASIINFEKLAAAAAAA',
                       mutant_amino_acid_start_offset=8, mutant_amino_acid_end_offset=9)
    candidate = _make_epitope('SIINFEKL', 25., source_sequence=fragment.amino_acids, offset=5)
    def predict(args, *, epitope_dataset, epitope_config_override):
        epitope_dataset.add_mutation(fragment, [candidate])
        return SimpleNamespace(isovar_results=[], variant_properties=lambda **kwargs: [])
    monkeypatch.setattr(entry, 'run_vaxrank_from_parsed_args', predict)
    args = parse_vaxrank_args(inputs(tmp_path, direct))
    options = VaccineConfig(preferred_peptide_length=12, min_peptide_length=12,
                            max_peptide_length=12, max_vaccine_peptides_per_variant=3)
    ranked, frame, _, _, _ = load_external_ranked(args, vaccine_config=options)
    assert len(ranked) == 1
    assert len(ranked[0][1]) == 3
    native = tmp_path / 'alternatives.tsv'
    frame.attrs['epitope_dataset'].save(native)
    again, _, _, _, _ = load_external_ranked(
        parse_vaxrank_args(['--input-epitopes', str(native)]), vaccine_config=options)
    assert [p.amino_acids for p in ranked[0][1]] == [p.amino_acids for p in again[0][1]]
