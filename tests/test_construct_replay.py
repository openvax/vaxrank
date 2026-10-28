"""General native designs retain actual products, settings and source graphs."""

from dataclasses import replace
import gzip
from pathlib import Path

import pytest

from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import ANTIGEN_COLUMN, EpitopeDataset
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation, VaccineAntigen,
)

from .test_epitope_dataset import evidence, forbidden
from .test_external_rescoring import ContextPredictor


@pytest.fixture(autouse=True)
def assembly_reports(monkeypatch):
    # These tests inspect the real modality writers and JSON/native contracts.
    # Report rendering has separate coverage and an end-to-end CLI smoke.
    for name in ('make_ascii_report', 'make_html_report', 'make_pdf_report'):
        monkeypatch.setattr('vaxrank.cli.entry_point.' + name, lambda *args, **kwargs: None)


@pytest.fixture
def source(tmp_path):
    sequences = ['AASIINFEKLKKKKKKK', 'AAGILGFVFTLKKKKKK']
    antigens = [VaccineAntigen(
        kind='viral', amino_acids=sequence, gene_name='SOURCE_%d' % i,
        source_identifier='source-%d' % i,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture',
            evidence_source='test_construct_replay', patient_specific=True,
            rationale_code='test_only')) for i, sequence in enumerate(sequences)]
    path = tmp_path / 'source.tsv'
    evidence(**{ANTIGEN_COLUMN: [to_native_json(a) for a in antigens],
                'source_sequence': sequences, 'peptide_offset': [2, 2],
                'sample_name': ['patient', 'patient']}).to_tsv(path)
    return path


def design(source, directory, native, *options):
    run_cli(['--input-topiary', str(source), '--output-dir', str(directory),
             '--output-epitopes', str(native), '--no-processing-aware-annotation',
             '--config-text', 'epitopes.score_expr=1.0', *options])


def replay(native, directory, *options):
    run_cli(['--input-epitopes', str(native), '--output-dir', str(directory),
             '--no-processing-aware-annotation', *options])


def product_files(directory):
    return {str(path.relative_to(directory)): path.read_bytes()
            for path in directory.rglob('*') if path.is_file() and path.suffix in ('.fasta', '.csv', '.json')
            and path.parent.name in ('peptide', 'mrna')}


@pytest.mark.parametrize('suffix', ['csv', 'tsv'])
def test_empty_native_result_replays_scope_settings_and_products(tmp_path, monkeypatch, suffix):
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_rescoring import prepare_reports
    from vaxrank.input_scope import InputProvenance, InputScope

    scope = InputScope(patient_id='empty-patient', reference_assembly='GRCh38',
                       mhc_alleles=('HLA-A*02:01',))
    original = EpitopeDataset.from_predictions((), config=EpitopeConfig(score_expr='1.0'))
    original.provenance = (InputProvenance('empty-input', 'vcf_bam', 'input.vcf', 'a' * 64, scope),)
    native = tmp_path / ('empty.' + suffix)
    original.save(native)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.mhc_binding_predictor_from_args', forbidden)
    _, frame, candidates = prepare_reports([('epitopes', native)])
    loaded = frame.attrs['epitope_dataset']
    assert not candidates and loaded.result.df.empty
    assert loaded.provenance == original.provenance
    assert loaded.config == original.config
    assert not loaded.antigens and not loaded.selection.get('source_evidence')

    first, second = tmp_path / 'first', tmp_path / 'second'
    replay(native, first, '--peptide-max-constructs', '3', '--mrna-max-constructs', '2')
    saved = EpitopeDataset.load(first / 'candidate_predictions.tsv')
    assembly = saved.selection['construct_assembly']
    assert set(assembly) == {'peptide', 'mrna'}
    assert assembly['peptide']['configuration']['max_constructs'] == 3
    assert assembly['mrna']['configuration']['max_constructs'] == 2
    for record in assembly.values():
        assert from_native_json(record['products'], list) == []
        assert from_native_json(record['source_graph'], dict) == {'inputs': [], 'products': {}}
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    replay(first / 'candidate_predictions.tsv', second)
    restored = EpitopeDataset.load(second / 'candidate_predictions.tsv')
    assert restored.provenance == original.provenance
    assert restored.config == original.config
    assert restored.selection['construct_assembly'] == assembly
    assert restored.result.df.empty and not restored.epitopes and not restored.antigens
    files = product_files(first)
    assert files and files == product_files(second)
    assert all(not contents for path, contents in files.items() if path.endswith('.fasta'))


def test_pre_species_design_restores_both_products_without_reassembly(tmp_path, monkeypatch):
    # Actual 3.44.3 export (3eaa3d0), generated with this module's source/design
    # helpers and --mrna-no-optimize-linkers --mrna-codon-method use_best_codon.
    # Relative input paths keep the fixture independent of the original machine.
    native = tmp_path / 'legacy.tsv'
    native.write_bytes(gzip.decompress(
        (Path(__file__).parent / 'data/epitope_fixtures/native_design_3_44_3.tsv.gz').read_bytes()))
    saved = EpitopeDataset.load(native).selection['construct_assembly']
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    output = tmp_path / 'replay'
    replay(native, output)
    restored = EpitopeDataset.load(output / 'candidate_predictions.tsv').selection['construct_assembly']
    for modality, original in saved.items():
        assert restored[modality]['identity'] == original['identity']
        assert restored[modality]['payload_sha256'] == original['payload_sha256']
        assert from_native_json(restored[modality]['products'], list) == from_native_json(original['products'], list)
    assert {path.parent.name for path in output.rglob('*.fasta')} == {'peptide', 'mrna'}


def test_legacy_duplicate_slots_rebuild_once_then_replay_offline(tmp_path, monkeypatch):
    # Actual bf1634b / 3.47.3 export: different source windows and predicted
    # ligands both emit SIINFEKLA. Previously they occupied two design slots.
    native = tmp_path / 'legacy.tsv'
    native.write_bytes(gzip.decompress((Path(__file__).parent /
        'data/epitope_fixtures/native_duplicate_design_3_47_3.tsv.gz').read_bytes()))
    before = EpitopeDataset.load(native)
    assert len(from_native_json(before.selection['construct_assembly']['peptide']['products'], list)) == 2
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    first, second = tmp_path / 'corrected', tmp_path / 'again'
    replay(native, first)
    corrected = EpitopeDataset.load(first / 'candidate_predictions.tsv')
    assert corrected.provenance == before.provenance
    assert corrected.antigens == before.antigens
    assert len(corrected.epitopes) == len(before.epitopes) == 2
    for modality in ('peptide', 'mrna'):
        assembly = corrected.selection['construct_assembly'][modality]
        product, = from_native_json(assembly['products'], list)
        assert len(product.antigen_names) == 1
        sources, = product.antigen_sources.values()
        assert {s['source_identifier'] for s in sources} == {'source-0', 'source-1'}
        assert {s['gene_ids'][0] for s in sources} == {'gene-0', 'gene-1'}
        graph = from_native_json(assembly['source_graph'], dict)
        assert len(graph['inputs']) == 2
        assert all(c['emitted_antigen_names'] == product.antigen_names
                   for row in graph['inputs'] for c in row['candidates'])
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    replay(first / 'candidate_predictions.tsv', second)
    assert product_files(first) == product_files(second)


@pytest.mark.parametrize('optimized', [False, True])
@pytest.mark.parametrize('content', ['mutation_spanning', 'minimal_epitope'])
def test_actual_general_products_replay_offline(source, tmp_path, monkeypatch, optimized, content):
    model = ContextPredictor()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', lambda args: [model])
    monkeypatch.setattr('vaxrank.cli.entry_point.mhc_binding_predictor_from_args', lambda args: model)
    native, first, second = tmp_path / 'native.tsv', tmp_path / 'first', tmp_path / 'replay'
    options = ['--peptide-max-constructs', '1', '--peptide-n-terminal-acetyl',
               '--peptide-c-terminal-amide', '--mrna-poly-a-length', '73',
               '--mrna-poly-a-segmented', '--mrna-antigens-per-construct', '2',
               '--mrna-max-constructs', '1', '--mrna-codon-method', 'use_best_codon',
               '--peptide-antigen-content', content, '--mrna-antigen-content', content]
    if content == 'mutation_spanning':
        options += ['--config-value', 'peptide.window_selection={"min_target_fraction":0.95}']
    if optimized:
        options += ['--mrna-junction-predictor', 'random', '--mrna-junction-alleles',
                    'HLA-A*02:01,HLA-B*07:02']
    else:
        options += ['--mrna-no-optimize-linkers']
    design(source, first, native, *options)
    saved = EpitopeDataset.load(native)
    assert 'cta_expression_admission' not in saved.selection
    assembly = saved.selection['construct_assembly']
    peptide, = from_native_json(assembly['peptide']['products'], list)
    mrna, = from_native_json(assembly['mrna']['products'], list)
    assert peptide.components['n_terminal_acetylation']
    assert peptide.components['c_terminal_amidation']
    assert len(mrna.poly_a_nt) == 83  # 73 A bases plus the 10-base segmented linker
    graph = from_native_json(assembly['mrna']['source_graph'], dict)
    assert len(graph['inputs']) == 2
    assert graph['products'][mrna.name] == mrna.antigen_names
    assert set(mrna.antigen_names) == {name for row in graph['inputs'] for candidate in row['candidates']
                                     for name in candidate['emitted_antigen_names']}
    assert {row['candidates'][0]['antigen'].source_identifier for row in graph['inputs']} == {'source-0', 'source-1'}
    peptide_graph = from_native_json(assembly['peptide']['source_graph'], dict)
    assert set(peptide.antigen_names) == {name for row in peptide_graph['inputs'] for c in row['candidates']
                                        for name in c['emitted_antigen_names']}
    assert assembly['mrna']['configuration']['optimize_linkers'] == optimized
    if optimized:
        assert mrna.elements['junction_swap']['enabled']
        assert mrna.elements['junction_swap']['prediction']['scored_prediction_count'] > 0
    assert saved.selection['run_configuration']['effective_constructs']['mrna']['poly_a_length'] == 73
    assert bool(model.calls) == optimized
    source.unlink()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.mhc_binding_predictor_from_args', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    monkeypatch.setattr('vaxrank.mrna.reverse_translate', forbidden)
    replay(native, second)
    assert product_files(first) == product_files(second)


def test_current_yaml_and_explicit_default_cli_options_override_saved_settings(source, tmp_path):
    native, first = tmp_path / 'native.tsv', tmp_path / 'first'
    design(source, first, native, '--peptide-max-constructs', '1',
           '--peptide-n-terminal-acetyl', '--mrna-poly-a-length', '73',
           '--mrna-no-optimize-linkers')
    changed, derived = tmp_path / 'changed', tmp_path / 'derived.tsv'
    config = tmp_path / 'override.yaml'
    config.write_text('peptide:\n  max_constructs: 20\n  n_terminal_acetyl: false\n'
                      'mrna:\n  poly_a_length: 99\n  utr_3p: HBB\n')
    replay(native, changed, '--config', str(config), '--mrna-poly-a-length=120',
           '--output-epitopes', str(derived))
    saved = EpitopeDataset.load(derived)
    products = saved.selection['construct_assembly']
    peptides = from_native_json(products['peptide']['products'], list)
    assert len(peptides) == 2
    assert not any(p.components['n_terminal_acetylation'] for p in peptides)
    assert all(p.poly_a_nt == 'A' * 120 for p in from_native_json(products['mrna']['products'], list))
    assert products['mrna']['configuration']['utr_3p'] == 'HBB'
    original = EpitopeDataset.load(native).selection['construct_assembly']
    assert all(products[m]['identity'] != original[m]['identity'] for m in products)


@pytest.mark.parametrize('change', ['policy', 'source', 'capacity'])
def test_changed_design_inputs_do_not_reuse_products_or_window_audits(source, tmp_path, monkeypatch, change):
    import vaxrank.cli.entry_point as entry
    native = tmp_path / 'native.tsv'
    design(source, tmp_path / 'first', native, '--mrna-no-optimize-linkers',
           '--config-value', 'peptide.window_selection={"min_target_fraction":0.95}')
    if change == 'source':
        dataset = EpitopeDataset.load(native)
        dataset.antigens = {key: replace(antigen, gene_name='CHANGED_SOURCE')
                            for key, antigen in dataset.antigens.items()}
        dataset.save(native)
    calls = []
    for name in ('assemble_peptide_constructs', 'assemble_mrna_constructs'):
        original = getattr(entry, name)
        def capture(*args, _name=name, _original=original, **kwargs):
            calls.append(_name)
            return _original(*args, **kwargs)
        monkeypatch.setattr(entry, name, capture)
    options = (['--config-text', 'epitopes.score_expr=2.0 - 1.0'] if change == 'policy'
               else ['--peptide-max-constructs', '1', '--mrna-max-length-nt', '100']
               if change == 'capacity' else [])
    replay(native, tmp_path / 'changed', *options)
    assert calls == ['assemble_peptide_constructs', 'assemble_mrna_constructs']


def test_corrupt_saved_products_are_rejected(source, tmp_path):
    native = tmp_path / 'native.tsv'
    design(source, tmp_path / 'first', native, '--mrna-no-optimize-linkers')
    saved = EpitopeDataset.load(native)
    payload = saved.selection['construct_assembly']['peptide']
    products = from_native_json(payload['products'], list)
    products[0].components['counterion'] = 'changed'
    payload['products'] = to_native_json(products)
    saved.save(native)
    with pytest.raises(ValueError, match='checksum'):
        replay(native, tmp_path / 'replay')


def test_legacy_native_evidence_does_not_invent_saved_assembly_settings(source, tmp_path):
    dataset = EpitopeDataset.from_topiary(__import__('topiary').read_tsv(source))
    native, directory = tmp_path / 'legacy.tsv', tmp_path / 'replay'
    dataset.save(native)
    replay(native, directory, '--output-epitopes', str(tmp_path / 'modern.tsv'),
           '--mrna-no-optimize-linkers')
    saved = EpitopeDataset.load(tmp_path / 'modern.tsv')
    assert saved.selection['construct_assembly']['mrna']['configuration']['poly_a_length'] == 120
    assert saved.selection['run_configuration']['effective_constructs']
    assert not (directory / 'selection_policy.json').exists()


def test_saved_source_limit_replays_without_losing_unselected_evidence(source, tmp_path, monkeypatch):
    native, first = tmp_path / 'native.tsv', tmp_path / 'first'
    design(source, first, native, '--max-mutations-in-report', '1', '--mrna-no-optimize-linkers')
    assert len(EpitopeDataset.load(native).epitopes) == 2
    with monkeypatch.context() as offline:
        offline.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
        offline.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
        replay(native, tmp_path / 'replay')
    assert product_files(first) == product_files(tmp_path / 'replay')
    replay(native, tmp_path / 'changed', '--max-mutations-in-report', '2')
    assert product_files(first) != product_files(tmp_path / 'changed')


def test_saved_mrna_only_design_restores_its_modality_and_ranking_settings(source, tmp_path, monkeypatch):
    native, first, second = tmp_path / 'native.tsv', tmp_path / 'first', tmp_path / 'replay'
    design(source, first, native, '--vaccine-type', 'mrna', '--mrna-no-optimize-linkers',
           '--mrna-poly-a-length', '73')
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    replay(native, second)
    assert not (second / 'peptide').exists()
    for name in ('cds.fasta', 'no_polyA.fasta', 'full.fasta', 'manifest.json', 'mrna-sequence-parts.csv'):
        assert (first / name).read_bytes() == (second / name).read_bytes()


def test_current_cli_reconciles_pooled_saved_construct_defaults(source, tmp_path):
    from .input_scope_helpers import write_input_manifest
    from topiary import read_tsv
    natives = [tmp_path / ('native%d.tsv' % i) for i in range(2)]
    for i, native in enumerate(natives):
        table = read_tsv(source)
        table.df['prediction_id'] = ['pool%d-%d' % (i, j) for j in range(len(table.df))]
        separate_source = tmp_path / ('source%d.tsv' % i)
        table.to_tsv(separate_source)
        design(separate_source, tmp_path / ('first%d' % i), native, '--mrna-no-optimize-linkers',
               '--peptide-max-constructs', str(i + 1))
    manifest = write_input_manifest(tmp_path / 'inputs.json', [('epitopes', path) for path in natives],
                                    patient_id='patient')
    argv = ['--input-manifest', manifest, '--output-dir', str(tmp_path / 'pool'),
            '--no-processing-aware-annotation']
    with pytest.raises(ValueError, match='different construct configurations'):
        run_cli(argv)
    run_cli(argv + ['--peptide-max-constructs', '20', '--output-epitopes', str(tmp_path / 'pool.tsv')])
    products = from_native_json(EpitopeDataset.load(tmp_path / 'pool.tsv').selection[
        'construct_assembly']['peptide']['products'], list)
    assert len(products) == 2


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('missing', [False, True])
def test_actual_products_replay_after_custom_annotation_moves(tmp_path, monkeypatch, strand, missing):
    from types import SimpleNamespace
    import shutil
    import pyensembl.download_cache
    from .test_native_references import native_dataset, delete_original_reference
    from .test_isovar_support_contract import reference
    genome = reference.__wrapped__(tmp_path, SimpleNamespace(param=strand))
    _, path = native_dataset(tmp_path, genome)
    first, native = tmp_path / 'first', tmp_path / 'design' / 'native.tsv'
    replay(path, first, '--vaccine-peptide-length', '10',
           '--mrna-no-optimize-linkers', '--output-epitopes', str(native))
    saved = EpitopeDataset.load(native)
    assert from_native_json(saved.selection['construct_assembly']['mrna']['products'], list)
    shutil.move(str(native.parent), str(tmp_path / 'relocated'))
    native = tmp_path / 'relocated' / native.name
    delete_original_reference(tmp_path)
    if missing:
        shutil.rmtree(str(native) + '.references')
        # Transcript-effect report annotation requires the bundle. Exercise
        # the real scoring, native and modality writers independently of it.
        monkeypatch.setattr('vaxrank.report.TemplateDataCreator.compute_template_data', lambda self: {})
    monkeypatch.setattr(pyensembl.download_cache.DownloadCache, '_fetch', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_peptide_constructs', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.assemble_mrna_constructs', forbidden)
    replay(native, tmp_path / 'replay', '--output-epitopes', str(tmp_path / 'again.tsv'))
    assert product_files(first) == product_files(tmp_path / 'replay')


def test_original_variant_prediction_captures_full_evidence_without_deferring_selection(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from .test_native_references import native_dataset
    from .test_isovar_support_contract import reference
    from vaxrank.construct_replay import pipeline_native_dataset
    from vaxrank.core_logic import vaccine_peptides_for_variant
    from vaxrank.vaccine_config import VaccineConfig
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    genome = reference.__wrapped__(tmp_path, SimpleNamespace(param='+'))
    evidence, _ = native_dataset(tmp_path, genome)
    fragment = next(iter(evidence.mutation_fragments.values()))
    vcf = tmp_path / 'source.vcf'
    vcf.write_text('synthetic VCF')
    bam = tmp_path / 'source.bam'
    bam.write_bytes(b'synthetic BAM')
    args = parse_vaxrank_args(['--vcf', str(vcf), '--bam', str(bam), '--mhc-alleles', 'HLA-A*02:01',
                              '--mhc-predictor', 'random', '--output-patient-id', 'patient'])
    captured = pipeline_native_dataset(args, [fragment.variant], evidence.config)
    monkeypatch.setattr('vaxrank.core_logic.MutantProteinFragment.from_isovar_result', lambda result: fragment)
    monkeypatch.setattr('vaxrank.core_logic.predict_epitopes', lambda **kwargs: evidence.epitopes)
    kwargs = dict(isovar_result=SimpleNamespace(variant=fragment.variant, passes_all_filters=True),
                  mhc_predictor=None, vaccine_peptide_length=10, epitope_config=evidence.config,
                  vaccine_config=VaccineConfig(preferred_peptide_length=10, min_peptide_length=10,
                                               max_peptide_length=10))
    original = vaccine_peptides_for_variant(**kwargs)
    selected = vaccine_peptides_for_variant(**kwargs, epitope_dataset=captured, defer_window_selection=False)
    assert selected and [p.amino_acids for p in selected] == [p.amino_acids for p in original]
    assert [p.target_epitope_score for p in selected] == [p.target_epitope_score for p in original]
    assert set(captured.mutation_fragments) == {e.prediction_id for p in selected for e in p.epitopes}
    assert captured.provenance[0].source_format == 'vcf_bam'
    assert captured.provenance[0].scope.patient_id == 'patient'
