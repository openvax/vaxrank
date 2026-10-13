"""Unified antigen evidence, native occurrence prediction and product replay."""

from dataclasses import replace
import hashlib
import json
from types import SimpleNamespace

import pytest
from topiary import TopiaryPredictor

from vaxrank import AntigenDataset
from vaxrank.antigen_input import ANTIGEN_DATASET_SCHEMA
from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.input_scope import InputScope
from vaxrank.vaccine_antigen import (
    ANTIGEN_KINDS, AminoAcidInterval, TargetableMask,
    TumorSpecificityAttestation, TumorSpecificityEvidence, VaccineAntigen,
)

from .test_epitope_dataset import forbidden
from .test_external_rescoring import ALLELES, ContextPredictor
from . import test_isovar_support_contract

reference = test_isovar_support_contract.reference


@pytest.fixture
def inputs():
    records = {}
    for index, kind in enumerate(sorted(ANTIGEN_KINDS)):
        sequence = 'ACDEFGHIKLMNPQRSTVWY'[index:index + 10]
        records[kind] = VaccineAntigen(
            kind=kind, amino_acids=sequence, source_identifier='synthetic-' + kind,
            gene_id='fixture-gene-' + kind, gene_name='fixture-' + kind,
            transcript_ids=('transcript-' + kind,), protein_ids=('protein-' + kind,),
            species='fixture virus' if kind == 'viral' else 'Homo sapiens',
            source_metadata=(('provenance', 'synthetic integration fixture'),),
            targetable_mask=TargetableMask((AminoAcidInterval(4, 5),)),
            tumor_specificity=TumorSpecificityAttestation(
                status='admitted', evidence_kind='synthetic_fixture', evidence_source=__name__,
                patient_specific=True, evidence_records=(TumorSpecificityEvidence(
                    evidence_kind='synthetic_fixture', evidence_source=__name__,
                    patient_specific=True, sample_id='sample', passed=True),)))
    records['held'] = replace(records['viral'], source_identifier='held-virus',
                              tumor_specificity=replace(records['viral'].tumor_specificity, status='held_out'))
    return AntigenDataset(records, InputScope(patient_id='patient', sample_id='sample',
        reference_assembly='GRCh38', annotation='ensembl:93', mhc_alleles=tuple(ALLELES)))


@pytest.fixture
def dependencies(monkeypatch):
    model = ContextPredictor()
    genome = SimpleNamespace(reference_name='GRCh38', release=93)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', lambda args: [model])
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', genome))
    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches',
        lambda peptides, antigen, genome: {
            p: replace(antigen.self_reference_match(p, False), source_provenance_complete=True)
            for p in peptides})
    return model, genome


def test_antigen_json_roundtrip_and_shared_prediction(inputs, tmp_path, dependencies):
    path = tmp_path / 'antigens.json'
    inputs.save(path)
    document = json.loads(path.read_text())
    assert document['schema'] == ANTIGEN_DATASET_SCHEMA
    assert document['antigens']['viral']['species'] == 'fixture virus'
    assert '__class__' not in path.read_text()
    assert AntigenDataset.load(path) == inputs
    model, genome = dependencies
    dataset = inputs.predict(TopiaryPredictor(models=[model]), genome=genome,
                             config=EpitopeConfig(score_expr='1.0'))
    assert model.calls
    assert {a.kind for a in dataset.antigens.values()} == ANTIGEN_KINDS
    assert len(dataset.epitopes) == 14  # two occurrences per source, including held-out evidence
    assert len({e.prediction_id for e in dataset.epitopes}) == 14
    assert all(e.self_reference_match.source_provenance_complete for e in dataset.epitopes)
    assert all(e.source_sequence[e.offset:e.offset + len(e.sequence)] == e.sequence for e in dataset.epitopes)
    assert all(e.patient_alleles == tuple(ALLELES) for e in dataset.epitopes)
    exported = AntigenDataset.from_epitope_dataset(dataset)
    assert len(exported.antigens) == 7
    assert exported.antigens == inputs.antigens
    assert set(exported.antigens.values()) == set(inputs.antigens.values())
    exported.save(path)
    assert AntigenDataset.load(path) == exported


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_all_sources_use_shared_design_and_offline_replay(inputs, dependencies, tmp_path, monkeypatch, modality):
    path = tmp_path / 'input.json'
    inputs.save(path)
    output = tmp_path / 'design'
    native = output / 'candidate_predictions.tsv'
    exported = tmp_path / 'export.json'
    options = ['--input-antigens', str(path), '--mhc-predictor', 'random',
        '--ensembl-release', '93', '--output-dir', str(output),
        '--output-antigens', str(exported), '--vaccine-type', modality,
        '--no-processing-aware-annotation', '--pdf-backend', 'weasyprint',
        '--vaccine-peptide-length', '10', '--config-text', 'epitopes.score_expr=1.0',
        '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score']
    if modality == 'mrna':
        options += ['--mrna-antigens-per-construct', '6', '--mrna-max-constructs', '1',
                    '--mrna-no-optimize-linkers', '--mrna-codon-method', 'use_best_codon']
    run_cli(options)
    assert len(AntigenDataset.load(exported).antigens) == 7
    saved = EpitopeDataset.load(native)
    assert {a.kind for a in saved.antigens.values()} == ANTIGEN_KINDS
    manifest = json.loads((output / 'manifest.json').read_text())
    sequences = ([p['sequence'] for p in manifest] if modality == 'peptide' else
                 [a['aa'] for p in manifest for a in p['antigens']])
    assert len(sequences) == 6
    assert set(sequences) == {a.amino_acids for a in inputs.antigens.values()}
    path.unlink()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(native), '--output-dir', str(replay),
             '--pdf-backend', 'weasyprint', '--no-processing-aware-annotation'])
    for fasta in output.glob('*.fasta'):
        assert fasta.read_bytes() == (replay / fasta.name).read_bytes()
    assert (output / 'manifest.json').read_bytes() == (replay / 'manifest.json').read_bytes()


@pytest.mark.parametrize('field,value', [('reference_assembly', 'GRCh37'),
                                       ('annotation', 'ensembl:75')])
def test_conflicting_host_reference_fails_before_prediction(inputs, dependencies, field, value):
    model, genome = dependencies
    inputs.scope = replace(inputs.scope, **{field: value})
    with pytest.raises(ValueError, match='conflicting ' + field):
        inputs.predict(TopiaryPredictor(models=[model]), genome=genome)
    assert not model.calls


def test_cli_genotype_conflict_fails_before_models(inputs, dependencies, tmp_path, monkeypatch):
    path = tmp_path / 'input.json'
    inputs.save(path)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='conflicting mhc_alleles'):
        run_cli(['--input-antigens', str(path), '--mhc-predictor', 'random',
                 '--ensembl-release', '93', '--mhc-alleles', 'HLA-A*01:01',
                 '--output-epitopes', str(tmp_path / 'out.tsv')])


def test_unmeasured_and_held_out_antigens_survive_prediction(inputs, dependencies, tmp_path):
    inputs.antigens['short'] = replace(inputs.antigens['held'], amino_acids='ACD',
        source_identifier='unmeasured-short', targetable_mask=TargetableMask())
    model, genome = dependencies
    result = inputs.predict(TopiaryPredictor(models=[model]), genome=genome)
    path = tmp_path / 'predictions.tsv'
    result.save(path)
    exported = AntigenDataset.from_epitope_dataset(EpitopeDataset.load(path))
    assert set(exported.antigens.values()) == set(inputs.antigens.values())


@pytest.mark.parametrize('change,match', [
    (lambda d: d.update(schema='future'), 'Expected antigen schema'),
    (lambda d: d.update(unknown=True), 'Unknown antigen'),
    (lambda d: d['antigens']['viral'].update(kind='guessed'), 'Unknown antigen kind'),
    (lambda d: d['antigens']['viral']['targetable_mask']['intervals'][0].update(end=999), 'beyond antigen'),
    (lambda d: d['antigens']['viral']['tumor_specificity'].update(status='overridden'), 'override reason'),
    (lambda d: d['antigens']['viral']['tumor_specificity'].update(patient_specific='false'), 'Expected `bool`'),
    (lambda d: d['antigens']['viral']['targetable_mask']['intervals'][0].update(start=1.5), 'Expected `int`'),
    (lambda d: d['antigens']['viral']['tumor_specificity'].update(typo='ignored'), 'Unknown TumorSpecificityAttestation'),
])
def test_malformed_antigen_files_fail(inputs, tmp_path, change, match):
    # JSON roundtrip gives mutable nested values as supplied by external callers.
    document = json.loads(json.dumps(inputs.to_dict()))
    change(document)
    path = tmp_path / 'bad.json'
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match=match):
        AntigenDataset.load(path)


def test_same_sequence_keeps_distinct_record_identifiers(inputs, dependencies, tmp_path):
    inputs.antigens['second-viral-occurrence'] = inputs.antigens['viral']
    model, genome = dependencies
    result = inputs.predict(TopiaryPredictor(models=[model]), genome=genome)
    original = [e for e in result.epitopes if e.source_name == 'viral']
    second = [e for e in result.epitopes if e.source_name == 'second-viral-occurrence']
    assert len(original) == len(second) == 2
    assert {e.prediction_id for e in original}.isdisjoint(e.prediction_id for e in second)
    path = tmp_path / 'predictions.tsv'
    result.save(path)
    assert AntigenDataset.from_epitope_dataset(EpitopeDataset.load(path)).antigens == inputs.antigens


def test_host_identity_and_antigen_self_exclusions_stay_distinct(inputs, dependencies, monkeypatch):
    inputs.antigens = {'CTA': replace(inputs.antigens['CTA'], self_reference_excluded_gene_ids=('fixture-gene-CTA',))}
    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches',
        lambda peptides, antigen, genome: {
            p: antigen.self_reference_match(p, not antigen.self_reference_excluded_gene_ids)
            for p in peptides})
    model, genome = dependencies
    result = inputs.predict(TopiaryPredictor(models=[model]), genome=genome)
    assert all(e.occurs_in_reference and not e.self_reference_match.occurs for e in result.epitopes)
    assert result.result.df.occurs_in_reference.all()
    assert not result.result.df.occurs_in_non_CTA_reference.any()


@pytest.mark.parametrize('genotype_source', ['input', 'flag', 'file'])
def test_real_mhctools_cli_and_host_reference(inputs, reference, tmp_path, monkeypatch, genotype_source):
    # Use the installed MHCtools predictor factory and actual PyEnsembl indexes.
    # Only annotation resolution is redirected to the small local fixture.
    _, variant = reference
    genome = variant.ensembl
    genome.reference_name = 'GRCh38'
    inputs.scope = replace(inputs.scope, annotation='test:1')
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', genome))
    path, output = tmp_path / 'input.json', tmp_path / 'predictions.tsv'
    inputs.save(path)
    options = ['--input-antigens', str(path), '--mhc-predictor', 'random',
               '--ensembl-release', '93', '--output-epitopes', str(output),
               '--config-text', 'epitopes.score_expr=1.0', '--no-processing-aware-annotation']
    if genotype_source == 'flag':
        options += ['--mhc-alleles', ','.join(ALLELES)]
    elif genotype_source == 'file':
        allele_file = tmp_path / 'alleles.txt'
        allele_file.write_text('\n'.join(ALLELES))
        options += ['--mhc-alleles-file', str(allele_file)]
    run_cli(options)
    result = EpitopeDataset.load(output)
    assert result.epitopes
    assert {a.kind for a in result.antigens.values()} == ANTIGEN_KINDS
    assert all(e.self_reference_match.source_provenance_complete for e in result.epitopes)
    assert set(result.result.df.allele) == set(ALLELES)
    assert result.provenance[0].content_sha256 == hashlib.sha256(path.read_bytes()).hexdigest()


def test_pooled_native_sources_keep_independent_antigen_catalogs(inputs, dependencies, tmp_path):
    from vaxrank.external_rescoring import prepare_reports
    model, genome = dependencies
    paths = []
    for sample in ['sample-1', 'sample-2']:
        inputs.scope = replace(inputs.scope, sample_id=sample)
        dataset = inputs.predict(TopiaryPredictor(models=[model]), genome=genome,
                                 config=EpitopeConfig(score_expr='1.0'))
        path = tmp_path / (sample + '.tsv')
        dataset.save(path)
        paths.append(('epitopes', str(path)))
    _, frame, _ = prepare_reports(paths, genome=genome)
    exported = AntigenDataset.from_epitope_dataset(frame.attrs['epitope_dataset'])
    assert len(exported.antigens) == 14
    assert {p.scope.sample_id for p in exported.provenance} == {'sample-1', 'sample-2'}
    path = tmp_path / 'pooled.json'
    exported.save(path)
    assert AntigenDataset.load(path) == exported


def test_native_mutation_export_retains_target_coordinates_and_provenance(reference, tmp_path):
    from .test_native_references import native_dataset
    dataset, _ = native_dataset(tmp_path, reference)
    exported = AntigenDataset.from_epitope_dataset(dataset)
    assert set(exported.antigens.values()) == set(dataset.antigens.values())
    assert all(a.kind == 'mutation' and a.targetable_mask.intervals for a in exported.antigens.values())
    path = tmp_path / 'mutation.json'
    exported.save(path)
    assert AntigenDataset.load(path) == exported


def test_selected_annotation_is_recorded_even_if_absent_from_input(inputs, dependencies):
    model, genome = dependencies
    inputs.scope = replace(inputs.scope, annotation=None)
    first = inputs.predict(TopiaryPredictor(models=[model]), genome=genome)
    genome.release = 94
    second = inputs.predict(TopiaryPredictor(models=[model]), genome=genome)
    assert first.provenance[0].scope.annotation == 'ensembl:93'
    assert second.provenance[0].scope.annotation == 'ensembl:94'
    assert first.provenance[0].input_id != second.provenance[0].input_id
    assert inputs.scope.annotation is None
