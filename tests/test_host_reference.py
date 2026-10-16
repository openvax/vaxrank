"""Experimental host snapshots retain ORF sources through design and replay."""

from dataclasses import replace
import gzip
import json
import shutil
from types import SimpleNamespace

from orfdata import Proteome
from orfdata.proteome import Occurrence, ProteomeEntry, ViralProteome
from orfdata.reference import ReferenceContext
import pytest
from topiary import TopiaryPredictor

from vaxrank import HostReference
from vaxrank.cli.entry_point import run_cli
from vaxrank.construct_audit import audit_construct_sequence
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.external_rescoring import read_reports
from vaxrank.input_scope import InputScope, read_input_manifest
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.reference_proteome import self_reference_matches

from . import test_antigen_input
from .test_construct_audit import RecordingPredictor, jlf_construct, JLF
from .test_antigen_manifest import design_options
from .test_epitope_dataset import forbidden

inputs = test_antigen_input.inputs
dependencies = test_antigen_input.dependencies


def reference(*, sequence='SIINFEKLA', extra=True, species='human', assembly=None):
    taxon, name, default_assembly = ((9606, 'homo_sapiens', 'GRCh38') if species == 'human'
                             else (10090, 'mus_musculus', 'GRCm38'))
    context = ReferenceContext(taxon, name, assembly or default_assembly, 'ensembl', '93')
    def occurrence(transcript, gene):
        return Occurrence(transcript, gene, gene, 'protein_coding', '1', '+', 0,
                          len(sequence) * 3, ((1, len(sequence) * 3),), transcript == 'source')
    entries = [ProteomeEntry('protein', sequence, 'annotated', taxon, 'ensembl:93',
                             occurrences=(occurrence('annotated', 'CTA'),))]
    if extra:
        entries.append(ProteomeEntry('uORF', sequence, 'extra_orf', taxon, 'catalog:fixture',
            orf_type='uORF', evidence=('ribosome',), evidence_detail={'fixture': 'not biological evidence'},
            occurrences=(occurrence('source', 'CTA'), occurrence('other', 'OTHER')),
            source_transcript_id='source', source_binding='synthetic binding history'))
    layers = [dict(role='host', kind='annotated', source='ensembl:93', reference=context.key)]
    if extra:
        layers.append(dict(role='host', kind='extra_orf', source='catalog:fixture',
                           reference=context.key, evidence={'kind': 'ribosome'}))
    return Proteome(context, entries, layers)


def test_noncanonical_matches_preserve_all_occurrences_and_annotation_exclusions(inputs):
    host = HostReference(reference())
    antigen = replace(inputs.antigens['CTA'], self_reference_excluded_gene_ids=('CTA',))
    match, = host.matches(['SIINFEKLA'], antigen).values()
    assert match.occurs and match.source_provenance_complete
    assert match.reference_sha256 == host.identity
    assert [(s.gene_id, s.transcript_id, s.orf_id, s.reference_occurrence) for s in match.sources] == [
        ('CTA', 'source', 'uORF', 0), ('OTHER', 'other', 'uORF', 1)]
    assert from_native_json(to_native_json(match), type(match)) == match
    annotated = HostReference(reference(extra=False))
    assert not annotated.matches(['SIINFEKLA'], antigen)['SIINFEKLA'].occurs
    assert host.identity != annotated.identity


def test_viral_layers_never_count_as_host_self(inputs):
    proteome = reference().with_viral(ViralProteome(
        name='fixture virus', taxon_id=1234, host_taxon_ids=(9606,),
        proteins=(('viral', 'GILGFVFTL'),), source='synthetic', host_source='synthetic'))
    host = HostReference(proteome)
    assert any(e.kind == 'viral' for e in host.proteome)
    assert not host.matches(['GILGFVFTL'], inputs.antigens['viral'])['GILGFVFTL'].occurs


def test_snapshot_loading_pins_identity_and_detects_corruption(tmp_path):
    host = HostReference(reference())
    path = tmp_path / 'reference.json.gz'
    host.save(path)
    restored = HostReference.load({'path': path, 'identity': host.identity})
    assert restored.identity == host.identity
    assert restored.proteome.entries[1].evidence_detail == host.proteome.entries[1].evidence_detail
    assert restored.proteome.entries[1].source_binding == 'synthetic binding history'
    with pytest.raises(ValueError, match='declared identity'):
        HostReference.load({'path': path, 'identity': '0' * 64})
    payload = json.loads(gzip.decompress(path.read_bytes()))
    payload['content']['entries'][0]['sequence'] = 'GILGFVFTL'
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='identity'):
        HostReference.load(path)


@pytest.mark.parametrize('claim', [dict(host_taxon_id=10090), dict(reference_assembly='GRCm38'),
                                  dict(annotation='ensembl:94'), dict(host_reference_sha256='0' * 64)])
def test_reference_disagreement_fails_early(claim):
    with pytest.raises(ValueError, match='conflicting'):
        HostReference(reference()).validate(InputScope(**claim))


def test_candidate_and_final_product_audit_share_the_exact_reference(inputs):
    peptide = JLF[-9:]
    host = HostReference(reference(sequence=peptide))
    record = jlf_construct()
    candidate = self_reference_matches([peptide], record.assessment_antigen(), host_reference=host)[peptide]
    audit = audit_construct_sequence(record, RecordingPredictor(), host_reference=host)
    final = audit.mhc_assessment.ligands[-1].self_reference_match
    assert final == candidate
    assert final.reference_sha256 == host.identity
    assert {s.orf_id for s in final.sources} == {'', 'uORF'}


@pytest.mark.parametrize('snapshot', [False, True])
def test_final_product_reference_validates_species_and_allows_an_explicit_hla_mouse(monkeypatch, snapshot):
    host = HostReference(reference(species='mouse'))
    reference_args = dict(host_reference=host)
    if not snapshot:
        reference_args = dict(genome=SimpleNamespace(
            species='mus_musculus', reference_name='GRCm38', release=93))
        monkeypatch.setattr('vaxrank.safety_assessment.self_reference_matches',
            lambda peptides, antigen, genome: {p: replace(antigen.self_reference_match(p, False),
                genome_release='93') for p in peptides})
    with pytest.raises(ValueError, match='conflicts with host'):
        audit_construct_sequence(jlf_construct(), RecordingPredictor(), **reference_args)
    audit = audit_construct_sequence(jlf_construct(), RecordingPredictor(),
                                    **reference_args, mhc_species='human')
    assert audit.mhc_status == 'predictions_returned'
    with pytest.raises(ValueError, match='conflicts with mhc_species'):
        audit_construct_sequence(jlf_construct(), RecordingPredictor(),
                                 **reference_args, mhc_species='dog')


def test_manifest_resolves_and_validates_snapshot_declaration(tmp_path):
    path = tmp_path / 'inputs.json'
    path.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1',
        host_reference='host.json.gz', inputs=[dict(format='antigens', path='antigens.json')])))
    specs = read_input_manifest(path)
    assert specs[0][2]['host_reference'] == str(tmp_path / 'host.json.gz')
    document = json.loads(path.read_text())
    document['host_reference'] = {'path': 'host.json.gz', 'identity': '0' * 64}
    path.write_text(json.dumps(document))
    assert read_input_manifest(path)[0][2]['host_reference']['path'] == str(tmp_path / 'host.json.gz')


def test_cli_predicts_with_snapshot_and_replays_after_relocation(inputs, dependencies, tmp_path, monkeypatch):
    _, _ = dependencies
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', None))
    source = tmp_path / 'source'
    source.mkdir()
    host = HostReference(reference())
    host.save(source / 'host.json.gz')
    replace(inputs, antigens={'viral': inputs.antigens['viral']}).save(source / 'antigens.json')
    manifest = source / 'inputs.json'
    manifest.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1',
        host_reference='host.json.gz', inputs=[dict(format='antigens', path='antigens.json')])))
    original = tmp_path / 'original'
    original.mkdir()
    saved = original / 'native.tsv'
    run_cli(['--input-manifest', str(manifest), '--mhc-predictor', 'random',
             '--output-epitopes', str(saved), '--no-processing-aware-annotation'])
    dataset = EpitopeDataset.load(saved)
    assert dataset.host_reference.identity == host.identity
    assert all(p.scope.host_reference_sha256 == host.identity for p in dataset.provenance)
    assert all(e.self_reference_match.reference_sha256 == host.identity for e in dataset.epitopes)
    moved = tmp_path / 'moved'
    shutil.copytree(original, moved)
    shutil.rmtree(original)
    shutil.rmtree(source)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('pyensembl.Genome.transcripts', forbidden)
    replay = moved / 'replayed.tsv'
    run_cli(['--input-epitopes', str(moved / 'native.tsv'), '--output-epitopes', str(replay),
             '--no-processing-aware-annotation'])
    restored = EpitopeDataset.load(replay)
    assert restored.host_reference.identity == host.identity
    assert restored.provenance == dataset.provenance
    assert restored.epitopes == dataset.epitopes
    resources = list(moved.glob('replayed.tsv.references/*'))
    assert len(resources) == 1
    resources[0].write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='corrupt'):
        read_reports([('epitopes', replay)])


@pytest.mark.parametrize('reverse', [False, True])
def test_mixed_inputs_inherit_the_saved_snapshot_before_prediction(
        inputs, dependencies, tmp_path, reverse):
    model, _ = dependencies
    host = HostReference(reference())
    predictor = TopiaryPredictor(models=[model])
    saved = tmp_path / 'saved.tsv'
    replace(inputs, antigens={'viral': inputs.antigens['viral']}).predict(
        predictor, host_reference=host).save(saved)
    new = tmp_path / 'new.json'
    replace(inputs, antigens={'mutation': inputs.antigens['mutation']}).save(new)
    entries = [('antigens', new), ('epitopes', saved)]
    if reverse:
        entries.reverse()
    reports = read_reports(entries, antigen_predictor_factory=lambda scope: predictor)
    assert all(report.dataset.host_reference.identity == host.identity for report in reports)
    assert all(report.input_provenance.scope.host_reference_sha256 == host.identity for report in reports)


@pytest.mark.parametrize('snapshot', [False, True])
@pytest.mark.parametrize('saved', [False, True])
@pytest.mark.parametrize('reverse', [False, True])
def test_new_sequences_inherit_pooled_mhc_species_before_validation(
        inputs, dependencies, tmp_path, snapshot, saved, reverse):
    model, genome = dependencies
    genome.species, genome.reference_name = 'mus_musculus', 'GRCm38'
    genome.annotation_name, genome.annotation_version = 'ensembl', '93'
    reference_args = (dict(host_reference=reference(species='mouse')) if snapshot else dict(genome=genome))
    predictor = TopiaryPredictor(models=[model])
    entries = []
    for index, mhc_species in enumerate((None, 'human')):
        source = replace(inputs, antigens={str(index): inputs.antigens['viral']},
            scope=replace(inputs.scope, host_taxon_id=10090, reference_assembly='GRCm38',
                          mhc_species=mhc_species, sample_id=str(index)))
        fmt = 'epitopes' if saved and index else 'antigens'
        path = tmp_path / (str(index) + ('.tsv' if fmt == 'epitopes' else '.json'))
        if fmt == 'epitopes':
            source.predict(predictor, **reference_args).save(path)
        else:
            source.save(path)
        entries.append((fmt, path))
    if reverse:
        entries.reverse()
    reports = read_reports(entries, **reference_args, antigen_predictor_factory=lambda scope: predictor)
    for index, report in enumerate(reports):
        assert report.input_provenance.scope.mhc_species == ('Homo sapiens',)
        output = tmp_path / f'export-{index}.tsv'
        report.dataset.save(output)
        assert read_reports([('epitopes', output)])[0].input_provenance.scope.mhc_species == ('Homo sapiens',)


def test_snapshot_reference_supports_custom_assemblies_through_prediction_and_junctions(inputs, dependencies):
    from vaxrank.input_scope import validate_prediction_species
    model, _ = dependencies
    host = HostReference(reference(assembly='GCF_000001405.40'))
    source = replace(inputs, antigens={'viral': inputs.antigens['viral']},
                     scope=replace(inputs.scope, reference_assembly=host.scope.reference_assembly))
    dataset = source.predict(TopiaryPredictor(models=[model]), host_reference=host)
    args = SimpleNamespace(_design_dataset=dataset, _input_species_scope={})
    scope = validate_prediction_species(args, [], dataset.provenance[0].scope.mhc_alleles, 'Junction')
    assert scope.host_reference_sha256 == host.identity


def test_snapshot_supplies_species_metadata_without_an_alias_registry_entry(inputs, dependencies, tmp_path):
    from vaxrank.external_rescoring import prepare_reports
    model, _ = dependencies
    base = reference()
    context = ReferenceContext(9544, 'macaca_mulatta', 'Mmul_10', 'refseq', '103')
    host = HostReference(Proteome(context, [replace(e, taxon_id=9544) for e in base.entries]))
    model.alleles = ['Mamu-A*01:01']
    source = replace(inputs, antigens={'viral': inputs.antigens['viral']}, scope=replace(inputs.scope,
        reference_assembly='Mmul_10', annotation='refseq:103', mhc_alleles=tuple(model.alleles)))
    dataset = source.predict(TopiaryPredictor(models=[model]), host_reference=host)
    saved = tmp_path / 'macaque.tsv'
    dataset.save(saved)
    _, frame, _ = prepare_reports([('epitopes', saved)])
    assert frame.attrs['epitope_dataset'].host_reference.identity == host.identity
    predictor = RecordingPredictor()
    original = predictor.predict_from_named_sequences
    predictor.predict_from_named_sequences = lambda sequences: original(sequences).assign(allele='Mamu-A*01:01')
    assert audit_construct_sequence(jlf_construct(), predictor, host_reference=host).mhc_status == 'predictions_returned'


def test_conflicting_snapshots_fail_before_predictor(inputs, dependencies, tmp_path):
    model, _ = dependencies
    paths = []
    for index, extra in enumerate((False, True)):
        path = tmp_path / f'{index}.tsv'
        replace(inputs, antigens={'viral': inputs.antigens['viral']}).predict(
            TopiaryPredictor(models=[model]), host_reference=reference(extra=extra)).save(path)
        paths.append(('epitopes', path))
    with pytest.raises(ValueError, match='conflicting host_reference'):
        read_reports(paths, antigen_predictor_factory=forbidden)


def test_host_snapshot_rejects_wrong_species_before_predictor(inputs, tmp_path):
    path = tmp_path / 'input.json'
    inputs.save(path)
    with pytest.raises(ValueError, match='conflicting reference_assembly'):
        read_reports([('antigens', path)], host_reference=reference(species='mouse'),
                     antigen_predictor_factory=forbidden)


def test_unpinned_prediction_evidence_cannot_change_self_reference_on_reload(inputs, dependencies, tmp_path):
    model, genome = dependencies
    saved = tmp_path / 'old.tsv'
    inputs.predict(TopiaryPredictor(models=[model]), genome=genome).save(saved)
    with pytest.raises(ValueError, match='saved predictions have no host_reference'):
        read_reports([('epitopes', saved)], host_reference=reference())


@pytest.mark.parametrize('assembly', ['GRCh38', 'custom-human-assembly'])
@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_product_design_with_snapshot_replays_offline(
        inputs, dependencies, tmp_path, monkeypatch, modality, assembly):
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', None))
    source = tmp_path / 'viral.json'
    replace(inputs, antigens={'viral': inputs.antigens['viral']},
            scope=replace(inputs.scope, reference_assembly=assembly)).save(source)
    snapshot = tmp_path / 'host.json.gz'
    HostReference(reference(assembly=assembly)).save(snapshot)
    original = tmp_path / 'original'
    run_cli(['--input-antigens', str(source), '--host-reference', str(snapshot),
             '--mhc-predictor', 'random'] + design_options(original, modality))
    saved = original / 'candidate_predictions.tsv'
    assert EpitopeDataset.load(saved).host_reference is not None
    source.unlink()
    snapshot.unlink()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('pyensembl.Genome.transcripts', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(saved), '--output-dir', str(replay),
             '--no-processing-aware-annotation', '--pdf-backend', 'weasyprint'])
    assert (replay / 'manifest.json').read_bytes() == (original / 'manifest.json').read_bytes()
    for fasta in original.glob('*.fasta'):
        assert (replay / fasta.name).read_bytes() == fasta.read_bytes()


def test_reference_identity_invalidates_assembly_even_with_identical_sequences(tmp_path):
    from types import SimpleNamespace
    from vaxrank.construct_replay import _assembly_identity
    from vaxrank.peptide import PeptideConstructConfig
    from vaxrank.native_references import NativeReferences

    host = reference()
    changed = Proteome(host.reference, host.entries, [
        *host.layers, dict(role='filter', description='different explicit selection', removed=0)])
    dataset = SimpleNamespace(host_reference=HostReference(host),
                              _assembly_references=NativeReferences(tmp_path / 'native.tsv'))
    args = SimpleNamespace(_design_dataset=dataset, _construct_policy={})
    before = _assembly_identity(args, [], PeptideConstructConfig(), None)
    dataset.host_reference = HostReference(changed)
    assert _assembly_identity(args, [], PeptideConstructConfig(), None) != before
