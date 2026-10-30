"""Released canine evidence with synthetic sequences, never a validated panel."""

from dataclasses import replace
import hashlib
import json
import shutil
from types import SimpleNamespace

import oncoref
from orfdata import Occurrence, ProteinSequence, Proteome
from orfdata.reference import ReferenceContext
from packaging.version import Version
import pandas as pd
import pytest

from vaxrank.cta_admission import CTAAdmissionPolicy, CTAOverrideEvidence
from vaxrank.cta_expression import CTAExpressionResult, admit_cta_expression
from vaxrank.cta_identity import identity_digest
from vaxrank.host_reference import HostReference
from vaxrank.species_cta import SpeciesCTAExpressionResult
from vaxrank.native_serialization import to_native_json
from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.input_scope import InputScope

from . import test_antigen_input

from .test_cta_expression import contract, SEQUENCE
from .test_cta_manifest import options
from .test_epitope_dataset import forbidden

PRAME = '477552'
CTCFL = '477271'
REVIEW = CTAOverrideEvidence('synthetic review fixture', '1', 'Exercise explicit research selection', 'fixture')
inputs = test_antigen_input.inputs


def host(*, shared=False, extra=False):
    snapshot = oncoref.canine.cta_snapshot()
    identity = snapshot.as_dict()['references'][0]['identity']
    reference = ReferenceContext(
        9615, 'canis_lupus_familiaris', 'UU_Cfam_GSD_1.0', 'NCBI RefSeq', '106',
        assembly_accession=identity['assembly_accession'],
        source_sha256=identity_digest({k: identity['content_provenance'][k] for k in ('protein', 'feature_table')}))
    entries = []
    for candidate in snapshot.candidates():
        for target in candidate['mapping']['targets']:
            if target['id'] not in (PRAME, CTCFL, '403710', '492173'):
                continue
            for pair in target['sequences']:
                if any(e.accession == pair['protein_id'] for e in entries):
                    continue
                occurrence = Occurrence(pair['transcript_id'], target['id'], target['symbol'],
                    'protein_coding', '1', '+', 0, len(SEQUENCE) * 3, ((1, len(SEQUENCE) * 3),), True)
                entries.append(ProteinSequence(pair['protein_id'], SEQUENCE, 'annotated', 9615,
                                              'synthetic fixture', occurrences=(occurrence,)))
    if shared or extra:
        occurrence = replace(entries[0].occurrences[0], gene_id='OTHER', transcript_id='other')
        entries.append(ProteinSequence('other', SEQUENCE, 'extra_orf' if extra else 'annotated',
                                     9615, 'synthetic fixture', occurrences=(occurrence,)))
    layer = dict(role='host', kind='annotated', source='synthetic fixture',
        content_provenance={k: identity['content_provenance'][k] for k in ('protein', 'feature_table')},
        identifier_namespaces=identity['identifier_namespaces'], added=len(entries))
    reference = replace(reference, sequence_sha256=identity_digest(sorted((e.accession, e.sequence)
        for e in entries if e.kind == 'annotated')))
    layer['added'] = sum(e.kind == 'annotated' for e in entries)
    return HostReference(Proteome(reference, entries, (layer,)))


def admit(tmp_path, rows=((PRAME, 10),), *, reference=None, level='gene', **kwargs):
    path = tmp_path / 'expression.tsv'
    pd.DataFrame(rows, columns=['feature_id', 'patient_tpm']).to_csv(path, sep='\t', index=False)
    return admit_cta_expression(path, input_contract=contract(level),
        admission_policy=CTAAdmissionPolicy(2, 'TPM'), host_reference=reference or host(), **kwargs)


def test_species_cta_requires_review_and_preserves_counterevidence(tmp_path):
    result = admit(tmp_path, [(PRAME, 10), (CTCFL, 10), ('999', 1)])
    assert isinstance(result, SpeciesCTAExpressionResult)
    assert [d.status for d in result.decisions] == ['held_out', 'held_out', 'not_cta']
    assert result.admitted_antigens == ()
    assert result.evidence_snapshot['coverage']['status'] == 'partial'
    assert any(e['tissue'] == 'kidney' for e in result.evidence_snapshot['evidence'])
    assert all(not d.assessment.antigen.self_reference_excluded_gene_ids
               for d in result.decisions if d.assessment)
    assert result.input_contract.id_namespace == 'ncbi_gene_id'
    ctcfl = result.decisions[1].assessment.antigen
    captured = dict(ctcfl.tumor_specificity.evidence_records[0].details)
    assert any(row['tissue'] == 'kidney' for row in json.loads(captured['observations']))
    assert json.loads(captured['coverage'])['status'] == 'partial'


def test_species_cta_review_does_not_bypass_expression_or_hide_other_genes(tmp_path):
    reference = host(shared=True, extra=True)
    result = admit(tmp_path, [(PRAME, 10), (CTCFL, 0)], reference=reference,
                   overrides={PRAME: REVIEW, CTCFL: REVIEW})
    assert [d.status for d in result.decisions] == ['admitted', 'held_out']
    antigen, = result.admitted_antigens
    assert antigen.tumor_specificity.status == 'overridden'
    assert antigen.self_reference_excluded_gene_ids == (PRAME,)
    match = reference.matches([SEQUENCE[:9]], antigen)[SEQUENCE[:9]]
    assert {source.gene_id for source in match.sources} >= {CTCFL, 'OTHER'}
    assert PRAME not in {source.gene_id for source in match.sources if not source.orf_id}


def test_species_cta_transcript_measurement_selects_only_exact_version(tmp_path):
    reference = host()
    transcript = next(e.transcript_id for e in reference.proteome if e.gene_id == PRAME)
    result = admit(tmp_path, [(transcript, 7)], reference=reference, level='transcript', overrides={PRAME: REVIEW})
    assert result.admitted_antigens[0].transcript_ids == (transcript,)
    assert result.decisions[0].assessment.tumor_expression.transcript_id == transcript
    assert result.input_contract.id_namespace == 'refseq_transcript_id'
    with pytest.raises(ValueError, match='versioned RefSeq'):
        admit(tmp_path, [(transcript.split('.')[0], 7)], reference=reference, level='transcript')


@pytest.mark.parametrize('change', [dict(taxon_id=9606, species='homo_sapiens'),
    dict(assembly_accession='GCF_000001405.40'), dict(annotation_version='107'),
    dict(annotation_name='ensembl'), dict(source_sha256='0' * 64)])
def test_species_cta_wrong_reference_fails_before_prediction(tmp_path, change):
    original = host()
    context = replace(original.proteome.reference, **change)
    entries = [replace(e, taxon_id=context.taxon_id) for e in original.proteome]
    layers = [{k: v for k, v in layer.items() if k != 'reference'} for layer in original.proteome.layers]
    reference = HostReference(Proteome(context, entries, layers))
    with pytest.raises(ValueError, match='annotation|SHA-256|genome|host annotation'):
        admit(tmp_path, reference=reference)


def test_species_cta_requires_both_source_assets_and_their_id_namespaces(tmp_path):
    original = host()
    layer = dict(original.proteome.layers[0])
    for change in (dict(content_provenance={'protein': '0' * 64}), dict(identifier_namespaces={})):
        with pytest.raises(ValueError, match='source SHA-256'):
            admit(tmp_path, reference=HostReference(Proteome(original.proteome.reference,
                original.proteome.proteins, ({**layer, **change},))))


def test_species_cta_rejects_a_filtered_annotated_self_inventory(tmp_path):
    original = host()
    filtered = original.proteome.filter(lambda e: e.gene_id == PRAME, description='target only')
    with pytest.raises(ValueError, match='complete annotated protein inventory'):
        admit(tmp_path, reference=HostReference(filtered))


@pytest.mark.parametrize('value', ['477552.0', 'PRAME', 'ENSG00000185686', '', None])
def test_species_cta_invalid_gene_ids_fail_on_ingest(tmp_path, value):
    with pytest.raises(ValueError, match='NCBI Gene'):
        admit(tmp_path, [(value, 10)])


def test_species_cta_missing_expression_and_ineffective_overrides(tmp_path):
    result = admit(tmp_path, [(PRAME, None)])
    assert result.decisions[0].status == 'unknown_expression'
    assert admit(tmp_path, [(PRAME, None)], overrides={PRAME: REVIEW}).decisions[0].status == 'unknown_expression'
    with pytest.raises(ValueError, match='resolved expression targets'):
        admit(tmp_path, overrides={'999': REVIEW})
    with pytest.raises(ValueError, match='Duplicate normalized'):
        admit(tmp_path, [(PRAME, 3), (PRAME, 4)])


def test_species_cta_unresolved_residues_hold_the_selected_isoform_and_keep_other_targets(tmp_path):
    original = host()
    longest = next(e for e in original.proteome if e.gene_id == PRAME)
    uncertain = replace(longest, sequence=longest.sequence + 'X')
    entries = [uncertain if e.accession == longest.accession else e for e in original.proteome]
    context = replace(original.proteome.reference,
        sequence_sha256=identity_digest(sorted((e.accession, e.sequence) for e in entries)))
    reference = HostReference(Proteome(context, entries, original.proteome.layers))
    result = admit(tmp_path, [(PRAME, 10), (CTCFL, 10)], reference=reference,
                   overrides={PRAME: REVIEW, CTCFL: REVIEW})
    unresolved, valid = result.decisions
    assert unresolved.status == 'unsupported_protein_sequence'
    assert unresolved.assessment is None
    assert unresolved.transcript_id == longest.transcript_id
    assert unresolved.protein_id == longest.accession
    assert unresolved.protein_sequence_sha256 == hashlib.sha256(uncertain.sequence.encode()).hexdigest()
    assert unresolved.unsupported_amino_acids == 'X'
    assert valid.status == 'admitted'
    assert {a.gene_id for a in result.admitted_antigens} == {CTCFL}
    saved = tmp_path / 'unresolved.json'
    result.save(saved)
    assert CTAExpressionResult.load(saved) == result
    with pytest.raises(ValueError, match='frozen mapping'):
        replace(result, decisions=(replace(unresolved, protein_id='XP_999.1'), valid))


def test_species_cta_saved_result_replays_without_current_reference(tmp_path, monkeypatch):
    result = admit(tmp_path, overrides={PRAME: REVIEW})
    path = tmp_path / 'admission.json'
    result.save(path)
    monkeypatch.setattr('oncoref.species_antigens.cta_snapshot', forbidden)
    monkeypatch.setattr(oncoref, '__version__', 'future')
    assert CTAExpressionResult.load(path) == result
    decision = result.decisions[0]
    with pytest.raises(ValueError, match='frozen mapping or expression'):
        replace(result, decisions=(replace(decision, assessment=replace(decision.assessment,
            antigen=replace(decision.assessment.antigen, amino_acids=SEQUENCE[::-1]))),))
    with pytest.raises(ValueError, match='annotation disagree'):
        replace(result, annotation_version='107')
    with pytest.raises(ValueError, match='frozen mapping or expression'):
        replace(result, decisions=(replace(decision, gene_identity={}),))
    with pytest.raises(ValueError, match='frozen mapping or expression'):
        replace(result, input_contract=replace(result.input_contract, sample_id='another-patient'))
    with pytest.raises(ValueError, match='lost or changed'):
        antigen = decision.assessment.antigen
        replace(result, decisions=(replace(decision, assessment=replace(decision.assessment,
            antigen=replace(antigen, tumor_specificity=replace(antigen.tumor_specificity,
                evidence_records=antigen.tumor_specificity.evidence_records[1:])))),))


def test_species_cta_unavailable_evidence_is_explicit(tmp_path, monkeypatch):
    # Exercise unavailable evidence independently of additions to the owner's
    # published species inventory.
    snapshot = SimpleNamespace(taxon_id=10090, coverage={'status': 'unavailable'})
    monkeypatch.setattr('oncoref.species_antigens.cta_snapshot', lambda *args, **kwargs: snapshot)
    reference = HostReference(Proteome(ReferenceContext(10090, 'mouse', 'GRCm38', 'ensembl', '93'), []))
    with pytest.raises(ValueError, match='evidence is unavailable'):
        admit(tmp_path, reference=reference)


@pytest.mark.skipif(Version(oncoref.__version__) < Version('1.22.1'),
                    reason='Native mouse P1A evidence was added in OncoRef 1.22.1')
def test_species_cta_partial_mouse_reference_remains_unsupported(tmp_path):
    snapshot = oncoref.mouse.cta_snapshot()
    assert snapshot.coverage['status'] == 'partial'
    identity = snapshot.as_dict()['references'][0]['identity']
    context = ReferenceContext(10090, 'mouse', 'GRCm39', identity['annotation_source'],
        identity['annotation_release'], assembly_accession=identity['assembly_accession'])
    # A native nomination and a single-gene reference subset cannot substitute
    # for the supported expression-mapping and complete-host inventory contract.
    with pytest.raises(ValueError, match='candidate subsets are insufficient'):
        admit(tmp_path, reference=HostReference(Proteome(context, [])))


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
@pytest.mark.parametrize('mixed', [False, True])
def test_species_cta_manifest_design_and_offline_replay(tmp_path, monkeypatch, modality, mixed, inputs):
    source = tmp_path / 'source'
    source.mkdir()
    reference = host()
    # This smoke exercises selection, not self-sharing; other-locus matches
    # are independently checked above.
    entries = [e for e in reference.proteome if e.gene_id == PRAME]
    context = replace(reference.proteome.reference,
        sequence_sha256=identity_digest(sorted((e.accession, e.sequence) for e in entries)))
    layers = [{**dict(reference.proteome.layers[0]), 'added': len(entries)}]
    reference = HostReference(Proteome(context, entries, layers))
    reference.save(source / 'host.json.gz')
    pd.DataFrame([(PRAME, 9)], columns=['feature_id', 'patient_tpm']).to_csv(
        source / 'expression.tsv', sep='\t', index=False)
    opts = options(overrides={PRAME: REVIEW.to_dict()})
    manifest = dict(schema='vaxrank.input_manifest.v1', host_species='dog',
        patient_id='dog-001', sample_id='dog-tumor', mhc_alleles=['DLA-88*001:01'],
        host_reference='host.json.gz', inputs=[dict(format='cta_expression', path='expression.tsv',
                                                  cta_expression=opts)])
    previous = None
    if mixed:
        from mhctools import RandomBindingPredictor
        from topiary import TopiaryPredictor
        scope = InputScope.from_dict(dict(patient_id='dog-001', sample_id='dog-tumor',
            host_species='dog', mhc_alleles=['DLA-88*001:01']))
        viral = replace(inputs, scope=scope,
            antigens={'viral': replace(inputs.antigens['viral'], amino_acids='YYYYYYYYYY')})
        viral.save(source / 'viral.json')
        mutation = replace(inputs, scope=scope, antigens={'mutation': replace(
            inputs.antigens['mutation'], species='dog', amino_acids='WWWWWWWWWW')})
        previous = mutation.predict(TopiaryPredictor(models=[RandomBindingPredictor(['DLA-88*001:01'])]),
            host_reference=reference, config=EpitopeConfig(score_expr='1.0'))
        previous.save(source / 'mutation.tsv')
        manifest['inputs'] += [dict(format='antigens', path='viral.json'),
                               dict(format='epitopes', path='mutation.tsv')]
    path = source / 'inputs.json'
    path.write_text(json.dumps(manifest))
    output = tmp_path / 'design'
    args = ['--input-manifest', str(path), '--mhc-predictor', 'random',
        '--output-dir', str(output), '--output-epitopes', str(output / 'native.tsv'),
        '--vaccine-type', modality, '--pdf-backend', 'weasyprint',
        '--config-text', 'epitopes.score_expr=1.0', '--no-processing-aware-annotation']
    if modality == 'mrna':
        args += ['--mrna-codon-method', 'use_best_codon', '--mrna-no-optimize-linkers']
    run_cli(args)
    dataset = EpitopeDataset.load(output / 'native.tsv')
    assert dataset.host_reference.identity == reference.identity
    assert {a.species for a in dataset.antigens.values() if a.kind == 'CTA'} == {'canis_lupus_familiaris'}
    saved = CTAExpressionResult.load(output / 'cta_admission.json')
    assert saved.admitted_antigens[0].gene_id == PRAME
    payload = (next(iter(dataset.selection['source_evidence'].values()))['cta_expression_admission']
               if mixed else dataset.selection['cta_expression_admission'])
    assert to_native_json(saved) == payload
    exported = test_antigen_input.AntigenDataset.from_epitope_dataset(dataset)
    target = next(a for a in exported.antigens.values() if a.kind == 'CTA')
    facts = dict(target.tumor_specificity.evidence_records[0].details)
    assert json.loads(facts['observations'])
    assert json.loads(facts['candidate_records'])[0]['taxon_id'] == 9615
    if mixed:
        assert {a.kind for a in dataset.antigens.values()} == {'CTA', 'mutation', 'viral'}
        historical = {e.prediction_id: e.predictions_flat() for e in previous.epitopes}
        assert {e.prediction_id: e.predictions_flat() for e in dataset.epitopes
                if e.prediction_id in historical} == historical
    shutil.rmtree(source)
    monkeypatch.setattr('oncoref.species_antigens.cta_snapshot', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(output / 'native.tsv'), '--output-dir', str(replay),
             '--pdf-backend', 'weasyprint', '--no-processing-aware-annotation'])
    fastas = list(output.glob('*.fasta'))
    assert fastas
    for file in [*fastas, output / 'manifest.json']:
        assert file.read_bytes() == (replay / file.name).read_bytes()
