"""Synthetic complete panels exercise the public species-neutral admission API.

All sequences, source hashes and observations here are invented fixtures. They
test the contract and never stand in for a validated human, dog or mouse panel.
"""

from dataclasses import asdict, replace
import hashlib
import json
import shutil

from oncoref.cta_evidence_adapters import native_observation_facts
from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_antigens import CTAEvidenceSnapshot
from orfdata import Occurrence, ProteinSequence, Proteome
from orfdata.reference import ReferenceContext
import pandas as pd
import pytest

from vaxrank.cta_admission import CTAAdmissionPolicy
from vaxrank.cta_expression import CTAExpressionResult, admit_cta_expression
from vaxrank.cta_identity import identity_digest
from vaxrank.cta_restriction import CTARestrictionAssessment, host_occurrences
from vaxrank.host_reference import HostReference
from vaxrank.native_serialization import from_native_json, to_native_json

from .test_cta_expression import contract, SEQUENCE
from .test_cta_restriction import evidence, policy


def fixture(tmp_path, taxon=9615, namespace='ncbi', shared=False, counterevidence=False, native_rna=None):
    species = {9606: 'homo_sapiens', 9615: 'canis_lupus_familiaris', 10090: 'mus_musculus'}[taxon]
    if namespace == 'ncbi':
        namespaces = dict(gene='NCBI Gene', transcript='NCBI RefSeq RNA', protein='NCBI RefSeq protein')
        roles = ('protein', 'feature_table')
        genes, transcripts, proteins = ['1', '2'], ['NM_1.1', 'NM_2.1'], ['NP_1.1', 'NP_2.1']
        annotation = 'NCBI RefSeq'
    else:
        namespaces = dict(gene='Ensembl gene', transcript='Ensembl transcript', protein='Ensembl protein')
        roles = ('gtf', 'rna', 'ncrna', 'protein')
        genes = ['ENSG00000000001', 'ENSG00000000002']
        transcripts = ['ENST00000000001.1', 'ENST00000000002.1']
        proteins = ['ENSP00000000001.1', 'ENSP00000000002.1']
        annotation = 'Ensembl'
    assets = {role: hashlib.sha256(('synthetic:' + role).encode()).hexdigest() for role in roles}
    reference = ReferenceIdentity(taxon, 'GCA_000000001.1', annotation, '1',
                                  'synthetic-fixture', namespaces, assets)
    source = dict(id='synthetic', version='1', url='https://example.invalid/fixture',
        doi=None, pmid=None, sha256='c' * 64, scope='Synthetic test fixture',
        study_id='synthetic', lineage=['synthetic'])
    expression = ExpressionContract(taxon, 'synthetic', '1', 'normal', None, 'unresolved',
        namespaces['gene'], ExpressionMeasurement('gene', 'TPM', 'source_native', False),
        {'source': source['sha256']}, 'Synthetic nomination; independent RNA is in the supplemental bundle')
    if native_rna is not None:
        expression = replace(expression, reference=reference, reference_status='verified',
            feature_namespace=namespaces['transcript'],
            measurement=ExpressionMeasurement('transcript', 'TPM', 'source_native', False))
    candidates, observations, facts, entries = [], [], [], []
    for index in range(2 if shared else 1):
        gid, tid, pid = genes[index], transcripts[index], proteins[index]
        eid, cid = f'synthetic:{index}', f'candidate:{index}'
        candidates.append(dict(id=cid, label=f'FIXTURE{index}', nomination_kind='assay_followup',
            derivation=[dict(kind='assay_followup', source_id='synthetic', locator=eid, note='Synthetic fixture')],
            evidence_status='limited', evidence_ids=[eid], notes='Synthetic fixture', mappings=[dict(
                reference_id=reference.reference_id, gene_namespace=namespaces['gene'], status='unique',
                method='Synthetic identity', complete=True, notes='Synthetic fixture',
                targets=[dict(id=gid, symbol=f'FIXTURE{index}',
                    sequences=[dict(transcript_id=tid, protein_id=pid)])])]))
        contrary = counterevidence and index == 1
        row = dict(id=eid, candidate_ids=[cid], source_id='synthetic', study_id='synthetic',
            contract_id=expression.contract_id, locator=eid, assay='Synthetic nomination',
            tissue='heart' if contrary else None, tissue_role='somatic' if contrary else 'unknown',
            sample_ids=[], sample_count=None, donor_count=None, technical_replicates=None,
            observation=dict(status='qualitative' if contrary else 'not_reported', value=None,
                reported='Synthetic heart signal' if contrary else 'Not reported', positive_count=None),
            source_fields={}, notes='Synthetic fixture')
        if native_rna is not None:
            row.update(tissue='heart', tissue_role='somatic', observation=dict(
                status=native_rna[0], value=native_rna[1], reported='Synthetic native RNA', positive_count=None))
        observations.append(row)
        facts.extend(native_observation_facts(row, expression.as_dict(), source)[0])
        origin = Occurrence(tid, gid, f'FIXTURE{index}', 'protein_coding', '1', '+',
                            0, len(SEQUENCE) * 3, ((1, len(SEQUENCE) * 3),), True)
        entries.append(ProteinSequence(pid, SEQUENCE, 'annotated', taxon,
                                       'synthetic fixture', occurrences=(origin,)))
    snapshot = CTAEvidenceSnapshot.create(dict(schema='oncoref.cta-evidence', schema_version=2,
        taxon_id=taxon, source_version='synthetic-fixture',
        coverage=dict(status='partial', limitations=['Synthetic fixture only']),
        references=[dict(id='synthetic', identity=reference.as_dict())], sources=[source],
        contracts=[expression.as_dict()], candidates=candidates, evidence=observations, facts=facts, curation=[]))
    context = ReferenceContext(taxon, species, reference.assembly_accession, annotation, '1',
        assembly_accession=reference.assembly_accession, source_sha256=identity_digest(assets),
        sequence_sha256=identity_digest(sorted((e.accession, e.sequence) for e in entries)))
    host = HostReference(Proteome(context, entries, (dict(role='host', kind='annotated',
        source='synthetic fixture', content_provenance=assets, identifier_namespaces=namespaces, added=len(entries)),)))
    supplemental = replace(evidence(taxon, targets=host_occurrences(host, SEQUENCE)),
        reference_id=reference.reference_id, host_reference_sha256=host.identity, snapshot_id=snapshot.snapshot_id)
    path = tmp_path / 'snapshot.json'
    snapshot.write(path)
    return host, path, supplemental, genes, transcripts


def admit(tmp_path, host, path, supplemental, rows, level='gene'):
    table = tmp_path / 'patient.tsv'
    pd.DataFrame(rows, columns=['feature_id', 'patient_tpm']).to_csv(table, sep='\t', index=False)
    return admit_cta_expression(table, input_contract=contract(level), admission_policy=CTAAdmissionPolicy(2, 'TPM'),
        host_reference=host, evidence_snapshot=path, restriction_policy=policy(), restriction_evidence=supplemental)


@pytest.mark.parametrize('taxon', [9606, 9615, 10090])
@pytest.mark.parametrize('namespace', ['ensembl', 'ncbi'])
def test_common_api_admits_equivalent_complete_evidence_and_replays(tmp_path, monkeypatch, taxon, namespace):
    host, path, supplemental, genes, transcripts = fixture(tmp_path, taxon, namespace)
    result = admit(tmp_path, host, path, supplemental, [(genes[0], 10)])
    antigen, = result.admitted_antigens
    assert antigen.tumor_specificity.status == 'admitted'
    assert antigen.tumor_specificity.requires_review
    assert antigen.self_reference_excluded_gene_ids == (genes[0],)
    assert result.decisions[0].assessment.tumor_expression.measurement_level == 'gene'
    assert dict(antigen.source_metadata)['sequence_selection'] == 'annotation_longest_protein'
    assert result.input_contract.id_namespace == ('ncbi_gene_id' if namespace == 'ncbi' else 'ensembl_gene_id')
    saved = tmp_path / 'result.json'
    result.save(saved)
    path.unlink()
    def forbidden(*args, **kwargs):
        raise AssertionError('Replay queried live evidence')
    monkeypatch.setattr('oncoref.species_antigens.cta_snapshot', forbidden)
    assert CTAExpressionResult.load(saved) == result
    restriction = result.decisions[0].assessment.restriction
    assert from_native_json(to_native_json(restriction), CTARestrictionAssessment) == restriction
    host.validate_antigen(antigen)


def test_shared_native_heart_evidence_cannot_be_hidden_by_an_independent_negative_panel(tmp_path):
    host, path, supplemental, genes, _ = fixture(tmp_path, shared=True, counterevidence=True)
    result = admit(tmp_path, host, path, supplemental, [(genes[0], 10)])
    assert not result.admitted_antigens
    assessment = result.decisions[0].assessment.restriction
    assert any(f.code == 'somatic_fact_requires_review' and f.evidence_id == 'synthetic:1' for f in assessment.findings)
    native = json.loads(assessment.native_evidence_json)
    assert {c['id'] for c in native['candidates']} == {'candidate:0', 'candidate:1'}
    assert native['facts'][1]['tissue'] == 'heart'
    with pytest.raises(ValueError, match='findings disagree'):
        replace(assessment, findings=())


@pytest.mark.parametrize('taxon', [9606, 9615, 10090])
@pytest.mark.parametrize('status,value,expected', [('reported', 0.5, None), ('below_limit', 0.5, None),
    ('reported', 10, 'somatic_fact_above_bound'), ('below_limit', 10, 'somatic_fact_requires_review')])
def test_verified_native_rna_is_compared_without_inventing_coverage(tmp_path, taxon, status, value, expected):
    host, path, supplemental, genes, _ = fixture(tmp_path, taxon, native_rna=(status, value))
    result = admit(tmp_path, host, path, supplemental, [(genes[0], 10)])
    assessment = result.decisions[0].assessment.restriction
    assert {f.code for f in assessment.findings} == ({expected} if expected else set())
    # Compatible low native facts never replace independently allocated coverage.
    empty = replace(supplemental, observations=())
    held = admit(tmp_path, host, path, empty, [(genes[0], 10)])
    assert not held.admitted_antigens
    assert 'insufficient_independent_donors' in {f.code for f in held.decisions[0].assessment.restriction.findings}


@pytest.mark.parametrize('same_sequence', [False, True])
def test_unresolved_extra_orf_transcript_is_retained_and_only_holds_equivalent_target(tmp_path, same_sequence):
    host, path, supplemental, genes, _ = fixture(tmp_path)
    origin = replace(host.proteome.proteins[0].occurrences[0], gene_id='2', transcript_id='')
    extra = ProteinSequence('unresolved-extra', SEQUENCE if same_sequence else 'M' * 12,
        'extra_orf', 9615, 'synthetic-extra', occurrences=(origin,))
    host = HostReference(Proteome(host.proteome.reference, [*host.proteome, extra],
        (*host.proteome.layers, dict(role='host', kind='extra_orf', source='synthetic-extra'))))
    supplemental = replace(supplemental, host_reference_sha256=host.identity)
    result = admit(tmp_path, host, path, supplemental, [(genes[0], 10)])
    assessment = result.decisions[0].assessment.restriction
    assert bool(result.admitted_antigens) is not same_sequence
    assert ('unresolved_host_transcript' in {f.code for f in assessment.findings}) is same_sequence


def test_wrong_reference_and_sequence_rejected_on_ingest_and_replay(tmp_path):
    host, path, supplemental, genes, _ = fixture(tmp_path)
    for change in (dict(host_reference_sha256='f' * 64), dict(taxon_id=10090),
                   dict(reference_id='openvax:reference:v1:sha256:' + 'f' * 64)):
        with pytest.raises(ValueError, match='host/reference/snapshot'):
            admit(tmp_path, host, path, replace(supplemental, **change), [(genes[0], 10)])
    row = supplemental.observations[0]
    corrupt = replace(row, occurrences=(replace(row.occurrences[0], sequence_sha256='f' * 64),))
    with pytest.raises(ValueError, match='occurrence or sequence'):
        admit(tmp_path, host, path, replace(supplemental,
            observations=(corrupt,) + supplemental.observations[1:]), [(genes[0], 10)])
    result = admit(tmp_path, host, path, supplemental, [(genes[0], 10)])
    restriction = result.decisions[0].assessment.restriction
    with pytest.raises(ValueError, match='findings disagree'):
        replace(restriction, evidence=None)


def test_exact_transcript_versions_and_patient_expression_remain_independent(tmp_path):
    host, path, supplemental, genes, transcripts = fixture(tmp_path, namespace='ensembl')
    result = admit(tmp_path, host, path, supplemental, [(transcripts[0], 0)], level='transcript')
    assert not result.admitted_antigens
    assert result.decisions[0].assessment.restriction.passed
    assert result.decisions[0].assessment.antigen.tumor_specificity.rationale_code == 'cta_expression_below_threshold'
    with pytest.raises(ValueError, match='exact versioned Ensembl'):
        admit(tmp_path, host, path, supplemental, [(transcripts[0].split('.')[0], 10)], level='transcript')


def test_older_supported_oncoref_reports_common_policy_capability_before_source_access(tmp_path, monkeypatch):
    host, path, supplemental, genes, _ = fixture(tmp_path)
    monkeypatch.setattr('oncoref.__version__', '1.19.2')
    def forbidden(*args, **kwargs):
        raise AssertionError('Unsupported common policy accessed source data')
    monkeypatch.setattr('oncoref.species_antigens.cta_snapshot', forbidden)
    with pytest.raises(ValueError, match='requires OncoRef >=1.23.1'):
        admit(tmp_path, host, path, supplemental, [(genes[0], 10)])


def test_common_admission_excludes_only_its_target_and_checks_all_saved_occurrences(tmp_path):
    from vaxrank.species_cta import _restriction_digest, _restriction_record
    host, path, supplemental, genes, _ = fixture(tmp_path, shared=True)
    result = admit(tmp_path, host, path, supplemental, [(genes[0], 10)])
    antigen, = result.admitted_antigens
    match = host.matches([SEQUENCE[:9]], antigen)[SEQUENCE[:9]]
    assert {s.gene_id for s in match.sources} == {genes[1]}
    with pytest.raises(ValueError, match='complete host reference'):
        host.validate_antigen(replace(antigen, amino_acids='Y' * len(antigen.amino_acids)))
    original = result.decisions[0].assessment.restriction
    shortened = replace(original, occurrences=original.occurrences[:1])
    records = tuple(r for r in antigen.tumor_specificity.evidence_records if r.evidence_kind != 'cta_rna_restriction')
    altered = replace(antigen, tumor_specificity=replace(antigen.tumor_specificity,
        evidence_records=records + (_restriction_record(shortened, genes[0]),)),
        source_metadata=tuple((k, _restriction_digest(shortened) if k == 'cta_restriction_sha256' else v)
                              for k, v in antigen.source_metadata))
    with pytest.raises(ValueError, match='complete host reference'):
        host.validate_antigen(altered)
    with pytest.raises(ValueError, match='lost its frozen restriction'):
        host.validate_antigen(replace(antigen, tumor_specificity=replace(antigen.tumor_specificity,
            evidence_records=records)))


def test_common_dog_construct_checks_final_windows_and_keeps_unsupported_dla_unknown(tmp_path):
    from mhctools.unsupported_allele import UnsupportedAllele
    from vaxrank.construct_audit import ConstructAudit, audit_construct_sequence
    from vaxrank.construct_sequence import ConstructEvidence, ConstructSequence, ConstructSequenceEdit
    from .test_construct_audit import RecordingPredictor

    host, path, supplemental, genes, _ = fixture(tmp_path)
    antigen, = admit(tmp_path, host, path, supplemental, [(genes[0], 10)]).admitted_antigens
    origin = ConstructEvidence(source='synthetic software fixture', evidence_level='reported',
                               description='Invented terminal residues for boundary validation')
    n = len(antigen.amino_acids)
    record = ConstructSequence(name='synthetic dog', sequence=antigen.amino_acids + 'KK',
        modality='peptide', evidence=origin, native_antigen=antigen, native_start=0, native_end=n,
        sequence_edits=(ConstructSequenceEdit(n, n, 'KK', origin),))

    class DogPredictor(RecordingPredictor):
        def predict_from_named_sequences(self, sequences):
            if self.failure:
                raise UnsupportedAllele('Fixture model does not support DLA-88*001:01')
            frame = super().predict_from_named_sequences(sequences)
            frame['allele'] = 'DLA-88*001:01'
            return frame

    predictor = DogPredictor()
    audit = audit_construct_sequence(record, predictor, host_reference=host)
    assert predictor.calls == [{record.cache_identity: record.sequence}]
    assert len(audit.mhc_assessment.ligands) == len(record.sequence) - 8
    assert audit.mhc_assessment.ligands[-1].crosses_construct_boundary
    assert {p.allele for ligand in audit.mhc_assessment.ligands for p in ligand.predictions} == {'DLA-88*001:01'}
    failed = audit_construct_sequence(record, DogPredictor(failure=True), host_reference=host)
    assert failed.mhc_status == 'failed' and failed.mhc_assessment is None
    assert 'DLA-88*001:01' in failed.error_message
    assert from_native_json(to_native_json(failed), ConstructAudit) == failed


@pytest.mark.parametrize('taxon,host_species,allele', [
    (9606, 'human', 'HLA-A*02:01'), (9615, 'dog', 'DLA-88*001:01'), (10090, 'mouse', 'H-2-Kb')])
@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_common_manifest_native_constructs_and_offline_replay(tmp_path, monkeypatch, taxon, host_species, allele, modality):
    from vaxrank.antigen_input import AntigenDataset
    from vaxrank.cli.entry_point import run_cli
    from vaxrank.epitope_dataset import EpitopeDataset
    from .test_cta_manifest import options
    source = tmp_path / 'source'
    source.mkdir()
    host, snapshot_path, supplemental, genes, _ = fixture(source, taxon,
        'ensembl' if taxon == 9606 else 'ncbi')
    host.save(source / 'host.json.gz')
    (source / 'restriction.json').write_text(json.dumps(asdict(supplemental)))
    pd.DataFrame([(genes[0], 10)], columns=['feature_id', 'patient_tpm']).to_csv(
        source / 'expression.tsv', sep='\t', index=False)
    manifest = dict(schema='vaxrank.input_manifest.v1', host_species=host_species,
        patient_id='synthetic-patient', sample_id='synthetic-tumor', mhc_alleles=[allele],
        host_reference='host.json.gz', inputs=[dict(format='cta_expression', path='expression.tsv',
            cta_expression=options(evidence_snapshot=snapshot_path.name,
                restriction_policy=asdict(policy()), restriction_evidence='restriction.json'))])
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
    assert dataset.host_reference.identity == host.identity
    admission = CTAExpressionResult.load(output / 'cta_admission.json')
    assert admission.decisions[0].assessment.restriction.passed
    exported = AntigenDataset.from_epitope_dataset(dataset)
    antigen, = exported.antigens.values()
    host.validate_antigen(antigen)
    assert antigen.targetable_mask == admission.admitted_antigens[0].targetable_mask
    assert any(r.evidence_kind == 'cta_rna_restriction' for r in antigen.tumor_specificity.evidence_records)
    shutil.rmtree(source)
    def forbidden(*args, **kwargs):
        raise AssertionError('Replay acquired live evidence or models')
    monkeypatch.setattr('oncoref.species_antigens.cta_snapshot', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(output / 'native.tsv'), '--output-dir', str(replay),
             '--pdf-backend', 'weasyprint', '--no-processing-aware-annotation'])
    fastas = list(output.glob('*.fasta'))
    assert fastas
    for file in [*fastas, output / 'manifest.json']:
        assert file.read_bytes() == (replay / file.name).read_bytes()
