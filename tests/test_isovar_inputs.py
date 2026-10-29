"""Producer contract replay, evidence gates and the common CLI/design path."""

from dataclasses import asdict, replace
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import gzip
import json

import pytest
from isovar import AlleleRead, ReadEvidence
from orfdata import TranscriptORF
from pyensembl import Genome
from pyensembl.species import find_species_by_name
from varcode import Variant

from vaxrank import (
    AntigenDataset, IsovarAdmission, IsovarSourceOptions, antigens_from_isovar,
    load_isovar_input, reconstruct_isovar_antigens,
)
from vaxrank.cli.entry_point import run_cli
from vaxrank.config.loader import load_vaxrank_config
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.input_scope import InputScope, read_input_manifest
from vaxrank.isovar_input import ISOVAR_INPUT_SCHEMA
from vaxrank.vaccine_antigen import TumorSpecificityAttestation
from . import test_antigen_input


DATA = Path(__file__).parent / 'data/isovar'
dependencies = test_antigen_input.dependencies
inputs = test_antigen_input.inputs
MANIFEST = json.loads((DATA / 'manifest.json').read_text())


def result(name):
    raw = gzip.decompress((DATA / (name + '.json.gz')).read_bytes())
    assert sha256(raw).hexdigest() == MANIFEST['files'][name]['sha256']
    return json.loads(raw)


def attestation():
    return TumorSpecificityAttestation(status='admitted', evidence_kind='synthetic_test_assertion',
        evidence_source=__name__, patient_specific=True, rationale_code='test_only')


@pytest.mark.parametrize('name', MANIFEST['files'])
def test_every_released_result_preserves_source_and_all_decisions(name, tmp_path):
    producer = result(name)
    dataset = antigens_from_isovar(producer)
    assert all(not a.tumor_specificity.admits_construct for a in dataset.antigens.values())
    evidence, = dataset.source_evidence.values()
    assert evidence['isovar_result'] == producer
    assert len(evidence['decisions']) == len(dataset.antigens)
    identity, = dataset.source_evidence
    assert all(dict(a.source_metadata)['isovar_evidence_id'] == identity for a in dataset.antigens.values())
    assert all('isovar_result' not in dict(a.source_metadata) for a in dataset.antigens.values())
    path = tmp_path / 'antigens.json'
    dataset.save(path)
    assert AntigenDataset.load(path) == dataset
    producer['source'] = 'changed after adaptation'
    assert evidence['isovar_result'] != producer


def test_ncorf_masks_synonymous_and_tiled_support():
    dataset = antigens_from_isovar(result('ncorf'), tumor_specificity=attestation())
    antigen, = dataset.antigens.values()
    assert antigen.tumor_specificity.admits_construct
    assert antigen.amino_acids == 'MAAAAAVAAAA'
    assert [(i.start, i.end) for i in antigen.targetable_mask.intervals] == [(6, 7)]
    assert antigen.transcript_ids == ('T',) and antigen.protein_ids == ('noncoding_orf',)
    synonymous, = antigens_from_isovar(result('synonymous'), tumor_specificity=attestation()).antigens.values()
    assert not synonymous.targetable_mask.intervals and not synonymous.tumor_specificity.admits_construct
    tiled = antigens_from_isovar(result('tiled'), tumor_specificity=attestation())
    for a in tiled.antigens.values():
        coding = json.loads(dict(a.source_metadata)['isovar_candidate'])['coding_support']
        if json.loads(dict(a.source_metadata)['isovar_candidate'])['sequence_extent'] == 'complete_coding_candidate':
            assert coding['single_read_span_support']['fragments'] == 0
        assert coding['fully_tiled']


def test_partial_ncorf_opt_in_keeps_all_real_colo829_alternatives():
    producer = result('colo829_dynll2')
    default = antigens_from_isovar(producer, tumor_specificity=attestation())
    assert len(default.antigens) == len(producer['protein_candidates']) == 300
    assert not any(a.tumor_specificity.admits_construct for a in default.antigens.values())
    relaxed = antigens_from_isovar(producer, tumor_specificity=attestation(),
                                  options=IsovarSourceOptions(admission=IsovarAdmission(
                                      allow_partial=True)))
    assert set(relaxed.antigens) == set(default.antigens)
    assert not any(a.tumor_specificity.admits_construct for a in relaxed.antigens.values())
    relaxed = antigens_from_isovar(producer, tumor_specificity=attestation(),
                                  options=IsovarSourceOptions(admission=IsovarAdmission(
                                      allow_partial=True, allow_ambiguous=True)))
    assert any(a.tumor_specificity.admits_construct for a in relaxed.antigens.values())
    assert all(a.gene_name == 'DYNLL2' and a.transcript_ids for a in relaxed.antigens.values())


def test_ncorf_frame_identity_is_independent_of_window_cropping():
    producer = result('ncorf')
    original, = producer['protein_candidates']
    cropped = deepcopy(original)
    # The producer frames each candidate separately. Removing noncoding
    # prefix changes local offsets, but cannot change the observed codons.
    cropped['candidate_id'] += ':rebased-window'
    cropped['nucleotide_sequence'] = cropped['nucleotide_sequence'][1:]
    for key in ('coding_interval', 'variant_cdna_interval'):
        cropped[key] = [value - 1 for value in cropped[key]]
    producer['protein_candidates'].append(cropped)
    producer['models'][0]['candidate_ids'].append(cropped['candidate_id'])
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation())
    assert len(dataset.antigens) == 2
    assert {a.amino_acids for a in dataset.antigens.values()} == {original['protein_sequence']}
    assert all(a.tumor_specificity.admits_construct for a in dataset.antigens.values())
    phases = {json.loads(dict(a.source_metadata)['isovar_candidate'])['variant_codon_phase']
              for a in dataset.antigens.values()}
    assert phases == {1}


def test_sv_annotated_and_exploratory_frames_keep_junction_masks():
    producer = result('sv')
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation(),
        options={'admission': {'allow_partial': True, 'allow_ambiguous': True, 'min_fragments': 1}})
    assert len(dataset.antigens) == 6
    antigen, = [a for key, a in dataset.antigens.items() if ':translation:' in key]
    translation, = producer['paths'][0]['translations']
    junction, = producer['paths'][0]['junctions']
    for start in range(len(antigen.amino_acids) - 8):
        expected = any(3 * start < b - translation['translation_start'] < 3 * (start + 9)
                       for b in junction['query_interval'])
        assert antigen.interval_is_targetable(start, start + 9) == expected
    assert antigen.kind == 'SV' and antigen.transcript_ids == ('D-201',)


def test_exported_sv_occurrences_do_not_borrow_union_support():
    producer = result('sv_occurrences')
    assert producer['candidates'][0]['rna_support']['fragments'] == 1
    assert all(o['witnesses'] for o in producer['candidates'][0]['occurrences'])
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation(),
        options={'admission': {'allow_ambiguous': True, 'min_fragments': 1}})
    assert len(dataset.antigens) == 2
    assert all('rna_fragment_support_unavailable' in a.tumor_specificity.rationale_code
               for a in dataset.antigens.values())


def measured_occurrences(second=True):
    from isovar.rna_evidence import evidence_set

    producer = result('sv_occurrences')
    candidate, = producer['candidates']
    identities = [('rg', 'fragment-A', 64), ('rg', 'fragment-B', 64)]
    scope = [producer['sample_id'], producer['source']]
    candidate['rna_support'].update(evidence_set(scope, identities))
    for occurrence, keys in zip(candidate['occurrences'], [identities, identities[:1] if second else []]):
        witness = occurrence['witnesses'][0]
        occurrence['witnesses'] = [dict(witness, observation=key[1]) for key in keys]
        occurrence['rna_support'] = dict(evidence_set(scope, keys),
                                         scope='full_interval_junction_linked_witnesses_for_occurrence')
    return producer


def test_measured_sv_placements_use_their_own_support_and_keep_admission_gates():
    producer = measured_occurrences()
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation(),
        options={'admission': {'allow_ambiguous': True, 'min_fragments': 2}})
    supports = []
    admitted = []
    for antigen in dataset.antigens.values():
        support = json.loads(dict(antigen.source_metadata)['isovar_candidate'])['occurrence_rna_support']
        supports.append(support['fragments'])
        admitted.append(antigen.tumor_specificity.admits_construct)
    assert supports == [2, 1] and admitted == [True, False]
    assert not any(a.tumor_specificity.admits_construct
                   for a in antigens_from_isovar(producer).antigens.values())
    for occurrence in producer['candidates'][0]['occurrences']:
        occurrence['frame_status'] = 'no_gene_anchor'
    assert not any(a.tumor_specificity.admits_construct for a in
                   antigens_from_isovar(producer, tumor_specificity=attestation()).antigens.values())
    assert next(iter(dataset.source_evidence.values()))['isovar_result']['candidates'][0]['rna_support']['fragments'] == 2


def test_measured_empty_sv_support_is_zero_and_not_unknown():
    from vaxrank.isovar_antigens import _sv_export

    candidates = list(_sv_export(measured_occurrences(second=False)))
    assert [c.fragments for c in candidates] == [2, 0]


def test_partially_measured_export_preserves_unmeasured_alternative():
    from vaxrank.isovar_antigens import _sv_export

    producer = measured_occurrences()
    del producer['candidates'][0]['occurrences'][1]['rna_support']
    assert [c.fragments for c in _sv_export(producer)] == [2, None]
    producer['candidates'][0]['occurrences'][1]['rna_support'] = None
    with pytest.raises(ValueError, match='scope'):
        list(_sv_export(producer))


@pytest.mark.parametrize('field,value', [
    ('scope', 'union_of_full_interval_junction_linked_witnesses_within_input'),
    ('evidence_scope', ['another-patient', 'another-library']),
    ('reads', True), ('fragments', 99), ('evidence_set_id', 'different'),
    ('read_ids', ['duplicate', 'duplicate']),
])
def test_inconsistent_measured_support_fails_without_union_fallback(field, value):
    producer = measured_occurrences()
    producer['candidates'][0]['occurrences'][0]['rna_support'][field] = value
    with pytest.raises(ValueError, match='Isovar SV RNA'):
        antigens_from_isovar(producer)


def test_fresh_raw_and_offline_sv_exports_have_the_same_placement_support():
    from isovar import export_sv_rna_orfs
    from vaxrank.isovar_antigens import _sv_export

    raw = result('sv')
    export = export_sv_rna_orfs(raw)
    offline = list(_sv_export(export))
    original = list(_sv_export(export, raw['observations']))
    assert offline == original and all(c.fragments is not None for c in offline)
    export = measured_occurrences()
    observations = {name: {'identity': ['rg', name, 64]} for name in ('fragment-A', 'fragment-B')}
    assert list(_sv_export(export)) == list(_sv_export(export, observations))
    observations['fragment-A']['identity'][1] = 'mismatching-ledger'
    with pytest.raises(ValueError, match='original observations'):
        list(_sv_export(export, observations))


def test_isoform_reference_only_never_becomes_observed_protein():
    dataset = antigens_from_isovar(result('isoforms'), tumor_specificity=attestation())
    assert len(dataset.antigens) == 2
    assert sum(a.tumor_specificity.admits_construct for a in dataset.antigens.values()) == 1
    partial = antigens_from_isovar(result('isoforms_partial'), tumor_specificity=attestation(),
                                 options={'admission': {'allow_partial': True}})
    for a in partial.antigens.values():
        evidence = json.loads(dict(a.source_metadata)['isovar_candidate'])
        if evidence['origin'] == 'reference_model':
            assert not a.tumor_specificity.admits_construct
            assert 'coding_interval_not_observed' in a.tumor_specificity.rationale_code
        assert a.kind == 'isoform' and a.gene_id == 'G'


def test_viral_coding_keeps_foreign_source_and_partial_boundaries():
    dataset = antigens_from_isovar(result('viral'), tumor_specificity=attestation())
    e7, = [a for a in dataset.antigens.values() if 'NP_041326.1' in a.protein_ids]
    assert len(e7.amino_acids) == 98
    assert e7.kind == 'viral' and e7.species and e7.species != 'homo_sapiens'
    assert e7.tumor_specificity.admits_construct
    evidence = json.loads(dict(e7.source_metadata)['isovar_candidate'])
    assert evidence['reference_orfs'][0]['protein_matches_reference'] is True
    partial, = antigens_from_isovar(result('viral_partial'), tumor_specificity=attestation()).antigens.values()
    assert not partial.tumor_specificity.admits_construct
    assert partial.amino_acids.startswith('L') and len(partial.amino_acids) == 20


def test_invalid_translation_reference_attestations_and_controls_fail():
    damaged = result('ncorf')
    damaged['protein_candidates'][0]['protein_sequence'] = 'WRONG'
    with pytest.raises(ValueError, match='protein disagrees'):
        antigens_from_isovar(damaged)
    with pytest.raises(ValueError, match='reference assemblies'):
        antigens_from_isovar(result('ncorf'), scope=InputScope(reference_assembly='other'))
    with pytest.raises(ValueError, match='unknown Isovar candidates'):
        antigens_from_isovar(result('ncorf'), candidate_attestations={'typo': attestation()})
    with pytest.raises(ValueError, match='saved Isovar result'):
        antigens_from_isovar(result('ncorf'), options={'reads': {'min_mapping_quality': 10}})
    with pytest.raises(ValueError, match='positive integer'):
        IsovarAdmission(min_fragments=0)
    with pytest.raises(ValueError, match='unknown field'):
        antigens_from_isovar(result('ncorf'), options={'admission': {'allow_partal': True}})
    for malformed in (False, [], {'admission': {'min_fragments': True}}):
        with pytest.raises(ValueError):
            antigens_from_isovar(result('ncorf'), options=malformed)


@pytest.fixture
def variant_inputs(tmp_path, request):
    # Tiny noncoding transcript adapted from Isovar's Apache-2.0 regression.
    parameter = getattr(request, 'param', 'GCC')
    competing = parameter == 'competing_frames'
    sequence = ('C' * 15 + 'ATGAATG' + 'GCC' * 10 + 'TAACCTAA' + 'G' * 15 if competing else
                'C' * 15 + 'ATG' + 'GCC' * 3 + parameter + 'GCC' * 6 + 'TAA' + 'G' * 15)
    attrs = 'gene_id "G"; gene_name "NONCODING"; transcript_id "T"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";'
    gtf, fasta = tmp_path / 'reference.gtf', tmp_path / 'reference.fa'
    gtf.write_text(''.join(f'1\ttest\t{kind}\t100\t{99 + len(sequence)}\t.\t+\t.\t{attrs}{extra}\n'
                           for kind, extra in [('gene', ''), ('transcript', ''),
                                                ('exon', ' exon_id "E"; exon_number "1";')]))
    fasta.write_text('>T\n' + sequence + '\n')
    genome = Genome('test_assembly', 'test_annotation', '1', gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(fasta)], cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name('homo_sapiens')
    genome.index()
    reference = TranscriptORF.at(genome.transcript_by_id('T'), 15,
                                 orf_id='noncoding_orf', source='synthetic')
    references = [reference]
    if competing:
        references.append(TranscriptORF.at(genome.transcript_by_id('T'), 19,
            orf_id='shifted_orf', source='synthetic'))
    variant = Variant('1', 134, sequence[34], 'T', genome)
    reads = [AlleleRead(sequence[:34], 'T', sequence[35:], name='read' + str(i),
                       reference_blocks=[(0, len(sequence), 99, 99 + len(sequence))]) for i in range(2)]
    yield dict(variant=variant, read_evidence=ReadEvidence.from_variant_and_allele_reads(variant, reads),
               reference_orfs=references, sample_id='sample', source='synthetic')
    genome.close()


@pytest.mark.parametrize('variant_inputs', ['competing_frames'], indirect=True)
def test_competing_frames_remain_ambiguous_after_window_rebasing(variant_inputs):
    from isovar import reconstruct_variant_orfs
    producer = reconstruct_variant_orfs(**variant_inputs)
    groups = {}
    for row in producer['protein_candidates']:
        groups.setdefault(row['rna_candidate_id'], []).append(row)
    competing = next(rows for rows in groups.values() if len({r['model_id'] for r in rows}) == 2)
    assert len(competing) == 2
    # Same producer proteins and mutation positions, with all noncoding prefix
    # removed. Their local coding starts now coincide despite distinct frames.
    for row in producer['protein_candidates']:
        prefix = row['coding_interval'][0]
        row['nucleotide_sequence'] = row['nucleotide_sequence'][prefix:]
        for key in ('coding_interval', 'variant_cdna_interval'):
            row[key] = [value - prefix for value in row[key]]
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation(),
                                   options={'admission': {'allow_partial': True}})
    assert all('ambiguous_coding_candidate' in dataset.antigens[row['candidate_id']].tumor_specificity.rationale_code
               for row in competing)
    relaxed = antigens_from_isovar(producer, tumor_specificity=attestation(),
        options={'admission': {'allow_partial': True, 'allow_ambiguous': True}})
    assert all(relaxed.antigens[row['candidate_id']].tumor_specificity.admits_construct for row in competing)


def test_live_python_api_uses_owner_controls_and_records_versions(variant_inputs):
    dataset = reconstruct_isovar_antigens('variant_orfs', inputs=variant_inputs,
        options={'protein': {'min_variant_sequence_coverage': 1, 'protein_sequence_preference': 'support'}},
        tumor_specificity=attestation())
    antigen, = dataset.antigens.values()
    assert antigen.tumor_specificity.admits_construct and antigen.amino_acids == 'MAAAAAVAAAA'
    evidence, = dataset.source_evidence.values()
    assert evidence['reconstruction']['controls']['protein']['min_variant_sequence_coverage'] == 1
    assert evidence['reconstruction']['controls']['protein']['max_protein_sequences_per_variant'] == 0
    assert evidence['reconstruction']['versions']['isovar']
    reference, = (m['reference'] for m in evidence['isovar_result']['models'])
    assert reference['resolved_transcript']['strand'] == '+'
    assert json.loads(dict(antigen.source_metadata)['isovar_candidate'])['reference_strand'] == '+'
    with pytest.raises(ValueError, match='unknown or inapplicable'):
        reconstruct_isovar_antigens('variant_orfs', inputs=variant_inputs, options={'protein': {'bogus': 1}})
    with pytest.raises(ValueError, match='ReadEvidence'):
        reconstruct_isovar_antigens('variant_orfs', inputs=variant_inputs,
                                    options={'reads': {'min_mapping_quality': 10}})


def test_manifest_yaml_cli_override_contract_and_relative_result(tmp_path):
    producer = tmp_path / 'producer.json'
    producer.write_text(json.dumps(result('colo829_dynll2')))
    path = tmp_path / 'source.json'
    path.write_text(json.dumps(dict(schema=ISOVAR_INPUT_SCHEMA, result=producer.name,
                                  tumor_specificity=asdict(attestation()))))
    from vaxrank.cli.arg_parser import choose_arg_parser
    argv = ['--input-isovar', str(path), '--set', 'isovar.sources.variant_orfs.admission.allow_partial=true']
    args = choose_arg_parser(argv).parse_args(argv)
    config = load_vaxrank_config(args)
    adapted = load_isovar_input(path, configured_sources=config['isovar']['sources'])
    direct = antigens_from_isovar(result('colo829_dynll2'), tumor_specificity=attestation(),
                                 options={'admission': {'allow_partial': True}})
    assert adapted == direct
    manifest = tmp_path / 'manifest.yaml'
    manifest.write_text('schema: vaxrank.input_manifest.v1\ninputs:\n'
                        '  - format: isovar\n    path: source.json\n'
                        '    isovar:\n      admission:\n        allow_partial: false\n')
    fmt, resolved_path, scope = read_input_manifest(manifest)[0]
    assert fmt == 'isovar' and resolved_path == str(path)
    overridden = load_isovar_input(path, options=scope['isovar'], configured_sources=config['isovar']['sources'])
    assert not any(a.tumor_specificity.admits_construct for a in overridden.antigens.values())


@pytest.mark.parametrize('name', ['ncorf', 'unresolved'])
def test_shared_cli_native_export_report_and_offline_replay(name, tmp_path, dependencies, monkeypatch):
    path = tmp_path / (name + '.json')
    producer = result(name)
    # Bind the synthetic producer to the shared deterministic host fixture.
    producer['variant']['reference_name'] = 'GRCh38'
    path.write_text(json.dumps(dict(schema=ISOVAR_INPUT_SCHEMA, result=producer,
        tumor_specificity=asdict(attestation()), scope=dict(patient_id='patient', sample_id='sample',
        reference_assembly='GRCh38', annotation='ensembl:93', mhc_alleles=['HLA-A*02:01', 'HLA-B*07:02']))))
    native, exported, ascii_path = tmp_path / 'predictions.tsv', tmp_path / 'antigens.json', tmp_path / 'report.txt'
    run_cli(['--input-isovar', str(path), '--mhc-predictor', 'random', '--ensembl-release', '93',
             '--output-antigens', str(exported), '--output-epitopes', str(native),
             '--output-ascii-report', str(ascii_path), '--report-candidate-details', '--no-processing-aware-annotation',
             '--vaccine-peptide-length', '10', '--config-text', 'epitopes.score_expr=1.0',
             '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score'])
    dataset = AntigenDataset.load(exported)
    assert len(dataset.antigens) == (1 if name == 'ncorf' else 0)
    assert dataset.source_evidence
    saved = EpitopeDataset.load(native)
    assert AntigenDataset.from_epitope_dataset(saved).source_evidence == dataset.source_evidence
    assert ascii_path.exists()
    if name == 'unresolved':
        assert 'no_alt_reads' in ascii_path.read_text()
    else:
        assert 'RNA reconstruction support — fragments: 2' in ascii_path.read_text()
        assert 'Reference ORF: noncoding_orf' in ascii_path.read_text()
        assert 'Construct admission: admitted' in ascii_path.read_text()
    replay = tmp_path / 'replay.json'
    html = tmp_path / 'replay.html'
    def forbidden(*args, **kwargs):
        raise AssertionError('Native report replay must not predict or reconstruct')
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.isovar_input.reconstruct_isovar_antigens', forbidden)
    run_cli(['--input-epitopes', str(native), '--output-antigens', str(replay),
             '--output-html-report', str(html), '--report-candidate-details', '--no-processing-aware-annotation'])
    assert AntigenDataset.load(replay).source_evidence == dataset.source_evidence
    if name == 'ncorf':
        assert 'RNA reconstruction support — fragments</th><td>2</td>' in html.read_text()
        assert 'Reference ORF</th><td>noncoding_orf</td>' in html.read_text()


def test_live_bam_request_and_saved_replay_share_owner_contract(variant_inputs, tmp_path):
    import pysam
    variant = variant_inputs['variant']
    observed = variant_inputs['read_evidence'].alt_reads[0].sequence
    bam_path = tmp_path / 'input.bam'
    header = {'SQ': [{'SN': '1', 'LN': 1000}], 'RG': [{'ID': 'library', 'SM': 'sample'}]}
    with pysam.AlignmentFile(bam_path, 'wb', header=header) as bam:
        for index in range(2):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = 'fragment' + str(index)
            read.query_sequence = observed
            read.reference_id, read.reference_start = 0, 99
            read.cigarstring = str(len(observed)) + 'M'
            read.mapping_quality = 60
            read.query_qualities = [30] * len(observed)
            read.set_tag('RG', 'library')
            bam.write(read)
    pysam.index(str(bam_path))
    path = tmp_path / 'request.json'
    request = dict(schema=ISOVAR_INPUT_SCHEMA, tumor_specificity=asdict(attestation()),
        options=dict(protein=dict(min_variant_sequence_coverage=1, protein_sequence_preference='support')),
        reconstruct=dict(family='variant_orfs', inputs=dict(
            bam=bam_path.name, variant=dict(contig='1', start=134, ref='C', alt='T'),
            reference_orfs=[dict(transcript_id='T', start=15, orf_id='noncoding_orf', source='synthetic')],
            sample_id='sample', source='synthetic')))
    path.write_text(json.dumps(request))
    live = load_isovar_input(path, genome=variant.genome)
    antigen, = live.antigens.values()
    assert antigen.amino_acids == 'MAAAAAVAAAA' and antigen.tumor_specificity.admits_construct
    evidence, = live.source_evidence.values()
    request.pop('reconstruct')
    request['result'] = evidence['isovar_result']
    request.pop('options')
    path.write_text(json.dumps(request))
    replay = load_isovar_input(path)
    assert replay.antigens == live.antigens
    assert next(iter(replay.source_evidence.values()))['isovar_result'] == evidence['isovar_result']


def test_isovar_manifest_combines_fresh_sequences_and_saved_predictions(inputs, dependencies, tmp_path):
    from topiary import TopiaryPredictor
    from vaxrank.epitope_config import EpitopeConfig
    model, genome = dependencies
    old = replace(inputs, antigens={'viral': inputs.antigens['viral']})
    saved = old.predict(TopiaryPredictor(models=[model]), genome=genome,
                        config=EpitopeConfig(score_expr='1.0'))
    saved.save(tmp_path / 'saved.tsv')
    producer = result('ncorf')
    producer['variant']['reference_name'] = 'GRCh38'
    (tmp_path / 'source.json').write_text(json.dumps(dict(schema=ISOVAR_INPUT_SCHEMA,
        result=producer, tumor_specificity=asdict(attestation()))))
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1', **asdict(inputs.scope),
        inputs=[dict(format='isovar', path='source.json'), dict(format='epitopes', path='saved.tsv')])))
    output = tmp_path / 'combined.tsv'
    run_cli(['--input-manifest', str(manifest), '--mhc-predictor', 'random', '--ensembl-release', '93',
             '--output-epitopes', str(output), '--no-processing-aware-annotation',
             '--config-text', 'epitopes.score_expr=1.0'])
    combined = EpitopeDataset.load(output)
    recovered = AntigenDataset.from_epitope_dataset(combined)
    assert len(recovered.antigens) == 2
    assert {a.kind for a in recovered.antigens.values()} == {'mutation', 'viral'}
    assert next(iter(recovered.source_evidence.values()))['isovar_result'] == producer
    original = {e.prediction_id: e.predictions_flat() for e in saved.epitopes}
    assert {e.prediction_id: e.predictions_flat() for e in combined.epitopes
            if e.prediction_id in original} == original


def test_fusion_live_api_and_existing_owner_json_input(tmp_path):
    from .test_isovar_fusion import example
    fusion, references, reads = example()
    inputs = dict(fusion=fusion, references=references, reads=reads)
    live = reconstruct_isovar_antigens('fusion', inputs=inputs, tumor_specificity=attestation())
    assert len(live.antigens) == 1
    assert next(iter(live.antigens.values())).tumor_specificity.admits_construct
    source = tmp_path / 'fusion-input.json'
    source.write_text(json.dumps(dict(fusion=asdict(fusion), references=[asdict(r) for r in references],
                                     reads=[asdict(r) for r in reads])))
    envelope = tmp_path / 'source.json'
    envelope.write_text(json.dumps(dict(schema=ISOVAR_INPUT_SCHEMA,
        reconstruct=dict(family='fusion', inputs=source.name), tumor_specificity=asdict(attestation()))))
    assert load_isovar_input(envelope).antigens == live.antigens
    stricter = reconstruct_isovar_antigens('fusion', inputs=inputs, options={'reconstruction': {'min_fragments': 3}})
    assert not stricter.antigens and stricter.source_evidence


@pytest.mark.parametrize('variant_inputs,table,amino', [('ATA', 2, 'M'), ('CTG', 3, 'T')],
                         indirect=['variant_inputs'])
def test_nonstandard_genetic_code_and_source_species_alias(variant_inputs, table, amino):
    from isovar import reconstruct_variant_orfs
    # These internal codons differ from table 1 (ATA=I, CTG=L).
    reference = variant_inputs['reference_orfs'][0]
    alternative = TranscriptORF.at(reference.transcript, 15, orf_id='alternate_code', source='synthetic',
                                   translation_table=table)
    producer = reconstruct_variant_orfs(**{**variant_inputs, 'reference_orfs': [alternative]})
    dataset = antigens_from_isovar(producer, species='Homo sapiens')
    assert dataset.antigens
    for a in dataset.antigens.values():
        assert a.amino_acids[4] == amino
        assert json.loads(dict(a.source_metadata)['isovar_candidate'])['translation_table'] == table
    producer['models'][0]['reference']['translation_table'] = 1
    with pytest.raises(ValueError, match='protein disagrees'):
        antigens_from_isovar(producer)


def test_ordinary_protein_export_keeps_all_translation_and_support_evidence(tmp_path):
    from types import SimpleNamespace
    from . import test_isovar_support_contract as contract
    from isovar import export_protein_hypotheses
    reference = contract.reference.__wrapped__(tmp_path, SimpleNamespace(param='+'))
    raw = contract.reconstruct(tmp_path, reference,
                               [(contract.MUTANT, [(len(contract.MUTANT), 'M')], 0)] * 2)
    producer = export_protein_hypotheses([raw], sample_id='sample', source='synthetic')
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation())
    antigen, = dataset.antigens.values()
    assert antigen.amino_acids == 'MAAAAQGGGG'
    assert antigen.interval_is_targetable(4, 8)
    assert next(iter(dataset.source_evidence.values()))['isovar_result'] == producer


def test_isovar_source_outcome_is_visible_even_without_proteins(tmp_path, dependencies):
    from vaxrank.isovar_antigens import isovar_evidence_summary
    dataset = antigens_from_isovar(result('unresolved'))
    text, = isovar_evidence_summary(dataset.source_evidence)
    assert '0 coding candidates' in text and 'no_alt_reads' in text


def _alignment(path, sequence, contig='1', start=0, cigar=None, count=2):
    import pysam
    header = {'SQ': [{'SN': contig, 'LN': 10000}], 'RG': [{'ID': 'library', 'SM': 'sample'}]}
    with pysam.AlignmentFile(path, 'wb', header=header) as bam:
        for index in range(count):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = 'fragment' + str(index)
            read.query_sequence = sequence
            read.reference_id, read.reference_start = 0, start
            read.cigarstring = cigar or str(len(sequence)) + 'M'
            read.mapping_quality = 60
            read.query_qualities = [30] * len(sequence)
            read.set_tag('RG', 'library')
            bam.write(read)
    pysam.index(str(path))


def test_live_isoform_request(tmp_path):
    producer = result('isoforms')
    target = producer['targets'][0]
    reference = target['reference']
    exons = reference['exons']
    cigar = ''
    previous = exons[0][0]
    for start, end in exons:
        if start != previous:
            cigar += str(start - previous) + 'N'
        cigar += str(end - start) + 'M'
        previous = end
    bam = tmp_path / 'isoforms.bam'
    _alignment(bam, reference['sequence'], start=exons[0][0], cigar=cigar)
    path = tmp_path / 'input.json'
    path.write_text(json.dumps(dict(schema=ISOVAR_INPUT_SCHEMA,
        tumor_specificity=asdict(attestation()), options=dict(reconstruction=dict(max_paths=8)),
        reconstruct=dict(family='isoforms', inputs=dict(bam=bam.name, targets=producer['targets'],
                                                        sample_id='sample', source='synthetic')))))
    dataset = load_isovar_input(path)
    assert sum(a.tumor_specificity.admits_construct for a in dataset.antigens.values()) == 1
    evidence, = dataset.source_evidence.values()
    assert evidence['reconstruction']['controls']['reconstruction']['max_paths'] == 8


def test_live_viral_request_uses_exact_synthetic_reference(tmp_path, monkeypatch):
    from orfdata import ReferenceContext, SourceRNA
    # A synthetic foreign taxon exercises reference-qualified coding without
    # downloading a viral genome or changing a biological viral sequence.
    dna = 'ATG' + 'GCT' * 20 + 'TAA'
    context = ReferenceContext(999999, 'synthetic_fixture', 'fixture', 'fixture', '1',
                               sequence_sha256=sha256(dna.encode()).hexdigest())
    rna = SourceRNA('toy-rna', dna, context, 'synthetic test', contig='foreign',
                    strand='+', exons=((1, len(dna)),))
    model = TranscriptORF.at(rna, 0, orf_id='toy-protein', source=rna.source)
    monkeypatch.setattr('isovar.viral_sequences.viral_coding_references',
                        lambda genome, contigs: ([model], dict(provider='synthetic_fixture', contigs=[])))
    bam, fasta = tmp_path / 'foreign.bam', tmp_path / 'reference.fa'
    _alignment(bam, dna, contig='foreign')
    # BAM dictionary identity is explicit even for the unused trailing bases.
    fasta.write_text('>foreign\n' + dna + 'C' * (10000 - len(dna)) + '\n')
    path = tmp_path / 'input.json'
    path.write_text(json.dumps(dict(schema=ISOVAR_INPUT_SCHEMA, tumor_specificity=asdict(attestation()),
        options=dict(reconstruction=dict(reference_profile=None, max_paths=8)),
        reconstruct=dict(family='viral', inputs=dict(bam=bam.name, genome_fasta=fasta.name,
                         contigs=['foreign'], sample_id='sample', source='synthetic',
                         reference_cache='cache', scratch_dir='.')))))
    dataset = load_isovar_input(path)
    antigen, = dataset.antigens.values()
    assert antigen.kind == 'viral' and antigen.amino_acids == 'M' + 'A' * 20
    assert antigen.species == 'synthetic_fixture' and antigen.tumor_specificity.admits_construct
    evidence, = dataset.source_evidence.values()
    assert evidence['reconstruction']['controls']['reconstruction']['reference_profile'] is None


def test_live_sv_request_retains_absence_of_rna(tmp_path):
    import pysam
    producer = result('sv')
    bam = tmp_path / 'empty.bam'
    with pysam.AlignmentFile(bam, 'wb', header={'SQ': [{'SN': c, 'LN': 10000} for c in ('1', '2')]}):
        pass
    pysam.index(str(bam))
    event = dict(event_id=producer['event_id'], reference_name=producer['reference_name'],
                 donor=producer['donor'], acceptor=producer['acceptor'], sample_id='sample',
                 references=[],
                 regions=[], event_provenance={'source': 'synthetic'})
    event_path = tmp_path / 'event.json'
    event_path.write_text(json.dumps(event))
    path = tmp_path / 'input.json'
    path.write_text(json.dumps(dict(schema=ISOVAR_INPUT_SCHEMA,
        options=dict(reconstruction=dict(max_paths=4)),
        reconstruct=dict(family='sv', inputs=dict(bam=bam.name, event=event_path.name, source='synthetic')))))
    dataset = load_isovar_input(path)
    assert not dataset.antigens and dataset.source_evidence
