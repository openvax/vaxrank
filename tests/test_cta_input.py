"""Expression-only design and native replay through the shared CLI."""

from dataclasses import replace
import json
from pathlib import Path
import shutil

import pandas as pd
import pytest

from . import test_cta_expression
from .test_cta_expression import (PRAME, MAGEA3, MAGEA4, NOT_CTA,
                                  ALTERNATE_TRANSCRIPT, PRAME_TRANSCRIPT)
from .test_external_rescoring import ContextPredictor, ALLELES
from .test_epitope_dataset import forbidden
from vaxrank.cli.entry_point import run_cli
from vaxrank.cta_expression import CTAExpressionResult
from vaxrank.cta_input import hitlist_annotations
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.native_serialization import from_native_json

genome = test_cta_expression.genome


@pytest.fixture
def cli_dependencies(monkeypatch, genome):
    import mhctools.cli
    import vaxrank.cli.entry_point
    import vaxrank.cta_input
    import vaxrank.safety_assessment
    model = ContextPredictor()
    genome.to_dict = lambda: dict(reference_name='GRCh38', annotation_version=93)
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', lambda args: [model])
    monkeypatch.setattr(vaxrank.cli.entry_point, 'resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', genome))
    monkeypatch.setattr(vaxrank.cta_input, 'self_reference_matches',
        lambda peptides, antigen, genome: {
            p: replace(antigen.self_reference_match(p, False), source_provenance_complete=True)
            for p in peptides})
    monkeypatch.setattr(vaxrank.safety_assessment, 'self_reference_matches',
                        vaxrank.cta_input.self_reference_matches)
    return model


def argv(path, output, level='gene'):
    return ['--input-cta-expression', str(path), '--cta-expression-level', level,
            '--cta-id-column', 'feature_id', '--cta-value-column', 'patient_tpm',
            '--cta-expression-unit', 'TPM', '--cta-sample-id', 'patient-001',
            '--cta-expression-source', 'Salmon', '--cta-expression-version', '1.10.3',
            '--cta-expression-assay', 'bulk RNA-seq', '--cta-min-expression', '2',
            '--ensembl-release', '93', '--mhc-predictor', 'random',
            '--mhc-alleles', ','.join(ALLELES), '--output-dir', str(output),
            '--no-processing-aware-annotation', '--num-epitopes-per-vaccine-peptide', '1000']


@pytest.mark.parametrize('level', ['gene', 'transcript'])
def test_expression_only_cli_outputs_and_offline_replay(tmp_path, monkeypatch, cli_dependencies, level):
    path = tmp_path / 'expression.tsv'
    rows = ([(PRAME, 3), (MAGEA3, 50), (MAGEA4, 4), (NOT_CTA, 10)] if level == 'gene'
            else [(ALTERNATE_TRANSCRIPT + '.2', 5), (PRAME_TRANSCRIPT + '.3', 0)])
    pd.DataFrame(rows, columns=['feature_id', 'patient_tpm']).to_csv(path, sep='\t', index=False)
    first, native = tmp_path / 'first', tmp_path / 'native.tsv'
    admission_path = tmp_path / 'admission.json'
    run_cli(argv(path, first, level) + ['--output-epitopes', str(native),
                                     '--output-cta-admission', str(admission_path),
                                     '--mrna-poly-a-length', '73',
                                     '--peptide-max-constructs', '1',
                                     '--peptide-n-terminal-acetyl',
                                     '--config-value', 'peptide.window_selection={"min_target_fraction":0.95}',
                                     '--vaccine-type', 'peptide', 'mrna'] +
            (['--mrna-junction-predictor', 'random'] if level == 'gene' else []))
    assert cli_dependencies.calls
    saved = EpitopeDataset.load(native)
    admission = from_native_json(saved.selection['cta_expression_admission'], CTAExpressionResult)
    assert admission == CTAExpressionResult.load(admission_path)
    assert saved.result.df.cta_expression_level.eq(level).all()
    assert saved.result.df.cta_expression_unit.eq('TPM').all()
    assert 'n_rna_alt' not in saved.result.df
    assert {a.gene_id for a in saved.antigens.values()} == ({PRAME, MAGEA4} if level == 'gene' else {PRAME})
    assert all(e.source_sequence[e.offset:e.offset + len(e.sequence)] == e.sequence for e in saved.epitopes)
    assert all(e.self_reference_match.source_provenance_complete for e in saved.epitopes)
    assert saved.selection['cta_product_audit']
    design = from_native_json((first / 'cta_design.json').read_text(), dict)
    assert design['final_product_audit']['mhc_status'] == 'predictions_returned'
    assert len(design['final_product_audit']['unassessed_processing_compartments']) == 3
    assert any('unsupported_chemical_modifications' in a.reason_codes
               for a in design['final_product_audit']['product_audits'])
    for region in design['final_product_audit']['selected_regions']:
        assert region['antigen'].amino_acids[region['start']:region['end']] == region['sequence']
    assert len({e.prediction_group_source for e in saved.epitopes}) == len(saved.epitopes)
    assert saved.selection['run_configuration']['effective_vaccine_peptides']['window_selection']
    peptide_fasta = first / 'peptide' / 'vaccine.fasta'
    assert peptide_fasta.exists() and (first / 'mrna' / 'cds.fasta').exists()
    assert all(len(s) <= 25 for s in peptide_fasta.read_text().splitlines() if not s.startswith('>'))
    path.unlink()
    monkeypatch.setattr(mhctools_cli(), 'predictors_from_args', forbidden)
    import vaxrank.cta_input
    monkeypatch.setattr(vaxrank.cta_input, 'admit_cta_expression', forbidden)
    monkeypatch.setattr(vaxrank.cta_input, 'load_hitlist_bundle', forbidden)
    run_cli(['--input-epitopes', str(native), '--output-dir', str(tmp_path / 'replay'),
             '--vaccine-type', 'peptide', 'mrna', '--no-processing-aware-annotation'])
    for file in first.rglob('*.fasta'):
        assert file.read_bytes() == (tmp_path / 'replay' / file.relative_to(first)).read_bytes()
    assert (first / 'cta_design.json').read_bytes() == (tmp_path / 'replay' / 'cta_design.json').read_bytes()
    assert (first / 'peptide' / 'window_selection.json').read_bytes() == (
        tmp_path / 'replay' / 'peptide' / 'window_selection.json').read_bytes()
    if level == 'transcript':
        # A capacity override must invalidate its audit even when this one
        # source still produces exactly the same product sequences.
        limits = tmp_path / 'changed_limits'
        run_cli(['--input-epitopes', str(native), '--output-dir', str(limits),
                 '--peptide-max-constructs', '20', '--no-processing-aware-annotation'])
        assert (first / 'peptide' / 'vaccine.fasta').read_bytes() == (limits / 'peptide' / 'vaccine.fasta').read_bytes()
        audit = from_native_json((limits / 'cta_design.json').read_text(), dict)['final_product_audit']
        assert audit['mhc_status'] == 'unassessed'
    # Explicit new settings cannot inherit either old products or their audit.
    changed = tmp_path / 'changed'
    run_cli(['--input-epitopes', str(native), '--output-dir', str(changed),
             '--vaccine-type', 'peptide', 'mrna', '--peptide-max-constructs', '20',
             '--mrna-no-optimize-linkers'] +
            (['--mrna-poly-a-length', '120'] if level == 'gene' else ['--mrna-poly-a-length=120']) + [
             '--config-value', 'peptide.n_terminal_acetyl=false',
             '--no-processing-aware-annotation'])
    changed_design = from_native_json((changed / 'cta_design.json').read_text(), dict)
    assert changed_design['final_product_audit']['mhc_status'] == 'unassessed'
    assert changed_design['final_product_audit']['reason'] == 'mhc_predictor_not_requested'
    assert all(p.poly_a_nt == 'A' * 120 for p in changed_design['final_product_audit']['products']['mrna'])
    if level == 'gene':
        assert len(changed_design['final_product_audit']['products']['peptide']) > 1
    assert (first / 'peptide' / 'manifest.json').read_bytes() != (changed / 'peptide' / 'manifest.json').read_bytes()


def test_expression_dsl_changes_selected_targets(tmp_path, cli_dependencies):
    path = tmp_path / 'expression.tsv'
    pd.DataFrame([(PRAME, 3), (MAGEA4, 4)], columns=['feature_id', 'patient_tpm']).to_csv(path, sep='\t', index=False)
    out = tmp_path / 'selected'
    run_cli(argv(path, out) + ['--vaccine-type', 'peptide',
        '--config-text', 'epitopes.filter_expr=cta_gene_tpm > 3',
        '--config-text', 'epitopes.score_expr=cta_gene_tpm'])
    decisions = pd.read_csv(out / 'cta_target_decisions.csv').set_index('gene_id')
    assert decisions.loc[PRAME, 'admission_status'] == 'admitted'
    assert decisions.loc[PRAME, 'design_status'] == 'no_eligible_region'
    assert decisions.loc[MAGEA4, 'design_status'] == 'selected_region'
    assert 'Expression features' in (out / 'vaccine_report.txt').read_text()


@pytest.mark.parametrize('level', ['gene', 'transcript'])
def test_verified_alternate_sources_survive_cli_native_occurrences(tmp_path, monkeypatch, level):
    from .test_cta_identity import Annotation, ALTERNATE
    import oncoref
    import vaxrank.cli.entry_point
    annotation = Annotation(112)
    annotation.to_dict = lambda: dict(reference_name='GRCh38', annotation_version=112)
    # Exercise real admission and self matching on the pinned annotation, with
    # a deterministic predictor so both identical source occurrences survive.
    import mhctools.cli
    import vaxrank.cta_input
    predictor = ContextPredictor()
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', lambda args: [predictor])
    monkeypatch.setattr(vaxrank.cli.entry_point, 'resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', annotation))
    rows = ([(PRAME, 3), (ALTERNATE, 9)] if level == 'gene' else
            [('ENST00000398743', 3), ('ENST00000617728', 9)])
    path, native = tmp_path / 'expression.tsv', tmp_path / 'native.tsv'
    pd.DataFrame(rows, columns=['feature_id', 'patient_tpm']).to_csv(path, sep='\t', index=False)
    first = tmp_path / 'first'
    run_cli(argv(path, first, level) + ['--ensembl-release', '112', '--output-epitopes', str(native),
                                      '--vaccine-type', 'peptide'])
    saved = EpitopeDataset.load(native)
    assert {a.gene_id for a in saved.antigens.values()} == {PRAME, ALTERNATE}
    assert saved.result.df.canonical_gene_id.eq(PRAME).all()
    assert set(saved.result.df.cta_expression_value) == {3, 9}
    assert {e.source_name for e in saved.epitopes} == {'patient-001:' + r[0] for r in rows}
    assert len({e.prediction_group_source for e in saved.epitopes}) == len(saved.epitopes)
    alternate = next(a for a in saved.antigens.values() if a.gene_id == ALTERNATE)
    expected_transcript = 'ENST00000539862' if level == 'gene' else 'ENST00000617728'
    assert alternate.transcript_ids == (expected_transcript,)
    assert alternate.protein_ids == (annotation.transcript_by_id(expected_transcript).protein_id,)
    path.unlink()
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    monkeypatch.setattr(vaxrank.cta_input, 'admit_cta_expression', forbidden)
    monkeypatch.setattr(oncoref, 'resolve_gene_identity', forbidden)
    monkeypatch.setattr(oncoref, 'cta_annotation_gene_identities', forbidden)
    run_cli(['--input-epitopes', str(native), '--output-dir', str(tmp_path / 'replay'),
             '--vaccine-type', 'peptide', '--no-processing-aware-annotation'])
    assert (first / 'cta_design.json').read_bytes() == (tmp_path / 'replay' / 'cta_design.json').read_bytes()


def mhctools_cli():
    import mhctools.cli
    return mhctools.cli


def test_missing_contract_fails_before_prediction(tmp_path, monkeypatch):
    monkeypatch.setattr(mhctools_cli(), 'predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='cta-expression-level'):
        run_cli(['--input-cta-expression', 'missing.tsv', '--output-dir', str(tmp_path / 'out')])


def test_ms_annotations_never_imply_unobserved_nested_peptide_support():
    bundle = {'artifacts': {'peptides.parquet': {'rows': [dict(
        peptide='ABCDEFGHIJK', n_ms_observations=2, cta_specific=True,
        avoid_sequence=False, tissue_review_required=False, n_essential_tissue_donors=0)]}}}
    frame = hitlist_annotations(['ABCDEFGHIJK', 'ABCDEFGHI', 'OTHER'], bundle)
    assert frame.hitlist_query_status.tolist() == ['captured', 'unqueried', 'unqueried']
    assert frame.hitlist_ms_observations[0] == 2
    assert pd.isna(frame.hitlist_ms_observations[1])
    assert pd.isna(frame.hitlist_avoid_sequence[2])
    assert json.loads(frame.to_json(orient='records'))[1]['hitlist_ms_observations'] is None


def test_ms_restriction_categories_and_rejected_observations_are_distinct():
    sequence = 'ACDEFGHIK'
    observations = [dict(peptide=sequence, mhc_restriction=ALLELES[0], restriction_evidence=category)
                    for category in ('monoallelic', 'experimental', 'predicted', 'unknown')]
    observations.append(dict(peptide=sequence, mhc_restriction='HLA class I',
        restriction_evidence='unknown', mhc_allele_provenance='sample_allele_match',
        mhc_allele_set=json.dumps(ALLELES)))
    bundle = {'artifacts': {
        'peptides.parquet': {'rows': [dict(peptide=sequence, n_ms_observations=5)]},
        'presentation.parquet': {'rows': observations},
        'excluded_observations.parquet': {'rows': [dict(peptide=sequence, assay_method='fluorescence')]},
    }}
    frame = hitlist_annotations([sequence, sequence], bundle, ALLELES)
    assert frame.hitlist_ms_observations.tolist() == [5, 5]
    assert frame.hitlist_ms_monoallelic.tolist() == [1, 0]
    assert frame.hitlist_ms_experimental_restriction.tolist() == [1, 0]
    assert frame.hitlist_ms_predicted_restriction.tolist() == [1, 0]
    assert frame.hitlist_ms_donor_hla_candidate.tolist() == [1, 1]
    assert bundle['artifacts']['excluded_observations.parquet']['rows'][0]['assay_method'] == 'fluorescence'
    alias = hitlist_annotations([sequence], bundle, ['A*02:01'])
    assert alias.hitlist_ms_monoallelic.tolist() == [1]


def test_hitlist_bundle_contract_capture_and_mismatch(tmp_path, monkeypatch, genome):
    import hashlib
    from importlib.metadata import version
    evidence_bundle = pytest.importorskip('hitlist.evidence_bundle')
    from vaxrank.cta_input import load_hitlist_bundle
    from .test_cta_expression import run
    admission = run(tmp_path, genome, [(PRAME, 3)])
    bundle = tmp_path / 'bundle'
    bundle.mkdir()
    frame = pd.DataFrame([dict(peptide='ACDEFGHIK', n_ms_observations=1)])
    frame.to_parquet(bundle / 'peptides.parquet', index=False)
    pd.DataFrame([dict(peptide='ACDEFGHIK', assay_method='fluorescence')]).to_parquet(
        bundle / 'excluded_observations.parquet', index=False)
    (bundle / 'lineage.json').write_text('{"fixture": "all source records retained"}')
    manifest = dict(kind='cta_expression', expression=dict(
        input=dict(sha256=admission.input_sha256), level='gene', id_column='feature_id',
        tpm_column='patient_tpm', units='TPM', ensembl_release=93, oncoref_version=version('oncoref')),
        artifacts={p.name: dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in bundle.iterdir()})
    calls = []
    def verify(path):
        calls.append(path)
        return manifest
    monkeypatch.setattr(evidence_bundle, 'verify_evidence_bundle', verify)
    captured = load_hitlist_bundle(bundle, admission)
    assert len(calls) == 2
    assert captured['artifacts']['peptides.parquet']['rows'] == frame.to_dict('records')
    assert captured['artifacts']['excluded_observations.parquet']['rows'][0]['assay_method'] == 'fluorescence'
    assert json.loads(captured['artifacts']['lineage.json'])['fixture']
    manifest['expression']['ensembl_release'] = 112
    with pytest.raises(ValueError, match='contract differs'):
        load_hitlist_bundle(bundle, admission)
    manifest['expression']['ensembl_release'] = 93
    monkeypatch.setattr(evidence_bundle, 'verify_evidence_bundle',
                        lambda path: (_ for _ in ()).throw(ValueError('artifact mismatch')))
    with pytest.raises(ValueError, match='artifact mismatch'):
        load_hitlist_bundle(bundle, admission)


def test_real_verified_hitlist_bundle_cli_and_portable_replay(tmp_path, monkeypatch, cli_dependencies):
    """Exercise the public verifier on a real producer-exported synthetic bundle."""
    import base64
    import hashlib
    evidence_bundle = pytest.importorskip('hitlist.evidence_bundle')
    fixture = Path(__file__).parent / 'data' / 'cta_hitlist_bundle'
    bundle = tmp_path / 'moved-evidence'
    shutil.copytree(fixture, bundle)
    expression = tmp_path / 'expression.tsv'
    expression.write_bytes((fixture / 'input-expression.tsv').read_bytes())
    manifest = evidence_bundle.verify_evidence_bundle(bundle)
    # This frozen producer bundle records its original reference version. The
    # annotation/admission dependencies above are synthetic, so reproduce that
    # declared version for this portable-bundle test. Do not relabel old evidence
    # as having been generated by the current OncoRef release.
    import importlib.metadata
    actual_version = importlib.metadata.version
    monkeypatch.setattr(importlib.metadata, 'version', lambda name:
                        manifest['expression']['oncoref_version'] if name == 'oncoref' else actual_version(name))
    first, native = tmp_path / 'first', tmp_path / 'native.tsv'
    run_cli(argv(expression, first) + ['--hitlist-evidence-bundle', str(bundle),
        '--output-epitopes', str(native), '--vaccine-type', 'peptide', 'mrna'])
    saved = EpitopeDataset.load(native)
    captured = saved.selection['hitlist_evidence_bundle']
    assert captured['manifest'] == manifest
    for name, facts in manifest['artifacts'].items():
        raw = base64.b64decode(captured['files_base64'][name])
        assert raw == (bundle / name).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == facts['sha256']
    assert len(captured['artifacts']['mappings.parquet']['rows']) == 2
    assert len(captured['artifacts']['contributors.parquet']['rows']) == 3
    rejected = captured['artifacts']['excluded_observations.parquet']['rows']
    assert len(rejected) == 1 and rejected[0]['assay_method'] == 'fluorescence'
    exact = hitlist_annotations(['SLYNTVATL', 'SLYNTVAT'], captured)
    assert exact.hitlist_ms_observations[0] == 1
    assert exact.hitlist_avoid_sequence[0]  # Shared non-CTA and benign-tissue sources.
    assert pd.isna(exact.hitlist_ms_observations[1])
    assert saved.result.df.hitlist_ms_observations.isna().all()
    expression.unlink()
    shutil.rmtree(bundle)
    monkeypatch.setattr(mhctools_cli(), 'predictors_from_args', forbidden)
    monkeypatch.setattr(evidence_bundle, 'verify_evidence_bundle', forbidden)
    run_cli(['--input-epitopes', str(native), '--output-dir', str(tmp_path / 'replay'),
        '--vaccine-type', 'peptide', 'mrna', '--no-processing-aware-annotation'])
    assert (first / 'cta_design.json').read_bytes() == (tmp_path / 'replay' / 'cta_design.json').read_bytes()
    for file in first.rglob('*.fasta'):
        assert file.read_bytes() == (tmp_path / 'replay' / file.relative_to(first)).read_bytes()
