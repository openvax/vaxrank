"""Expression-first CTA inputs for the shared prediction/design pipeline."""

from dataclasses import replace
import base64
import hashlib
import io
import json
from pathlib import Path

import pandas as pd
from topiary import TopiaryPredictor, TopiaryResult

from .cta_admission import CTAAdmissionPolicy
from .cta_expression import CTAExpressionInput, CTATargetExclusionPolicy, admit_cta_expression
from .epitope_config import EpitopeConfig
from .epitope_dataset import ANTIGEN_COLUMN, EpitopeDataset
from .external_report import ExternalReport
from .input_scope import InputProvenance, InputScope, normalize_alleles
from .native_serialization import to_native_json
from .reference_proteome import self_reference_matches


def validate_cta_args(args):
    """Reject ambiguous input contracts before reference loading/prediction."""
    if not getattr(args, 'input_cta_expression', None):
        if getattr(args, 'hitlist_evidence_bundle', None):
            raise ValueError('--hitlist-evidence-bundle requires --input-cta-expression')
        return
    conflicts = ('vcf', 'bam', 'input_lens', 'input_pvacseq', 'input_epitopes',
                 'input_topiary', 'external_input', 'input_manifest')
    if any(getattr(args, name, None) for name in conflicts):
        raise ValueError('--input-cta-expression requires a single expression source without VCF/BAM')
    required = ('cta_expression_level', 'cta_id_column', 'cta_value_column',
                'cta_expression_unit', 'cta_sample_id', 'cta_expression_source',
                'cta_expression_version', 'cta_expression_assay', 'cta_min_expression',
                'ensembl_release', 'mhc_predictor')
    missing = ['--' + name.replace('_', '-') for name in required
               if getattr(args, name, None) is None]
    if missing:
        raise ValueError('CTA expression input requires ' + ', '.join(missing))
    if getattr(args, 'external_predictions', 'input') != 'input' or getattr(args, 'external_peptide_only', False):
        raise ValueError('CTA expression inputs predict native proteins; external prediction modes do not apply')
    if getattr(args, 'prediction_cache', None):
        raise ValueError('CTA expression context prediction caches are not yet supported')
    if getattr(args, 'cta_no_gene_exclusions', False) and args.cta_exclude_gene_pattern:
        raise ValueError('--cta-no-gene-exclusions conflicts with --cta-exclude-gene-pattern')


def load_hitlist_bundle(path, admission):
    """Capture the verified public contract, retaining rejected evidence too."""
    try:
        from importlib.metadata import version
        from packaging.version import Version
        from hitlist.evidence_bundle import verify_evidence_bundle
    except ImportError as error:
        raise ValueError('Hitlist evidence requires pip install "vaxrank[hitlist]" (Hitlist >=1.65.1)') from error
    if Version(version('hitlist')) < Version('1.65.1'):
        raise ValueError('Hitlist CTA evidence requires Hitlist >=1.65.1')
    directory = Path(path)
    manifest = verify_evidence_bundle(directory)
    if manifest['kind'] != 'cta_expression':
        raise ValueError('Expected a Hitlist CTA-expression evidence bundle')
    expression = manifest['expression']
    contract = admission.input_contract
    expected = dict(level=contract.measurement_level, id_column=contract.id_column,
                    tpm_column=contract.value_column, units=contract.expression_unit,
                    ensembl_release=int(admission.annotation_version))
    if expression.get('oncoref_version') != version('oncoref'):
        raise ValueError('Hitlist bundle and CTA admission use different OncoRef reference versions')
    if (expression['input']['sha256'] != admission.input_sha256
            or any(expression.get(k) != v for k, v in expected.items())):
        raise ValueError('Hitlist bundle expression/reference contract differs from this input')
    artifacts = {}
    files = {}
    for name in manifest['artifacts']:
        file = directory / name
        raw = file.read_bytes()
        if hashlib.sha256(raw).hexdigest() != manifest['artifacts'][name]['sha256']:
            raise ValueError('Hitlist evidence artifact changed during capture: ' + name)
        files[name] = base64.b64encode(raw).decode('ascii')
        if file.suffix == '.parquet':
            # Pandas JSON normalizes nullable scalars and list-valued columns.
            frame = pd.read_parquet(io.BytesIO(raw))
            artifacts[name] = dict(columns=list(frame.columns),
                                   rows=json.loads(frame.to_json(orient='records')))
        else:
            artifacts[name] = raw.decode('utf8')
    # Catch concurrent edits between verification and capture.
    if manifest != verify_evidence_bundle(directory):
        raise ValueError('Hitlist evidence changed during capture')
    return dict(manifest=manifest, artifacts=artifacts, files_base64=files)


def hitlist_annotations(sequences, bundle, alleles=None):
    """Exact public observations; an unqueried sequence has unknown support."""
    columns = ('hitlist_ms_observations', 'hitlist_cta_specific', 'hitlist_avoid_sequence',
               'hitlist_tissue_review_required', 'hitlist_essential_tissue_donors')
    summaries = ({row['peptide']: row for row in bundle['artifacts']['peptides.parquet']['rows']}
                 if bundle else {})
    source_fields = ('n_ms_observations', 'cta_specific', 'avoid_sequence',
                     'tissue_review_required', 'n_essential_tissue_donors')
    from mhcgnomes import ParseError
    from .allele_validation import parse_prediction_allele
    observations = {}
    if bundle:
        for row in bundle['artifacts'].get('presentation.parquet', {}).get('rows', []):
            observations.setdefault(row['peptide'], []).append(row)
    rows = []
    for sequence, allele in zip(sequences, alleles if alleles is not None else [''] * len(sequences)):
        allele = parse_prediction_allele(allele) if allele else ''
        summary = summaries.get(sequence)
        row = dict(zip(columns, (summary.get(k) if summary is not None else None
                                 for k in source_fields)))
        row['hitlist_query_status'] = 'captured' if summary is not None else 'unqueried'
        counts = dict(hitlist_ms_monoallelic=0, hitlist_ms_experimental_restriction=0,
                      hitlist_ms_predicted_restriction=0, hitlist_ms_donor_hla_candidate=0,
                      hitlist_ms_coarse_or_untyped=0, hitlist_ms_unknown_restriction_evidence=0)
        for observation in observations.get(sequence, ()):
            try:
                restriction = parse_prediction_allele(observation.get('mhc_restriction', ''))
            except (ParseError, TypeError, ValueError):
                restriction = None  # Coarse/untyped restriction remains in captured raw evidence.
            evidence = observation.get('restriction_evidence')
            counts['hitlist_ms_coarse_or_untyped'] += int(restriction is None)
            counts['hitlist_ms_unknown_restriction_evidence'] += int(evidence not in {'monoallelic', 'experimental', 'predicted'})
            if restriction == allele:
                key = {'monoallelic': 'hitlist_ms_monoallelic',
                       'experimental': 'hitlist_ms_experimental_restriction',
                       'predicted': 'hitlist_ms_predicted_restriction'}.get(evidence)
                if key:
                    counts[key] += 1
            if observation.get('mhc_allele_provenance') in {'sample_allele_match', 'sample_locus_match'}:
                candidates = observation.get('mhc_allele_set') or []
                if isinstance(candidates, str):
                    candidates = json.loads(candidates)
                counts['hitlist_ms_donor_hla_candidate'] += int(allele in candidates)
        row.update({k: v if summary is not None else None for k, v in counts.items()})
        rows.append(row)
    return pd.DataFrame(rows)


def prepare_cta_report(args, genome, epitope_config=None):
    """Predict all admitted protein occurrences before applying shared DSL."""
    from mhctools.cli import mhc_alleles_from_args, predictors_from_args
    validate_cta_args(args)
    contract = CTAExpressionInput(
        args.cta_expression_level, args.cta_id_column, args.cta_value_column,
        args.cta_expression_unit, args.cta_sample_id, args.cta_expression_source,
        args.cta_expression_version, args.cta_expression_assay)
    exclusions = CTATargetExclusionPolicy(
        () if args.cta_no_gene_exclusions else
        tuple(args.cta_exclude_gene_pattern) if args.cta_exclude_gene_pattern is not None else ('MAGE*',),
        tuple(args.cta_allow_gene) if args.cta_allow_gene is not None else ('MAGEA4',))
    admission = admit_cta_expression(
        args.input_cta_expression, input_contract=contract,
        admission_policy=CTAAdmissionPolicy(args.cta_min_expression, contract.expression_unit),
        exclusion_policy=exclusions, genome=genome)
    bundle = (load_hitlist_bundle(args.hitlist_evidence_bundle, admission)
              if args.hitlist_evidence_bundle else None)
    if args.output_cta_admission:
        Path(args.output_cta_admission).parent.mkdir(parents=True, exist_ok=True)
        admission.save(args.output_cta_admission)
    if not admission.admitted_antigens:
        raise ValueError('No admitted CTA proteins; inspect --output-cta-admission decisions')
    alleles = normalize_alleles(mhc_alleles_from_args(args), 'CTA patient genotype')
    provenance = InputProvenance(
        input_id='cta:' + hashlib.sha256((admission.input_sha256 + to_native_json(contract)).encode()).hexdigest(),
        source_format='cta_expression', path=str(args.input_cta_expression),
        content_sha256=admission.input_sha256,
        scope=InputScope(patient_id=args.output_patient_id or contract.sample_id,
                         reference_assembly=admission.reference_assembly,
                         annotation=admission.annotation_name + ':' + admission.annotation_version,
                         mhc_alleles=alleles, sample_id=contract.sample_id))
    models = predictors_from_args(args)
    if not models:
        raise ValueError('CTA inputs require an MHC predictor')
    predictor = TopiaryPredictor(models=models)
    antigens = {a.source_identifier: a for a in admission.admitted_antigens}
    frame = predictor.predict_from_named_sequences({name: a.amino_acids for name, a in antigens.items()})
    if frame.empty:
        raise ValueError('CTA predictor returned no peptide measurements')
    frame = frame.reset_index(drop=True)
    from .allele_validation import validate_prediction_frame, parse_prediction_allele
    validate_prediction_frame(frame, 'CTA native protein prediction')
    if any(a and parse_prediction_allele(a) not in alleles for a in frame.allele):
        raise ValueError('CTA predictor returned an allele outside the declared patient genotype')
    decisions = {d.assessment.antigen.source_identifier: d for d in admission.decisions if d.status == 'admitted'}
    for column, values in {
        'source_sequence': [antigens[n].amino_acids for n in frame.source_sequence_name],
        ANTIGEN_COLUMN: [to_native_json(antigens[n]) for n in frame.source_sequence_name],
        'gene_id': [antigens[n].gene_id for n in frame.source_sequence_name],
        'gene_name': [antigens[n].display_gene_name for n in frame.source_sequence_name],
        'transcript_id': [decisions[n].transcript_id for n in frame.source_sequence_name],
        'cta_expression_value': [decisions[n].value for n in frame.source_sequence_name],
    }.items():
        frame[column] = values
    frame['cta_expression_unit'] = contract.expression_unit
    frame['cta_expression_level'] = contract.measurement_level
    if contract.expression_unit == 'TPM':
        frame['cta_' + contract.measurement_level + '_tpm'] = frame.cta_expression_value
    frame['sample_name'] = contract.sample_id
    frame['source_class'] = 'self'
    frame['source_species'] = 'Homo sapiens'
    annotations = hitlist_annotations(frame.peptide, bundle, frame.allele)
    for name in annotations:
        frame[name] = annotations[name]
    adapted = EpitopeDataset.from_topiary(TopiaryResult(frame), label=provenance.input_id)
    def occurrence_id(name, peptide, offset):
        return provenance.input_id + ':' + hashlib.sha256(
            json.dumps([name, peptide, int(offset)]).encode()).hexdigest()
    # As in direct variant prediction, each native protein occurrence is a
    # design input. Shared peptide strings never choose an arbitrary gene or
    # erase the distinct flank context/expression of another occurrence.
    frame['prediction_id'] = [occurrence_id(r.source_sequence_name, r.peptide, r.peptide_offset)
                              for r in frame.itertuples()]
    candidates = tuple(replace(e, prediction_id=occurrence_id(e.source_name, e.sequence, e.offset))
                       for e in adapted.epitopes)
    dataset = EpitopeDataset.from_predictions(candidates, result=TopiaryResult(frame))
    dataset.antigens = {e.prediction_id: antigens[e.source_name] for e in candidates}
    # Exact source-level self matches use the unchanged full candidate CTA
    # reference, independently of which targets were admitted/excluded.
    matches = {name: self_reference_matches(frame.loc[frame.source_sequence_name.eq(name), 'peptide'], antigen, genome)
               for name, antigen in antigens.items()}
    dataset.epitopes = tuple(replace(
        e, input_provenance=provenance, patient_alleles=alleles,
        occurs_in_reference=True,
        occurs_in_non_CTA_reference=matches[e.source_name][e.sequence].occurs,
        self_reference_match=matches[e.source_name][e.sequence]) for e in dataset.epitopes)
    dataset.result.df['occurs_in_reference'] = True
    dataset.result.df['occurs_in_non_CTA_reference'] = [
        matches[r.source_sequence_name][r.peptide].occurs for r in dataset.result.df.itertuples()]
    dataset.provenance = (provenance,)
    dataset.config = epitope_config or EpitopeConfig()
    dataset.selection['cta_expression_admission'] = to_native_json(admission)
    if bundle is not None:
        dataset.selection['hitlist_evidence_bundle'] = bundle
    report_frame = dataset.report_frame()
    report_frame.attrs['cta_predictor'] = predictor
    return ExternalReport('cta_expression', str(args.input_cta_expression),
                          report_df=report_frame, epitopes=dataset.epitopes,
                          dataset=dataset, input_provenance=provenance)
