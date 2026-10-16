"""Expression-first CTA inputs for the shared prediction/design pipeline."""

from dataclasses import asdict, dataclass, fields, replace
import base64
import hashlib
import io
import json
from pathlib import Path

import pandas as pd
from topiary import TopiaryPredictor

from .cta_admission import CTAAdmissionPolicy
from .cta_expression import CTAExpressionInput, CTAExpressionResult, CTATargetExclusionPolicy, admit_cta_expression
from .antigen_input import predict_antigens
from .external_report import ExternalReport
from .input_scope import InputProvenance, InputScope, combine_scopes, normalize_alleles, scope_dict, validate_input_scopes
from .native_serialization import to_native_json
from .reference_proteome import self_reference_matches


def validate_cta_args(args):
    """Reject ambiguous input contracts before reference loading/prediction."""
    standalone = bool(getattr(args, 'input_cta_expression', None))
    if not standalone:
        if getattr(args, 'hitlist_evidence_bundle', None):
            raise ValueError('--hitlist-evidence-bundle requires --input-cta-expression')
        explicit = getattr(args, '_explicit_cli_args', ())
        if any(name.startswith('cta_') for name in explicit):
            raise ValueError('Declare CTA expression options per input in --input-manifest')
    required = ('cta_expression_level', 'cta_id_column', 'cta_value_column',
                'cta_expression_unit', 'cta_sample_id', 'cta_expression_source',
                'cta_expression_version', 'cta_expression_assay', 'cta_min_expression',
                'ensembl_release', 'mhc_predictor')
    missing = ['--' + name.replace('_', '-') for name in required
               if getattr(args, name, None) is None]
    if standalone and missing:
        raise ValueError('CTA expression input requires ' + ', '.join(missing))
    from .external_rescoring import external_input_specs
    if not any(fmt == 'cta_expression' for fmt, _, _ in external_input_specs(args)):
        return
    if any(getattr(args, name, None) for name in ('vcf', 'bam')):
        raise ValueError('CTA expression inputs can join saved VCF/BAM evidence using epitopes inputs')
    if not getattr(args, 'mhc_predictor', None) or getattr(args, 'ensembl_release', None) is None:
        raise ValueError('CTA expression inputs require --mhc-predictor and --ensembl-release')
    if getattr(args, 'external_predictions', 'input') != 'input' or getattr(args, 'external_peptide_only', False):
        raise ValueError('CTA expression inputs predict native proteins; external prediction modes do not apply')
    if getattr(args, 'prediction_cache', None):
        raise ValueError('CTA expression context prediction caches are not yet supported')
    if getattr(args, 'cta_no_gene_exclusions', False) and args.cta_exclude_gene_pattern:
        raise ValueError('--cta-no-gene-exclusions conflicts with --cta-exclude-gene-pattern')


def cta_options_from_args(args):
    """Translate the existing single-source flags into the manifest contract."""
    names = dict(measurement_level='expression_level', id_column='id_column',
                 value_column='value_column', expression_unit='expression_unit',
                 evidence_source='expression_source', evidence_version='expression_version',
                 assay='expression_assay', min_expression='min_expression')
    options = {key: getattr(args, 'cta_' + flag, None) for key, flag in names.items()}
    exclusions = {}
    if getattr(args, 'cta_no_gene_exclusions', False):
        exclusions['gene_patterns'] = []
    elif getattr(args, 'cta_exclude_gene_pattern', None) is not None:
        exclusions['gene_patterns'] = args.cta_exclude_gene_pattern
    if getattr(args, 'cta_allow_gene', None) is not None:
        exclusions['gene_exceptions'] = args.cta_allow_gene
    if exclusions:
        options['exclusion_policy'] = exclusions
    if getattr(args, 'hitlist_evidence_bundle', None):
        options['hitlist_evidence_bundle'] = str(args.hitlist_evidence_bundle)
    return options


def cta_expression_contract(options, scope, location):
    """Validate per-source options without annotation, admission or prediction."""
    names = {f.name for f in fields(CTAExpressionInput)} - {'sample_id'}
    allowed = names | {'min_expression', 'exclusion_policy', 'hitlist_evidence_bundle'}
    if not isinstance(options, dict) or set(options) - allowed:
        raise ValueError(f'{location}: invalid or unknown cta_expression options')
    missing = sorted(names - set(options))
    if missing:
        raise ValueError(f'{location}: CTA expression requires {", ".join(missing)}')
    minimum = options.get('min_expression')
    if isinstance(minimum, bool) or not isinstance(minimum, (int, float)):
        raise ValueError(f'{location}: CTA min_expression must be an explicit number')
    contract = CTAExpressionInput(sample_id=scope.sample_id, **{k: options[k] for k in names})
    policy = CTAAdmissionPolicy(minimum, contract.expression_unit)
    excluded = options.get('exclusion_policy', {})
    if not isinstance(excluded, dict) or set(excluded) - {'gene_patterns', 'gene_exceptions'}:
        raise ValueError(f'{location}: invalid or unknown CTA exclusion_policy')
    if any(not isinstance(v, (list, tuple)) for v in excluded.values()):
        raise ValueError(f'{location}: CTA exclusions must be lists')
    exclusions = CTATargetExclusionPolicy(**excluded)
    bundle = options.get('hitlist_evidence_bundle')
    if 'hitlist_evidence_bundle' in options and (not isinstance(bundle, str) or not bundle.strip()):
        raise ValueError(f'{location}: hitlist_evidence_bundle must be a nonempty path')
    return contract, policy, exclusions, bundle


@dataclass(frozen=True)
class CTAExpressionSource:
    """Experimental prepared admission evidence, before model acquisition."""

    admission: CTAExpressionResult
    provenance: InputProvenance
    bundle: dict | None = None

    @classmethod
    def load(cls, path, *, options, scope, genome, manifest_path=None, declarations=None):
        contract, policy, exclusions, bundle_path = cta_expression_contract(options, scope, str(path))
        if genome is None:
            raise ValueError('CTA expression requires a host annotation')
        if not scope.patient_id or not scope.mhc_alleles:
            raise ValueError('CTA expression requires declared patient_id and mhc_alleles')
        scope = combine_scopes(scope, InputScope(
            reference_assembly=str(genome.reference_name),
            annotation=str(genome.annotation_name) + ':' + str(genome.annotation_version)), str(path))
        validate_input_scopes([InputProvenance('validation', 'cta_expression', str(path), '', scope)], genome=genome)
        admission = admit_cta_expression(path, input_contract=contract, admission_policy=policy,
                                        exclusion_policy=exclusions, genome=genome)
        bundle = load_hitlist_bundle(bundle_path, admission) if bundle_path else None
        identity = dict(admission=to_native_json(admission), scope=scope_dict(scope),
                        bundle_manifest=bundle['manifest'] if bundle is not None else None)
        provenance = InputProvenance(
            input_id='cta:' + hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
            source_format='cta_expression', path=str(path), content_sha256=admission.input_sha256,
            scope=scope, manifest_path=str(Path(manifest_path).resolve()) if manifest_path else None,
            manifest_declarations=dict(declarations or {}))
        return cls(admission, provenance, bundle)


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


def write_cta_admissions(args, sources):
    """Write source-specific admission records without merging their decisions."""
    if not sources:
        return
    explicit = getattr(args, 'output_cta_admission', None)
    directory = getattr(args, 'output_dir', None)
    if explicit and len(sources) != 1:
        raise ValueError('--output-cta-admission requires one CTA expression input; '
                         'use --output-dir for separate per-source files')
    if not explicit and not directory:
        return
    for source in sources:
        if explicit:
            path = Path(explicit)
        elif len(sources) == 1:
            path = Path(directory) / 'cta_admission.json'
        else:
            digest = source.provenance.input_id.split(':', 1)[1]
            path = Path(directory) / 'cta_admissions' / (digest + '.json')
        path.parent.mkdir(parents=True, exist_ok=True)
        source.admission.save(path)


def prepare_cta_report(args, genome, epitope_config=None):
    """Compatibility entry point using the shared input preparation API."""
    from mhctools.cli import mhc_alleles_from_args, predictors_from_args
    from .external_rescoring import external_input_specs, read_reports
    validate_cta_args(args)
    specs = external_input_specs(args)
    if len(specs) != 1 or specs[0][0] != 'cta_expression':
        raise ValueError('prepare_cta_report requires one CTA expression input')

    def factory(scope):
        args.mhc_alleles = ','.join(scope.mhc_alleles)
        args.mhc_alleles_file = None
        return TopiaryPredictor(models=predictors_from_args(args))

    return read_reports(
        [(fmt, path) for fmt, path, _ in specs], scopes=[s for _, _, s in specs],
        genome=genome, epitope_config=epitope_config, antigen_predictor_factory=factory,
        antigen_alleles=normalize_alleles(mhc_alleles_from_args(args), 'CTA patient genotype'),
        cta_admission_writer=lambda sources: write_cta_admissions(args, sources))[0]


def predict_cta_source(source, predictor, *, genome, epitope_config=None, host_reference=None):
    """Predict admitted occurrences, retaining every source admission decision."""
    admission, provenance, bundle = source.admission, source.provenance, source.bundle
    contract = admission.input_contract
    antigens = {a.source_identifier: a for a in admission.admitted_antigens}
    decisions = {d.assessment.antigen.source_identifier: d for d in admission.decisions if d.status == 'admitted'}

    def annotate(frame):
        frame['canonical_gene_id'] = [decisions[n].reference_resolution.evidence.gene_id
                                      for n in frame.source_sequence_name]
        frame['transcript_id'] = [decisions[n].transcript_id for n in frame.source_sequence_name]
        frame['cta_expression_value'] = [decisions[n].value for n in frame.source_sequence_name]
        frame['cta_expression_unit'] = contract.expression_unit
        frame['cta_expression_level'] = contract.measurement_level
        if contract.expression_unit == 'TPM':
            frame['cta_' + contract.measurement_level + '_tpm'] = frame.cta_expression_value
        frame['source_class'] = 'self'
        annotations = hitlist_annotations(frame.peptide, bundle, frame.allele)
        for name in annotations:
            frame[name] = annotations[name]

    if antigens:
        dataset = predict_antigens(antigens, predictor, genome=genome, provenance=provenance,
                                   config=epitope_config, annotate=annotate,
                                   self_matcher=self_reference_matches, host_reference=host_reference)
    else:
        from .epitope_dataset import EpitopeDataset
        dataset = EpitopeDataset.from_predictions((), config=epitope_config)
        dataset.host_reference = host_reference
        dataset.provenance = (provenance,)
        dataset.selection['antigen_inputs'] = {provenance.input_id: {}}
    # CTA proteins are encoded by the host; the antigen-aware self match
    # separately retains matches outside the declared CTA source exclusions.
    dataset.epitopes = tuple(replace(e, occurs_in_reference=True) for e in dataset.epitopes)
    dataset.result.df['occurs_in_reference'] = True
    dataset.selection['antigen_inputs'][provenance.input_id].update({
        d.assessment.antigen.source_identifier: asdict(d.assessment.antigen)
        for d in admission.decisions if d.assessment is not None})
    dataset.selection['cta_expression_admission'] = to_native_json(admission)
    if bundle is not None:
        dataset.selection['hitlist_evidence_bundle'] = bundle
    report_frame = dataset.report_frame()
    report_frame.attrs['cta_predictor'] = predictor
    return ExternalReport('cta_expression', provenance.path,
                          report_df=report_frame, epitopes=dataset.epitopes,
                          dataset=dataset, input_provenance=provenance)
