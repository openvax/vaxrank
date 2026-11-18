"""Adapt direct variant/RNA prediction to shared, reloadable input evidence."""
from dataclasses import replace
from pathlib import Path

from topiary import TopiaryResult, combine_sources

from .epitope_dataset import EpitopeDataset
from .epitope_dsl import epitopes_to_topiary_df
from .epitope_io import normalize_hla_allele
from .file_digest import file_sha256
from .external_report import ExternalReport
from .input_scope import (
    InputProvenance, combine_scopes, normalize_alleles, provenance_columns,
    read_input_manifest, resolved_reference_assembly, scope_from_mapping,
    scope_dict, validate_input_scopes,
    species_scope_from_args, with_species_declarations, load_reference_context,
)


def has_direct_input(args):
    return any(getattr(args, name, None) for name in ('vcf', 'maf', 'variant', 'json_variants', 'bam'))


def direct_input_paths(args):
    """File-backed variant and RNA evidence; inline alleles are hashed separately."""
    paths = []
    for name in ('vcf', 'maf', 'json_variants', 'germline_vcf', 'bam'):
        value = getattr(args, name, None)
        paths.extend([value] if isinstance(value, str) else value or ())
    return list(dict.fromkeys(paths))


def direct_content_identity(args, variants):
    from .selection_policy import canonical_digest
    from .native_serialization import to_native_json
    hashes = {str(Path(path).resolve()): file_sha256(path) for path in direct_input_paths(args)}
    # File bytes include metadata; normalized variants also identify inline input
    # and the actual alleles/reference used by every direct input format.
    content = dict(file_sha256=list(hashes.values()), variants=[to_native_json(v) for v in variants])
    return canonical_digest(content), hashes


def direct_input_scope(args, epitope_config, *, pooled=False):
    """Resolve actual variant/RNA scope before any source starts prediction."""
    from .epitope_config import EpitopeConfig
    epitope_config = epitope_config or EpitopeConfig()
    from mhctools.cli import mhc_alleles_from_args
    from varcode.cli import variant_collection_from_args

    manifest = getattr(args, 'input_manifest', None)
    if pooled and not manifest:
        raise ValueError('Combining direct variants/RNA with other sources requires --input-manifest '
                         'to declare their shared patient, reference_assembly and mhc_alleles')
    declarations = read_input_manifest(manifest, include_direct=True)[1] if manifest else {
        'patient_id': getattr(args, 'output_patient_id', None) or None}
    declarations = with_species_declarations(declarations, species_scope_from_args(args), 'Direct variant/RNA input')
    scope = scope_from_mapping(declarations, 'Direct variant/RNA input')
    from .patient_background import direct_patient_background, direct_rna_scope
    background = direct_patient_background(args, epitope_config)
    scope = combine_scopes(scope, direct_rna_scope(args, background), 'Direct RNA')
    alleles = normalize_alleles(mhc_alleles_from_args(args), 'Direct prediction')
    if not alleles:
        raise ValueError('Direct variant/RNA predictions require explicit MHC alleles')
    if scope.mhc_alleles is None and not manifest:
        scope = combine_scopes(scope, scope_from_mapping({'mhc_alleles': alleles}, 'Direct alleles'), 'Direct alleles')
    variants = variant_collection_from_args(args)
    context = load_reference_context(getattr(args, 'reference_context', None))
    if context is not None:
        scope = combine_scopes(scope, scope_from_mapping({
            'host_taxon_id': context.taxon_id, 'reference_assembly': context.genome_build,
            'annotation': f'{context.annotation_name}:{context.annotation_version}',
        }, 'Direct reference_context'), 'Direct reference_context')
    # Validate actual variant reference resolution as well as explicit CLI options.
    for variant in variants:
        genome = variant.genome
        scope = combine_scopes(scope, scope_from_mapping({
            'reference_assembly': resolved_reference_assembly(genome),
            'annotation': 'ensembl:' + str(genome.release) if hasattr(genome, 'release') else None,
        }, 'Direct variant reference'), 'Direct variant reference')
    paths = direct_input_paths(args)
    # Check scope before hashing large BAMs, instantiating predictors, or outputs.
    provisional = InputProvenance(
        input_id='direct', source_format='vcf_bam', path='; '.join(paths) or 'inline variants', scope=scope,
        content_sha256='', report_declarations={'reference_context': context.to_dict()} if context else {})
    validate_input_scopes([provisional], prediction_alleles=alleles,
        retained_genomes={'direct': tuple({id(v.genome): v.genome for v in variants}.values())},
        reference_context=context)
    return provisional, variants, context, declarations, paths, alleles


def prepare_direct_report(args, reports, epitope_config, *, prepared_input=None):
    """Validate input scope, then predict only the direct candidates."""
    from .cli.entry_point import run_vaxrank_from_parsed_args
    from .vaf import extract_dna_vaf_by_variant
    provisional, variants, context, declarations, paths, alleles = (
        prepared_input or direct_input_scope(args, epitope_config, pooled=bool(reports)))
    scope = provisional.scope
    manifest = getattr(args, 'input_manifest', None)
    shared = validate_input_scopes([provisional] + [r.input_provenance for r in reports],
                          prediction_alleles=alleles,
                          retained_genomes={'direct': tuple({id(v.genome): v.genome for v in variants}.values())},
                          reference_context=context)
    args._input_species_scope = scope_dict(shared)
    from .selection_policy import canonical_digest
    content_hash, hashes = direct_content_identity(args, variants)
    provenance = replace(
        provisional, input_id='vcf_bam:' + content_hash + ':' + canonical_digest(scope_dict(scope)),
        content_sha256=content_hash, manifest_path=str(Path(manifest).resolve()) if manifest else None,
        manifest_declarations=declarations,
        report_declarations=dict(provisional.report_declarations, file_sha256=hashes))
    dataset = EpitopeDataset.from_predictions((), config=epitope_config)
    dataset.provenance = (provenance,)
    results = run_vaxrank_from_parsed_args(
        args, epitope_dataset=dataset, epitope_config_override=epitope_config, loaded_variants=variants)
    dna_vaf = extract_dna_vaf_by_variant(
        variants, tumor_sample_name=getattr(args, 'tumor_sample_name', None))
    from .construct_replay import capture_direct_sources
    capture_direct_sources(dataset, results, dna_vaf_by_variant=dna_vaf)
    dataset.mutation_fragments = {
        key: replace(fragment, dna_vaf=dna_vaf.get(fragment.variant, fragment.dna_vaf))
        for key, fragment in dataset.mutation_fragments.items()}
    def canonical_predictions(peptide):
        return tuple(replace(p, allele=normalize_hla_allele(p.allele) if p.allele else '')
                     for p in peptide.predictions_flat())
    observed = sorted({normalize_hla_allele(p.allele)
                       for e in dataset.epitopes for p in e.predictions_flat() if p.allele})
    provenance = replace(provenance, observed_mhc_alleles=tuple(observed))
    dataset.provenance = (provenance,)
    dataset.epitopes = tuple(replace(
        e, predictions=canonical_predictions(e), input_provenance=provenance,
        comparators={name: replace(c, predictions=canonical_predictions(c))
                     for name, c in e.comparators.items()},
        patient_alleles=tuple(normalize_hla_allele(a) for a in e.patient_alleles),
        per_allele_scores={}, allele_attributions=()) for e in dataset.epitopes)
    frame = epitopes_to_topiary_df(dataset.epitopes)
    if frame.empty:
        frame = dataset.result.df.copy()
    fields = ('n_rna_alt', 'n_rna_ref', 'n_rna_overlapping',
              'n_rna_supporting_protein_sequence', 'rna_evidence_method',
              'rna_evidence_subject', 'gene_name', 'sequence_source',
              'sequence_source_version', 'dna_vaf')
    for field in fields:
        frame[field] = frame.prediction_id.map(
            {key: getattr(fragment, field) for key, fragment in dataset.mutation_fragments.items()})
    for name, value in provenance_columns(provenance).items():
        frame[name] = value
    for field in ('gene_id', 'transcript_ids', 'protein_ids', 'species', 'source_identifier'):
        frame[field] = frame.prediction_id.map(
            {key: getattr(antigen, field) for key, antigen in dataset.antigens.items()})
    frame['input_source'] = provenance.input_id
    frame['input_format'] = 'vcf_bam'
    frame['input_path'] = provenance.path
    dataset.result = combine_sources({provenance.input_id: TopiaryResult(frame)},
                                     sample_name=scope.sample_id or scope.patient_id or provenance.input_id)
    identities = dataset.result.df[['prediction_id', 'source_observation_id']].drop_duplicates()
    if identities.prediction_id.duplicated().any():
        raise ValueError('Direct prediction rows disagree about source occurrence identity')
    ids = dict(identities.itertuples(index=False, name=None))
    dataset.epitopes = tuple(replace(e, prediction_id=ids[e.prediction_id]) for e in dataset.epitopes)
    dataset.antigens = {ids[key]: value for key, value in dataset.antigens.items()}
    dataset.mutation_fragments = {ids[key]: value for key, value in dataset.mutation_fragments.items()}
    return ExternalReport('vcf_bam', provenance.path, report_df=dataset.report_frame(),
                          epitopes=dataset.epitopes, dataset=dataset,
                          input_provenance=provenance)
