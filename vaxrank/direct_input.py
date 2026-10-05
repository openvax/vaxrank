"""Adapt the existing VCF/BAM pipeline to shared, reloadable input evidence."""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

from topiary import TopiaryResult, combine_sources

from .epitope_dataset import EpitopeDataset
from .epitope_dsl import epitopes_to_topiary_df
from .epitope_io import normalize_hla_allele
from .external_report import ExternalReport
from .input_scope import (
    InputProvenance, combine_scopes, normalize_alleles, provenance_columns,
    read_input_manifest, scope_from_mapping, validate_input_scopes,
)


def _file_digest(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_direct_report(args, reports, epitope_config):
    """Validate pooled scope, then predict only the direct candidates."""
    from mhctools.cli import mhc_alleles_from_args
    from varcode.cli import variant_collection_from_args
    from .cli.entry_point import run_vaxrank_from_parsed_args
    from .vaf import extract_dna_vaf_by_variant

    manifest = getattr(args, 'input_manifest', None)
    if not manifest:
        raise ValueError('Combining VCF/BAM with tables requires --input-manifest '
                         'to declare their shared patient, reference_assembly and mhc_alleles')
    _, declarations = read_input_manifest(manifest, include_direct=True)
    scope = scope_from_mapping(declarations, 'Direct VCF/BAM input')
    alleles = normalize_alleles(mhc_alleles_from_args(args), 'Direct prediction')
    if not alleles:
        raise ValueError('Direct VCF/BAM predictions require explicit MHC alleles')
    variants = variant_collection_from_args(args)
    # Validate actual VCF reference resolution as well as explicit CLI options.
    for variant in variants:
        genome = variant.ensembl
        scope = combine_scopes(scope, scope_from_mapping({
            'reference_assembly': genome.reference_name,
            'annotation': 'ensembl:' + str(genome.release) if hasattr(genome, 'release') else None,
        }, 'Direct VCF reference'), 'Direct VCF reference')
    paths = list(getattr(args, 'vcf', None) or [])
    if not paths:
        raise ValueError('Combined direct inputs require --vcf')
    if getattr(args, 'bam', None):
        paths.append(args.bam)
    paths.extend(getattr(args, 'germline_vcf', None) or [])
    # Check scope before hashing large BAMs, instantiating predictors, or outputs.
    provisional = InputProvenance(
        input_id='direct', source_format='vcf_bam', path='; '.join(paths), scope=scope,
        content_sha256='')
    validate_input_scopes([provisional] + [r.input_provenance for r in reports],
                          prediction_alleles=alleles)
    hashes = {str(Path(path).resolve()): _file_digest(path) for path in paths}
    content_hash = hashlib.sha256(json.dumps(list(hashes.values())).encode()).hexdigest()
    scope_hash = hashlib.sha256(json.dumps(asdict(scope), sort_keys=True).encode()).hexdigest()
    provenance = replace(
        provisional, input_id='vcf_bam:' + content_hash + ':' + scope_hash,
        content_sha256=content_hash, manifest_path=str(Path(manifest).resolve()),
        manifest_declarations=declarations, report_declarations={'file_sha256': hashes})
    dataset = EpitopeDataset.from_predictions((), config=epitope_config)
    dataset.provenance = (provenance,)
    results = run_vaxrank_from_parsed_args(
        args, epitope_dataset=dataset, epitope_config_override=epitope_config)
    dna_vaf = extract_dna_vaf_by_variant(
        variants, tumor_sample_name=getattr(args, 'tumor_sample_name', None))
    from .native_serialization import to_native_json
    from .gene_pathway_check import GenePathwayCheck
    gene_pathway_check = (GenePathwayCheck()
                          if getattr(args, 'output_passing_variants_csv', None) else None)
    dataset.direct_sources = [
        {'variant': to_native_json(result.variant), 'properties': properties}
        for result, properties in zip(results.isovar_results, results.variant_properties(
            dna_vaf_by_variant=dna_vaf, gene_pathway_check=gene_pathway_check))]
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
                                     sample_name=scope.patient_id)
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
