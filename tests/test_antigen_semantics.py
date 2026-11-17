"""Producer-independent biological facts and actual prediction routing."""
from types import SimpleNamespace

import pandas as pd
import pytest
from pyensembl import EnsemblRelease
from topiary import TopiaryPredictor
from varcode import StructuralVariant

from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_logic import predict_epitopes
from vaxrank.mutant_protein_fragment import MutantProteinFragment


@pytest.mark.parametrize('sv_type', ['DEL', 'DUP', 'INV', 'INS', 'BND'])
@pytest.mark.parametrize('reference', ['SIINFEKA', None])
def test_direct_sv_prediction_keeps_variant_overlap_and_supported_comparators(
        monkeypatch, sv_type, reference):
    variant = StructuralVariant('1', 100, sv_type, end=200, genome=EnsemblRelease(110),
        **({'mate_contig': '2', 'mate_start': 300} if sv_type == 'BND' else {}))
    fragment = MutantProteinFragment(variant=variant, gene_name='SV', amino_acids='AASIINFEKLKK',
        mutant_amino_acid_start_offset=5, mutant_amino_acid_end_offset=6,
        supporting_reference_transcripts=[], n_overlapping_reads=20, n_alt_reads=10,
        n_ref_reads=10, n_alt_reads_supporting_protein_sequence=10)
    lookups, predictions = [], []
    def comparator(self, start, length):
        lookups.append((start, length))
        return reference
    monkeypatch.setattr(MutantProteinFragment, 'wildtype_peptide_at', comparator)
    monkeypatch.setattr('vaxrank.epitope_logic.ReferenceProteome', type('Reference', (), {
        '__init__': lambda self, *a, **k: None, 'contains': lambda self, p: False,
        'from_genome': staticmethod(lambda *a, **k: SimpleNamespace(contains=lambda p: False))}))
    monkeypatch.setattr('vaxrank.epitope_logic.self_reference_matches',
        lambda peptides, antigen, genome: {p: antigen.self_reference_match(p, False) for p in peptides})
    class Predictor(TopiaryPredictor):
        def __init__(self):
            pass
        def frame(self, named, offset):
            return pd.DataFrame([dict(source_sequence_name=name, peptide=peptide,
                peptide_offset=offset, peptide_length=len(peptide), allele='HLA-A*02:01',
                prediction_method_name='random', predictor_version='', kind='pMHC_affinity',
                affinity=20., value=20., score=1., percentile_rank=.1) for name, peptide in named.items()])
        def predict_from_named_sequences(self, named):
            return self.frame({n: s[2:10] for n, s in named.items()}, 2)
        def predict_from_named_peptides(self, named):
            predictions.append(named)
            return self.frame(named, 0)
    epitopes = predict_epitopes(Predictor(), protein_fragment=fragment,
        genome=EnsemblRelease(110), epitope_config=EpitopeConfig(score_expr='1.0'))
    assert len(epitopes) == 1
    epitope, = epitopes
    assert epitope.source_class == 'mutation' and epitope.overlaps_mutation and epitope.is_targetable
    assert lookups == [(2, 8)]
    assert bool(predictions) == (reference is not None)
    assert (epitope.reference_comparator.sequence if epitope.reference_comparator else None) == reference


def test_lens_and_isovar_share_junction_fields_without_inventing_dna_or_translation():
    from vaxrank.antigen_origin import origin_from_lens, origin_from_isovar
    lens = origin_from_lens('fusion', {'left_breakpoint': 'chr1:100:+', 'right_breakpoint': 'chr2:200:-'})
    isovar = origin_from_isovar({'donor': {'contig': 'chr1', 'position': 100, 'strand': '+'},
        'acceptor': {'contig': 'chr2', 'position': 200, 'strand': '-'}}, 'fusion', {})
    assert lens == isovar
    assert not lens.variants and not lens.junctions[0].protein_intervals
    splice = origin_from_lens('splice', {'splice_coords': 'chr1:10-20(+)'})
    assert splice.rna_product == 'splice' and not splice.variants
    assert not splice.junctions[0].protein_intervals


@pytest.mark.parametrize('ref,alt,expected', [('A', 'T', 'SNV'), ('AA', 'TC', 'MNV'),
    ('AA', 'A', 'INDEL'), ('A', 'ATG', 'INDEL')])
def test_direct_and_isovar_alleles_have_identical_classification(ref, alt, expected):
    from varcode import Variant
    from vaxrank.antigen_origin import VariantOrigin
    variant = Variant('1', 100, ref, alt, genome=EnsemblRelease(110))
    direct = VariantOrigin.from_variant(variant)
    producer = VariantOrigin.from_alleles(dict(contig=variant.contig, start=variant.start,
        ref=variant.ref, alt=variant.alt, reference_name=variant.genome.reference_name))
    assert direct == producer and direct.variant_type == expected


def test_native_origin_roundtrip_and_slicing(tmp_path):
    from dataclasses import replace
    from vaxrank.antigen_input import AntigenDataset
    from vaxrank.native_serialization import to_native_json, from_native_json
    from vaxrank.vaccine_antigen import VaccineAntigen
    from vaxrank.antigen_origin import origin_from_lens
    from .test_assembly_antigens import candidate
    antigen = replace(candidate(0)[0], kind='FUSION',
        origin=origin_from_lens('fusion', {'left_breakpoint': 'chr1:100:+'}))
    assert from_native_json(to_native_json(antigen), VaccineAntigen) == antigen
    assert antigen.sliced(1, len(antigen.amino_acids)).origin == antigen.origin
    path = tmp_path / 'input.json'
    AntigenDataset({'fusion': antigen}).save(path)
    assert AntigenDataset.load(path).antigens['fusion'] == antigen


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_variant_subtype_and_context_filters_apply_to_both_products(modality):
    from dataclasses import replace
    from vaxrank.antigen_origin import AntigenOrigin, CodingReference, VariantOrigin
    from vaxrank.antigen_allocation import AntigenAllocation, select_ranked_sources
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
    from .test_assembly_antigens import candidate
    def entry(index, variant, context, sequence):
        antigen, peptides = candidate(index, sequence)
        antigen = replace(antigen, origin=AntigenOrigin(variants=(variant,), coding_references=(context,)))
        return antigen, [replace(peptides[0], antigen=antigen)]
    snv = entry(0, VariantOrigin('SNV'), CodingReference(kind='extra_orf', orf_type='uORF'), 'SIINFEKLA')
    indel = entry(1, VariantOrigin('INDEL'), CodingReference(kind='annotated'), 'GILGFVFTL')
    deletion = entry(2, VariantOrigin('SV', sv_type='DEL'), CodingReference(kind='annotated'), 'NLVPMVATV')
    policy = AntigenAllocation(variant_types=('snv', 'sv:del'), variant_contexts=('uORF',))
    ranked = [snv, indel, deletion]
    audit = {}
    assert select_ranked_sources(ranked, policy, audit=audit) == [snv]
    assert [row['reason'] for row in audit['ranked_source_selection']] == [
        'selected', 'variant_filter_mismatch', 'variant_filter_mismatch']
    assert select_ranked_sources(ranked, AntigenAllocation(variant_types=('sv:del',))) == [deletion]
    options = dict(min_antigen_length_aa=8, max_antigen_length_aa=13, antigen_allocation=policy)
    if modality == 'peptide':
        products = assemble_peptide_constructs(ranked, PeptideConstructConfig(**options))
    else:
        products = assemble_mrna_constructs(ranked, RNAConstructConfig(**options, optimize_linkers=False,
            codon_method='use_best_codon'), host_taxon_id=9606)
    assert len(products) == 1
    source, = products[0].antigens[0]['sources']
    assert source['origin']['variants'][0]['variant_type'] == 'SNV'
    assert source['coding_contexts'][0]['orf_type'] == 'uORF'


def test_selectors_normalize_cli_and_configuration_consistently():
    from vaxrank.cli.arg_parser import external_input_arg_parser
    from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
    from vaxrank.antigen_allocation import AntigenAllocation
    from vaxrank.vaccine_config import VaccineConfig
    args = external_input_arg_parser().parse_args(['--antigen-types', 'VARIANT', 'cta',
        '--variant-types', 'snv', 'Sv:del', '--variant-contexts', 'uORF', 'ANNOTATED',
        '--include-lens-source', 'splice'])
    policy = AntigenAllocation.resolve(vaccine_config_from_args(args).antigen_allocation)
    assert not any(policy.max_per_kind.get(kind) == 0 for kind in ('mutation', 'ncORF', 'SV', 'CTA'))
    assert policy.variant_types == ('SNV', 'SV:DEL')
    assert policy.variant_contexts == ('annotated', 'uorf')
    assert VaccineConfig(included_lens_sources=('snv', 'splice')).included_lens_sources == ('SNV', 'SPLICE')
    assert AntigenAllocation(max_per_kind={'sv': 0, 'cta': 2}).max_per_kind == {'SV': 0, 'CTA': 2}


def test_isovar_ncorf_retains_independent_event_and_context_fields():
    from vaxrank.isovar_antigens import antigens_from_isovar
    from vaxrank.antigen_allocation import AntigenAllocation, variant_enabled, selection_type
    from vaxrank.antigen_biology import antigen_biology
    from .test_isovar_inputs import result, attestation
    antigen, = antigens_from_isovar(result('ncorf'), tumor_specificity=attestation()).antigens.values()
    assert antigen.origin.variants[0].variant_type == 'SNV'
    reference, = antigen.origin.coding_references
    assert reference.kind == 'extra_orf' and reference.orf_id == 'noncoding_orf'
    assert selection_type(antigen) == 'ncORF'
    assert variant_enabled(antigen, AntigenAllocation(variant_types=('snv',), variant_contexts=('ncORF',)))
    assert not variant_enabled(antigen, AntigenAllocation(variant_types=('indel',)))
    assert dict(antigen_biology(antigen)['facts'])['Variant type'] == 'SNV'


def test_lens_real_inputs_populate_shared_origin():
    from pathlib import Path
    from vaxrank.external_input import lens_ranking_result, ExternalConstructOptions
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.antigen_biology import antigen_biology
    from vaxrank.vaccine_antigen import VaccineAntigen
    report = read_lens_report(Path(__file__).parent / 'data/epitope_fixtures/real_lens_subsets/lens_v1.9_real_subset.tsv')
    ranked = lens_ranking_result(report, list(report.epitopes), options=ExternalConstructOptions(
        included_lens_sources=('SNV', 'INDEL', 'FUSION', 'SPLICE'))).ranked
    fusion = next(source for source, _ in ranked if isinstance(source, VaccineAntigen) and source.kind == 'fusion')
    assert fusion.origin.rna_product == 'fusion' and not fusion.origin.variants
    assert fusion.origin.junctions[0].left.contig
    assert antigen_biology(fusion)['junctions'][0]['left'] == fusion.origin.junctions[0].left.label()
    # DNA-backed LENS entries use the same Varcode adapter as direct VCF input.
    mutations = [p.antigen for _, peptides in ranked for p in peptides if p.antigen.kind == 'mutation']
    assert mutations and all(a.origin and a.origin.variants for a in mutations)


def test_native_origin_rejects_unknown_nested_fields(tmp_path):
    from dataclasses import replace
    import json
    from vaxrank.antigen_input import AntigenDataset
    from vaxrank.antigen_origin import AntigenOrigin, VariantOrigin
    from .test_assembly_antigens import candidate
    antigen = replace(candidate(0)[0], origin=AntigenOrigin(variants=(VariantOrigin('SNV'),)))
    data = AntigenDataset({'test': antigen}).to_dict()
    data['antigens']['test']['origin']['variants'][0]['made_up'] = True
    path = tmp_path / 'invalid.json'
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='Unknown VariantOrigin fields'):
        AntigenDataset.load(path)


@pytest.mark.parametrize('ref,alt', [('A', 'T'), ('A', 'AT'), ('AT', 'A')])
def test_product_event_identity_is_shared_by_direct_and_isovar_origins(ref, alt):
    from dataclasses import replace
    from varcode import Variant
    from vaxrank.antigen_origin import AntigenOrigin, VariantOrigin
    from vaxrank.assembly_antigens import summarize_products
    from vaxrank.peptide import assemble_peptide_constructs, PeptideConstructConfig
    from .test_assembly_antigens import candidate
    variant = Variant('1', 100, ref, alt, genome=EnsemblRelease(110))
    direct = VariantOrigin.from_variant(variant)
    exported = VariantOrigin.from_alleles(dict(contig='1', start=variant.start, ref=variant.ref, alt=variant.alt, reference_name='GRCh38'))
    ranked = []
    for index, origin in enumerate((direct, exported)):
        antigen, peptides = candidate(index)
        antigen = replace(antigen, origin=AntigenOrigin(variants=(origin,)))
        ranked.append((variant if index == 0 else antigen, [replace(peptides[0], antigen=antigen)]))
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    assert len(products) == 1
    sources = products[0].antigens[0]['sources']
    assert len(sources) == 2 and sources[0]['source_event'] == sources[1]['source_event']
    assert sources[0]['source_event']['ref'] == variant.ref
    assert sources[0]['source_event']['alt'] == variant.alt
    assert summarize_products(products)['source_events'] == 1


@pytest.mark.parametrize('recorded_origin', [False, True])
def test_native_fragment_ranking_preserves_saved_origin(tmp_path, recorded_origin):
    from dataclasses import replace
    from varcode import Variant
    from vaxrank.antigen_origin import AntigenOrigin, CodingReference, VariantOrigin
    from vaxrank.epitope_dataset import EpitopeDataset, dataset_ranking_result
    from vaxrank.external_input import ExternalConstructOptions
    from vaxrank.vaccine_config import VaccineConfig
    from vaxrank.vaccine_antigen import VaccineAntigen
    from .test_core_logic_config import _make_epitope
    variant = Variant('1', 100, 'A', 'T', genome=EnsemblRelease(110))
    fragment = MutantProteinFragment(variant=variant, gene_name='RNA', amino_acids='AASIINFEKLKK',
        mutant_amino_acid_start_offset=5, mutant_amino_acid_end_offset=6,
        supporting_reference_transcripts=[], n_overlapping_reads=20, n_alt_reads=10,
        n_ref_reads=10, n_alt_reads_supporting_protein_sequence=10)
    origin = AntigenOrigin(variants=(VariantOrigin.from_variant(variant),),
        coding_references=(CodingReference(kind='extra_orf', orf_type='uORF'),)) if recorded_origin else None
    antigen = replace(VaccineAntigen.from_mutant_protein_fragment(fragment), origin=origin)
    data = EpitopeDataset.from_predictions((), config=EpitopeConfig(score_expr='1.0'))
    from vaxrank.input_scope import InputProvenance, InputScope
    data.provenance = (InputProvenance('direct', 'vcf_bam', 'fixture.vcf', 'test', InputScope()),)
    data.add_mutation(fragment, [_make_epitope('SIINFEKL', 20., offset=2, source_sequence=fragment.amino_acids)])
    data.antigens = {key: antigen for key in data.antigens}
    path = tmp_path / 'recorded.tsv'
    data.save(path)
    saved = EpitopeDataset.load(path)
    options = ExternalConstructOptions.from_configs(VaccineConfig(preferred_peptide_length=10,
        min_peptide_length=8, max_peptide_length=12, combined_score_expr='target_epitope_score'))
    outcome = dataset_ranking_result(SimpleNamespace(dataset=saved, source_format='epitopes'), saved.epitopes, options=options)
    assert outcome.ranked
    for _, peptides in outcome.ranked:
        assert all(p.antigen.origin == origin for p in peptides)
        assert all(p.amino_acids == antigen.amino_acids[p.antigen.source_protein_offset:
            p.antigen.source_protein_offset+len(p.amino_acids)] for p in peptides)
    saved.antigens = {key: replace(value, tumor_specificity=replace(value.tumor_specificity,
        status='held_out', rationale_code='not_tumor_specific')) for key, value in saved.antigens.items()}
    assert not dataset_ranking_result(SimpleNamespace(dataset=saved, source_format='epitopes'),
        saved.epitopes, options=options).ranked


@pytest.mark.parametrize('start,mate', [(0, 100), (100, 0)])
def test_origin_retains_vcf_telomeric_breakend_coordinates(start, mate):
    from vaxrank.antigen_origin import VariantOrigin, GenomicEndpoint
    variant = StructuralVariant('1', start, 'BND', genome=EnsemblRelease(110),
        mate_contig='2', mate_start=mate, mate_orientation='++')
    origin = VariantOrigin.from_variant(variant)
    assert (origin.start, origin.mate_start) == (start, mate)
    assert GenomicEndpoint.from_record({'contig': '1', 'position': 0, 'strand': '+'}).position == 0
