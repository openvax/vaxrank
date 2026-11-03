"""Host defaults must not be inferred from antigen origin or MHC repertoire."""
from types import SimpleNamespace

import pytest
from pyensembl import EnsemblRelease

from vaxrank.input_scope import InputScope, host_scope_from_args, resolve_host_scope, validate_prediction_species
from vaxrank.mrna import RNAConstructConfig, resolve_codon_species
from vaxrank.cta_input import cta_expression_contract, CTAExpressionSource


@pytest.mark.parametrize('species,taxon,codon', [('human', 9606, 'h_sapiens'), ('mouse', 10090, 'm_musculus')])
def test_reference_resolves_host_and_codon_default(species, taxon, codon):
    genome = EnsemblRelease(102, species=species)
    scope = validate_prediction_species(SimpleNamespace(), [genome], [], 'test')
    assert scope.host_taxon_id == taxon
    assert resolve_codon_species(RNAConstructConfig().codon_species, scope.host_taxon_id) == codon


def test_host_declaration_and_saved_reference_use_same_resolution():
    assert host_scope_from_args(SimpleNamespace(host_species='mouse')).host_taxon_id == 10090
    args = SimpleNamespace(_input_species_scope={'reference_assembly': 'GRCm38'})
    assert host_scope_from_args(args).host_taxon_id == 10090
    args.host_species = 'human'
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        host_scope_from_args(args)


def test_mhc_repertoire_and_source_species_do_not_infer_host():
    args = SimpleNamespace(mhc_species=['human'], source_species='human')
    assert host_scope_from_args(args).host_taxon_id is None
    with pytest.raises(ValueError, match='resolved host species'):
        resolve_codon_species()


def test_transgenic_mouse_keeps_mouse_codons_and_explicit_target_wins():
    args = SimpleNamespace(host_species='mouse', mhc_species=['human'])
    scope = validate_prediction_species(args, [EnsemblRelease(102, species='mouse')], ['HLA-A*02:01'], 'test')
    assert resolve_codon_species(host_taxon_id=scope.host_taxon_id) == 'm_musculus'
    assert resolve_codon_species('human', scope.host_taxon_id) == 'h_sapiens'
    with pytest.raises(ValueError, match='No bundled codon table'):
        resolve_codon_species(host_taxon_id=9615)


def test_conflicting_reference_hosts_fail():
    args = SimpleNamespace()
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        host_scope_from_args(args, [EnsemblRelease(102, species='mouse'), EnsemblRelease(102)])


def cta_options():
    return dict(measurement_level='gene', id_column='id', value_column='TPM',
                expression_unit='TPM', evidence_source='fixture', evidence_version='1',
                assay='RNA-seq', min_expression=1)


def test_cta_defaults_wait_for_reference_and_nonhuman_needs_mapped_evidence(tmp_path):
    unknown = InputScope(sample_id='sample')
    assert cta_expression_contract(cta_options(), unknown, 'test')[2] is None
    human = resolve_host_scope(unknown, EnsemblRelease(102))
    mouse = resolve_host_scope(unknown, EnsemblRelease(102, species='mouse'))
    assert cta_expression_contract(cta_options(), human, 'test')[2].gene_patterns
    assert not cta_expression_contract(cta_options(), mouse, 'test')[2].gene_patterns
    with pytest.raises(ValueError, match='mapped host_reference'):
        CTAExpressionSource.load(tmp_path / 'missing.tsv', options=cta_options(),
                                 scope=unknown, genome=EnsemblRelease(102, species='mouse'))


@pytest.mark.parametrize('host', [None, 'mouse', 'dog'])
def test_human_only_processing_rejects_unknown_and_nonhuman_hosts(host):
    from vaxrank.input_scope import validate_processing_host
    args = SimpleNamespace(host_species=host, pepsickle_human_only=True)
    with pytest.raises(ValueError, match='resolved human host'):
        validate_processing_host(args)
    args.pepsickle_human_only = False
    validate_processing_host(args)
    args.host_species, args.pepsickle_human_only = 'human', True
    validate_processing_host(args)
