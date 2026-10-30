"""Registered source capabilities and total budgets are checked before transfer."""

import oncoref
import pytest

from vaxrank.cta_reference_adapter import _reference_context, build_cta_host_reference

def test_owner_reference_binds_the_matching_variant_annotation_without_indexing(monkeypatch):
    from orfdata import Proteome
    from pyensembl import EnsemblRelease
    from vaxrank.host_reference import HostReference
    from vaxrank.input_scope import InputScope

    genome = EnsemblRelease(115)
    def forbidden(*args, **kwargs):
        raise AssertionError('Reference metadata tried to download or index an annotation')
    monkeypatch.setattr(type(genome), 'download', forbidden)
    monkeypatch.setattr(type(genome), 'index', forbidden)
    identity = oncoref.human.cta_snapshot().as_dict()['references'][0]['identity']
    context = _reference_context(identity, 'homo_sapiens')
    # No proteins are needed to exercise the scope/annotation binding boundary.
    host = HostReference(Proteome(context, []))
    scope = InputScope(reference_assembly=genome.assembly, annotation=f'ensembl:{genome.release}')
    assert host.validate(scope, genome=genome).reference_assembly == 'GRCh38'
    with pytest.raises(ValueError, match='evidence annotation'):
        _reference_context({**identity, 'assembly_accession': 'GCA_000001405.28'}, 'homo_sapiens')


@pytest.fixture
def no_acquisition(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('A rejected source request tried to acquire an asset')
    monkeypatch.setattr('oncoref.species_assets.ensure_asset', forbidden)


@pytest.mark.parametrize('species,budget', [('human', 120_000_000), ('dog', 10_000_000)])
def test_total_missing_source_budget_is_checked_before_any_download(tmp_path, no_acquisition, species, budget):
    # Each required individual file fits; the complete source set does not.
    with pytest.raises(ValueError, match='needs .* bytes; transfer budget'):
        build_cta_host_reference(species, cache_root=tmp_path, local_overrides={}, max_bytes=budget)


def test_mouse_candidate_seed_is_rejected_as_incomplete_host_before_acquisition(tmp_path, no_acquisition):
    with pytest.raises(ValueError, match='candidate subsets are insufficient'):
        build_cta_host_reference('mouse', cache_root=tmp_path, local_overrides={}, max_bytes=1_000_000)


@pytest.mark.parametrize('budget', [True, -1, 1.5, None])
def test_builder_requires_explicit_integral_nonnegative_budget(tmp_path, no_acquisition, budget):
    with pytest.raises(ValueError, match='nonnegative integer'):
        build_cta_host_reference('human', cache_root=tmp_path, max_bytes=budget)
