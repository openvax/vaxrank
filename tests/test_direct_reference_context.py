"""Direct custom-genome assertions use a declared reference, never a build-name guess."""

import json
from types import SimpleNamespace

from orfdata import ReferenceContext
import pytest

from vaxrank.input_scope import (
    InputProvenance, InputScope, load_reference_context, validate_input_scopes,
    validate_prediction_species,
)


def custom_genome(**changes):
    values = dict(reference_name='CanFam3.1', annotation_name='Ensembl', annotation_version=98)
    values.update(changes)
    return SimpleNamespace(**values)


def context_file(tmp_path):
    path = tmp_path / 'reference-context.json'
    path.write_text(json.dumps(ReferenceContext(
        9615, 'dog', 'CanFam3.1', 'Ensembl', '98').to_dict()))
    return path


def args(context):
    return SimpleNamespace(_input_species_scope={}, host_species='dog', host_taxon_id=None,
                           mhc_species=None, reference_context=context)


def test_direct_custom_genome_requires_and_accepts_explicit_context(tmp_path):
    genome = custom_genome()
    with pytest.raises(ValueError, match='provide an explicit ReferenceContext'):
        validate_prediction_species(args(None), [genome], ['DLA-88*034:01'], 'Direct prediction')
    path = context_file(tmp_path)
    scope = validate_prediction_species(args(path), [genome], ['DLA-88*034:01'], 'Direct prediction')
    assert scope == InputScope(host_taxon_id=9615, reference_assembly='CanFam3.1',
                               annotation='Ensembl:98')
    assert validate_prediction_species(args(load_reference_context(path)), [genome],
                                       ['DLA-88*034:01'], 'Junction prediction') == scope
    assert validate_prediction_species(args(path), [], ['DLA-88*034:01'],
                                       'Junction prediction') == scope


@pytest.mark.parametrize('change,reason', [
    (dict(reference_name='GRCh38'), 'reference_name'),
    (dict(annotation_name='RefSeq'), 'annotation_name'),
    (dict(annotation_version=99), 'annotation_version'),
    (dict(species='human'), 'Genome taxon'),
    (dict(taxon_id=9606), 'Genome taxon'),
])
def test_declared_context_checks_all_available_genome_metadata(tmp_path, change, reason):
    with pytest.raises(ValueError, match=reason):
        validate_prediction_species(args(context_file(tmp_path)), [custom_genome(**change)],
                                    ['DLA-88*034:01'], 'Direct prediction')


def test_context_rejects_conflicting_host_and_invalid_json(tmp_path):
    path = context_file(tmp_path)
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        validate_prediction_species(replace_args(args(path), host_species='human'),
                                    [custom_genome()], ['DLA-88*034:01'], 'Direct prediction')
    record = json.loads(path.read_text())
    record['invented_taxon'] = 9615
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='serialized ORFdata ReferenceContext'):
        load_reference_context(path)
    with pytest.raises(ValueError, match='requires assembly, annotation name and version'):
        load_reference_context(ReferenceContext(9615, 'dog', 'CanFam3.1', 'Ensembl'))


def replace_args(namespace, **changes):
    return SimpleNamespace(**{**vars(namespace), **changes})


def test_mixed_direct_scope_uses_context_for_retained_custom_genome(tmp_path):
    context = load_reference_context(context_file(tmp_path))
    scope = InputScope(patient_id='dog-1', reference_assembly='CanFam3.1',
                       annotation='Ensembl:98', mhc_alleles=('DLA-88*034:01',), host_taxon_id=9615)
    direct = InputProvenance('direct', 'vcf_bam', 'variants.json', '', scope)
    assert validate_input_scopes([direct], prediction_alleles=('DLA-88*034:01',),
                                 retained_genomes={'direct': (custom_genome(),)},
                                 reference_context=context) == scope
    with pytest.raises(ValueError, match='annotation_version'):
        validate_input_scopes([direct], prediction_alleles=('DLA-88*034:01',),
                              retained_genomes={'direct': (custom_genome(annotation_version=99),)},
                              reference_context=context)
