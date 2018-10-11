from vaxrank.gene_pathway_check import (
    GenePathwayCheck,
    _IFNG_RESPONSE_COLUMN_NAME,
    _CLASS_I_MHC_COLUMN_NAME,
    _DRIVER_GENE_COLUMN_NAME,
    _DRIVER_VARIANT_COLUMN_NAME
)
from varcode import Variant
from nose.tools import eq_


def test_HRAS_G13C_in_cancer_driver_genes():
    HRAS_G13C = Variant("11", 534286, "C", "A", "GRCh37")
    effect = HRAS_G13C.effects().top_priority_effect()
    eq_(effect.gene.name, "HRAS")
    eq_(effect.short_description, "p.G13C")
    gene_pathway_check = GenePathwayCheck()
    variant_info = gene_pathway_check.make_variant_dict(HRAS_G13C)
    assert not variant_info[_IFNG_RESPONSE_COLUMN_NAME]
    assert not variant_info[_CLASS_I_MHC_COLUMN_NAME]
    # even though it's a RAS G13 variant, it's not actually that common
    # and thus didn't make the threshold for our source dataset
    assert not variant_info[_DRIVER_VARIANT_COLUMN_NAME]
    assert variant_info[_DRIVER_GENE_COLUMN_NAME]


def test_HRAS_G13V_in_cancer_driver_genes_and_variants():
    HRAS_G13V = Variant("11", 534285, "C", "A", "GRCh37")
    effect = HRAS_G13V.effects().top_priority_effect()
    eq_(effect.gene.name, "HRAS")
    eq_(effect.short_description, "p.G13V")
    gene_pathway_check = GenePathwayCheck()
    variant_info = gene_pathway_check.make_variant_dict(HRAS_G13V)
    assert not variant_info[_IFNG_RESPONSE_COLUMN_NAME]
    assert not variant_info[_CLASS_I_MHC_COLUMN_NAME]
    assert variant_info[_DRIVER_VARIANT_COLUMN_NAME]
    assert variant_info[_DRIVER_GENE_COLUMN_NAME]
