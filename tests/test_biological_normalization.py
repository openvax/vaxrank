import pytest
from topiary import AlleleError, canonical_allele_name

from vaxrank.allele_validation import parse_prediction_allele

from vaxrank import (
    has_only_standard_amino_acids,
    normalize_ensembl_gene_id,
    normalize_mhc_allele,
    validate_amino_acid_sequence,
)


def test_normalize_ensembl_gene_id_removes_version_and_whitespace():
    assert normalize_ensembl_gene_id(" ENSG00000185686.4 ") == "ENSG00000185686"
    assert normalize_ensembl_gene_id("ENSG00000185686") == "ENSG00000185686"
    assert normalize_ensembl_gene_id("") == ""
    assert normalize_ensembl_gene_id(None) == ""


def test_normalize_mhc_allele_uses_mhcgnomes_for_human_and_mouse_alleles():
    assert normalize_mhc_allele("A*02:01") == "HLA-A*02:01"
    assert normalize_mhc_allele("H2-Kb") == "H2-K*b"
    with pytest.raises(AlleleError, match="not-an-allele"):
        normalize_mhc_allele("not-an-allele")


@pytest.mark.parametrize("allele", [
    "A*02:01", "H-2-Kb", "DLA-88*001:01", "Mamu-A*001:01",
    "HLA-DPA1*01:03-DPB1*04:01", "HLA-DRB1*08:01", "HLA-A*02:01L",
])
def test_identifier_and_prediction_share_allele_names(allele):
    expected = canonical_allele_name(allele)
    assert normalize_mhc_allele(allele) == parse_prediction_allele(allele) == expected


@pytest.mark.parametrize("allele", [
    "A2", "HLA-A", "HLA-A*02", "HLA-A*02:01N",
    "HLA-A*02:01/HLA-B*07:02", "HLA-A*01:01:01G", "", None, 201,
])
def test_identifier_and_prediction_share_allele_rejections(allele):
    for canonicalize in (canonical_allele_name, normalize_mhc_allele, parse_prediction_allele):
        with pytest.raises(AlleleError):
            canonicalize(allele)


def test_canonical_amino_acid_validation_has_boolean_and_raising_apis():
    assert has_only_standard_amino_acids("ACDEFGHIKLMNPQRSTVWY")
    assert has_only_standard_amino_acids("")
    assert not has_only_standard_amino_acids("ACDEX")
    validate_amino_acid_sequence("SIINFEKL", "Target peptide")
    with pytest.raises(ValueError, match="sequence is required"):
        validate_amino_acid_sequence("", "Target peptide")
    with pytest.raises(ValueError, match="non-canonical amino acids: J, X"):
        validate_amino_acid_sequence("AXJ", "Target peptide")
