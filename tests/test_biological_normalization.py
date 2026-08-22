import pytest

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
    with pytest.raises(ValueError, match="Invalid MHC allele"):
        normalize_mhc_allele("not-an-allele")


def test_canonical_amino_acid_validation_has_boolean_and_raising_apis():
    assert has_only_standard_amino_acids("ACDEFGHIKLMNPQRSTVWY")
    assert has_only_standard_amino_acids("")
    assert not has_only_standard_amino_acids("ACDEX")
    validate_amino_acid_sequence("SIINFEKL", "Target peptide")
    with pytest.raises(ValueError, match="sequence is required"):
        validate_amino_acid_sequence("", "Target peptide")
    with pytest.raises(ValueError, match="non-canonical amino acids: J, X"):
        validate_amino_acid_sequence("AXJ", "Target peptide")
