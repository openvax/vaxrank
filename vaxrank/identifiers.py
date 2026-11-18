# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Canonical forms for biological identifiers used across vaxrank."""

from typing import Optional

from topiary import canonical_allele_name


def normalize_ensembl_gene_id(value: Optional[str]) -> str:
    """Strip whitespace and an optional Ensembl stable-ID version suffix."""
    if value is None:
        return ""
    return str(value).strip().split(".")[0]


def normalize_mhc_allele(value: str) -> str:
    """Compatibility alias for Topiary's expressed-allele canonicalizer.

    Uses mhcgnomes through Topiary, matching prediction and exposure inputs.
    Genes, serotypes, allele groups and null alleles are not prediction alleles.
    """
    return canonical_allele_name(value)
