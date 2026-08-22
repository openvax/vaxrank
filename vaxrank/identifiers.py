# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Canonical forms for biological identifiers used across vaxrank."""

from typing import Optional

from mhcgnomes import parse as parse_mhc


def normalize_ensembl_gene_id(value: Optional[str]) -> str:
    """Strip whitespace and an optional Ensembl stable-ID version suffix."""
    if value is None:
        return ""
    return str(value).strip().split(".")[0]


def normalize_mhc_allele(value: str) -> str:
    """Return the mhcgnomes canonical representation of an MHC allele."""
    try:
        allele = parse_mhc(str(value), raise_on_error=True)
    except Exception as error:
        raise ValueError(f"Invalid MHC allele {value!r}") from error
    if allele is None:
        raise ValueError(f"Invalid MHC allele {value!r}")
    return allele.to_string()
