# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Validation of amino-acid sequences used in vaccine safety models."""

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def has_only_standard_amino_acids(sequence: str) -> bool:
    """Return whether every residue belongs to the standard amino-acid set."""
    return all(residue in STANDARD_AMINO_ACIDS for residue in sequence)


def validate_amino_acid_sequence(sequence: str, label: str) -> None:
    """Reject empty sequences and residues outside the standard alphabet."""
    if not sequence:
        raise ValueError(f"{label} amino-acid sequence is required")
    invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
    if invalid:
        raise ValueError(
            f"{label} contains non-canonical amino acids: {', '.join(invalid)}"
        )
