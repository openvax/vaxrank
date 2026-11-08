"""Actionable allele checks at prediction input boundaries, not rendering."""

from functools import lru_cache

from mhcgnomes import ParseError
from topiary import canonical_allele_name
from topiary.ranking import KIND_ALIASES

from . import cells
from .allele_evidence import is_allele_scoped_kind


@lru_cache(maxsize=4096)
def parse_prediction_allele(allele):
    return canonical_allele_name(allele)


def validate_prediction_allele(kind, allele, location):
    """Require an allele for MHC-scoped evidence, leaving peptide data free."""
    kind = KIND_ALIASES.get(str(kind).lower(), kind)
    text = cells.text(allele)
    if not text:
        if is_allele_scoped_kind(kind):
            raise ValueError(f"{location}: {kind} prediction requires a patient allele")
        return
    try:
        parse_prediction_allele(text)
    except (ParseError, TypeError, ValueError) as error:
        raise ValueError(f"{location}: invalid allele {text!r} for {kind}: {error}") from error


def validate_peptide_alleles(peptide, location):
    """Validate all leaves and comparator contexts without changing objects."""
    for prediction in peptide.predictions_flat():
        validate_prediction_allele(
            prediction.kind, prediction.allele,
            f"{location}, peptide {peptide.sequence!r}, predictor {prediction.predictor_name!r}")
    for name, comparator in getattr(peptide, "comparators", {}).items():
        validate_peptide_alleles(comparator, f"{location}, comparator {name!r}")


def validate_prediction_frame(frame, location):
    """Validate normalized predictor output, preserving its row provenance."""
    for row_number, (_, row) in enumerate(frame.iterrows(), start=2):
        validate_prediction_allele(
            row.get("kind"), row.get("allele"),
            f"{location}, normalized row {row_number}, "
            f"variant {cells.text(row.get('variant'))!r}, peptide {row.get('peptide')!r}")


def validate_allele_arguments(args):
    """Reject an unparseable ``--mhc-alleles`` or ``--mhc-alleles-file`` entry as typed.

    mhctools normalizes names it cannot parse (it strips the ``*`` separator),
    so by the time a later check fails the message shows a spelling the user
    never wrote. Checking the raw tokens first names the flag and the entry.
    """
    entries = [("--mhc-alleles", token)
               for group in (getattr(args, "mhc_alleles", None) or "").split(",")
               for token in group.split()]
    path = getattr(args, "mhc_alleles_file", None)
    if path:
        with open(path) as stream:
            entries += [("--mhc-alleles-file", line.strip()) for line in stream if line.strip()]
    for flag, token in entries:
        try:
            parse_prediction_allele(token)
        except (ParseError, TypeError, ValueError) as error:
            raise ValueError(f"{flag}: cannot parse {token!r} as an MHC allele") from error
