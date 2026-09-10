"""Actionable allele checks at prediction input boundaries, not rendering."""

from functools import lru_cache

from mhcgnomes import Allele, Pair, ParseError, parse
from topiary.ranking import KIND_ALIASES

from . import cells
from .allele_evidence import is_allele_scoped_kind


@lru_cache(maxsize=4096)
def parse_prediction_allele(allele):
    return parse(allele, infer_class2_pairing=False,
                 required_result_types=(Allele, Pair)).to_string()


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
