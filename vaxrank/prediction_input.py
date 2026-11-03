# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Validation and normalization for tabular prediction input values."""

import math
from collections.abc import Mapping

from mhctools.pred import value_unit


def finite_prediction_value(value, label: str = "Prediction value"):
    """Return a finite float or ``None`` for a missing/non-finite value.

    Prediction frames use NaN for fields that do not apply to a prediction
    kind.  Converting those values to ``None`` keeps serialized prediction
    evidence standards-compliant while rejecting non-numeric input.
    """
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} {value!r} is not numeric") from error
    return result if math.isfinite(result) else None


def physical_prediction_value(kind, value, measurement_context=None):
    """Return a finite value only when its physical unit is known.

    ``mhctools.Prediction.value`` is a linear physical quantity, while
    predictor-native dimensionless outputs belong in ``score``.  Older
    prediction frames often duplicated those scores into a generic ``value``
    column even for kinds such as ``pMHC_presentation``.  Keep values for
    kinds with a canonical unit, or when the producer supplied an explicit
    units-bearing measurement context; otherwise discard the ambiguous
    duplicate rather than inventing a unit.
    """
    result = finite_prediction_value(value)
    if result is None:
        return None
    if isinstance(measurement_context, Mapping):
        unit = measurement_context.get("unit")
    else:
        unit = getattr(measurement_context, "unit", None)
    return result if unit or value_unit(kind) is not None else None


def prediction_row_fields(row, default_kind: str = "") -> dict:
    """Return the evidence fields shared by one topiary prediction row.

    Parameters
    ----------
    row : mapping-like
        A topiary prediction row (``dict`` or ``pandas.Series``) with optional
        ``kind``, ``prediction_method_name``, ``predictor_version``,
        ``allele``, ``affinity``, ``value``, ``score`` and ``percentile_rank``.
    default_kind : str, optional
        Kind used when the row's ``kind`` is missing or empty.

    Returns
    -------
    dict
        Keyword arguments ``kind``, ``predictor_name``, ``predictor_version``,
        ``allele``, ``score``, ``value`` and ``percentile_rank``. ``value``
        prefers the ``affinity`` column over the generic ``value`` column;
        numeric fields are finite floats or ``None``.

    Raises
    ------
    ValueError
        If a numeric field is present but not numeric.
    """
    value = finite_prediction_value(row.get("affinity"))
    if value is None:
        value = finite_prediction_value(row.get("value"))
    return dict(
        kind=str(row.get("kind") or default_kind),
        predictor_name=str(row.get("prediction_method_name") or ""),
        predictor_version=str(row.get("predictor_version") or ""),
        allele=str(row.get("allele") or ""),
        score=finite_prediction_value(row.get("score")),
        value=value,
        percentile_rank=finite_prediction_value(row.get("percentile_rank")),
    )


def prediction_integer(value, label: str) -> int:
    """Return an exactly integral prediction-frame value.

    Floating point representations such as ``9.0`` are accepted because
    pandas commonly widens integer columns, while booleans and fractional
    values are rejected.
    """
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer")
    try:
        result = int(value)
        equivalent = float(value) == result
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer") from error
    if not equivalent:
        raise ValueError(f"{label} must be an integer")
    return result
