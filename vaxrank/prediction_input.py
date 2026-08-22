# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Validation and normalization for tabular prediction input values."""

import math


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
