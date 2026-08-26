# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The single normalization vocabulary for tabular input cells.

Every reader of an external report (LENS, pVACseq) and of vaxrank's own
native epitope files funnels cell values through these functions. There is
deliberately *one* answer to "is this cell missing", "what is this cell as
text", and "what is this cell as a number", because a prediction identity
built with one answer must join against an identity built with another.

The historical failure mode this replaces: ``epitope_io._safe_str`` did not
strip whitespace while ``external_prediction.external_text`` did, so a single
trailing space in a LENS ``pep_context`` produced two different identities for
the same row and aborted the run with an opaque duplicate-identity error.
"""

from __future__ import annotations

import pandas as pd


# Cell text that means "no value" regardless of the producing tool. Compared
# case-insensitively after stripping.
MISSING_TEXT = frozenset({"", "na", "n/a", "nan", "none", "null", "."})

_TRUE_TEXT = frozenset({"true", "t", "yes", "y", "1"})
_FALSE_TEXT = frozenset({"false", "f", "no", "n", "0"})


def missing(value) -> bool:
    """True when a cell carries no value.

    Covers ``None``, pandas/numpy NA sentinels, and the textual spellings in
    :data:`MISSING_TEXT` that TSV producers use interchangeably.
    """
    if value is None:
        return True
    try:
        if bool(pd.isna(value)):
            return True
    except (TypeError, ValueError):
        # Array-likes and exotic objects are not scalar-missing.
        pass
    if isinstance(value, str):
        return value.strip().lower() in MISSING_TEXT
    return False


def text(value) -> str:
    """Return a cell as stripped text, or ``""`` when it is missing."""
    if missing(value):
        return ""
    return str(value).strip()


def values(*cells) -> tuple[str, ...]:
    """Return the sorted unique comma/semicolon-delimited values of *cells*.

    External reports pack multi-valued annotations (gene IDs, transcript IDs)
    into one cell with either separator, and different rows list them in
    different orders. Sorting and de-duplicating here is what lets those rows
    resolve to one identity.
    """
    result = set()
    for cell in cells:
        cell_text = text(cell)
        if not cell_text:
            continue
        for comma_part in cell_text.split(","):
            for part in comma_part.split(";"):
                part = part.strip()
                if part:
                    result.add(part)
    return tuple(sorted(result))


def number(value, default=None):
    """Return a cell as ``float``, or *default* when missing/unparseable."""
    if missing(value):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def integer(value, default=0) -> int:
    """Return a cell as a rounded ``int``, or *default* when unusable.

    Read counts and depths arrive as floats, as float-formatted strings, or
    not at all. A missing count is *default* (normally ``0``, an honest "no
    signal") rather than a fabricated stand-in.
    """
    parsed = number(value)
    if parsed is None:
        return default
    try:
        return int(round(parsed))
    except (TypeError, ValueError, OverflowError):
        return default


def boolean(value, default=False) -> bool:
    """Return a cell as ``bool``, honoring the usual TSV spellings."""
    if missing(value):
        return default
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in _TRUE_TEXT:
            return True
        if lowered in _FALSE_TEXT:
            return False
        return default
    return bool(value)


def first(row, *names):
    """Return the first present, non-missing value among *names* in *row*.

    The ordering of *names* is a documented precedence, not an accident:
    when one logical field is spelled differently across report flavors,
    every reader must consult the spellings in the same order or the two
    readers disagree about what the row says.
    """
    for name in names:
        value = row.get(name)
        if not missing(value):
            return value
    return None


def first_text(row, *names) -> str:
    """:func:`first` coerced through :func:`text`."""
    return text(first(row, *names))


def first_number(row, *names, default=None):
    """:func:`first` coerced through :func:`number`."""
    return number(first(row, *names), default=default)
