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

"""One parse of an external report, shared by every downstream consumer.

An external report is read exactly once. That read produces four things that
must agree with each other:

``report_df``  the user-facing per-(peptide, allele) neoepitope table
``epitopes``   grouped, DSL-scored :class:`CandidateEpitope` objects
``records``    every source row paired with the identity derived from it
``rows``       the source rows themselves, for format-specific evidence

The ``records`` list is what makes the rest safe. Construct ranking used to
re-read the report from disk and re-derive each row's
:class:`ExternalPredictionKey` from a *different* column vocabulary than the
loader had used; whenever the two derivations disagreed — a pVACseq file
carrying both ``Index`` and genomic coordinates, a LENS ``pep_context`` with a
stray trailing space — the join silently produced zero rows and the run
emitted no vaccine peptides. Deriving the identity once and carrying it
forward removes that entire class of failure by construction: there is no
second derivation left to disagree.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import pandas as pd

from .external_prediction import ExternalPredictionKey

# Column added to a normalized row when the source file carried real genomic
# coordinates that the normalization step did not preserve. It is deliberately
# separate from the row's *identity*: the two answer different questions and
# conflating them is what makes a report load but rank nothing.
GENOMIC_VARIANT_COLUMN = "vaxrank_genomic_variant"

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ExternalRecord:
    """One source row of an external report and the identity it produced.

    ``row`` is the row the identity was derived from, kept so format-specific
    evidence (read counts, VAFs, alleles, expression) stays reachable without
    a second pass over the file. ``epitope`` is attached later, once the DSL
    has decided which candidates are eligible for ranking.
    """

    key: ExternalPredictionKey
    row: dict
    epitope: object = None

    def with_epitope(self, epitope) -> "ExternalRecord":
        """Return a copy bound to a scored candidate epitope."""
        return ExternalRecord(key=self.key, row=self.row, epitope=epitope)


@dataclass(frozen=True)
class ExternalReport:
    """Everything one read of an external report produced."""

    source_format: str
    path: str
    report_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    epitopes: tuple = ()
    records: tuple = ()
    rows: tuple = ()

    @property
    def loaded(self):
        """Back-compat view matching the historical loader return shape."""
        return self.report_df, list(self.epitopes)

    def records_with_epitopes(self, epitopes):
        """Bind *epitopes* to the records whose identity selected them.

        *epitopes* is the ranking-eligible subset (post-DSL-filter,
        post-``min_epitope_score``). Records whose identity is absent from it
        are dropped: their candidate did not survive scoring, so they carry no
        evidence a construct may be built from.
        """
        by_id = {}
        for epitope in epitopes:
            if not epitope.prediction_id:
                raise ValueError(
                    "External construct ranking requires prediction "
                    "provenance, but the candidate for peptide %r in %s "
                    "carries none. Candidates must come from a %s reader, "
                    "not be constructed directly."
                    % (epitope.sequence, self.path, self.source_format))
            if epitope.prediction_id in by_id:
                # Reachable only if two candidates normalized to one identity
                # while grouping kept them apart — i.e. a reader is using a
                # different vocabulary for the identity than for the group
                # key. Name both so the offending field is findable.
                other = by_id[epitope.prediction_id]
                raise ValueError(
                    "Two %s candidates share prediction identity %s but were "
                    "grouped separately: %r at offset %d of %r, and %r at "
                    "offset %d of %r. The reader is normalizing the identity "
                    "and the grouping key differently."
                    % (self.source_format, epitope.prediction_id,
                       other.sequence, other.offset, other.source_sequence,
                       epitope.sequence, epitope.offset,
                       epitope.source_sequence))
            by_id[epitope.prediction_id] = epitope
        bound = []
        for record in self.records:
            epitope = by_id.get(record.key.identifier)
            if epitope is not None:
                bound.append(record.with_epitope(epitope))
        if epitopes and not bound:
            # Every candidate was eligible and none joined: the identities the
            # reader produced do not match the identities on the candidates.
            # Silence here is what F1 looked like in production.
            logger.warning(
                "None of the %d ranking-eligible %s candidate(s) matched a "
                "row of %s, so this input will produce no vaccine peptides. "
                "This means prediction identities were derived inconsistently "
                "— it is a bug, not a property of the data.",
                len(epitopes), self.source_format, self.path)
        return tuple(bound)
