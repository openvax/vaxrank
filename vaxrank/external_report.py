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

from dataclasses import dataclass, field

import pandas as pd

from .external_prediction import ExternalPredictionKey

# Column added to a normalized row when the source file carried real genomic
# coordinates that the normalization step did not preserve. It is deliberately
# separate from the row's *identity*: the two answer different questions and
# conflating them is what makes a report load but rank nothing.
GENOMIC_VARIANT_COLUMN = "vaxrank_genomic_variant"


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
                    "provenance")
            if epitope.prediction_id in by_id:
                raise ValueError(
                    "Duplicate external prediction identity after loading")
            by_id[epitope.prediction_id] = epitope
        bound = []
        for record in self.records:
            epitope = by_id.get(record.key.identifier)
            if epitope is not None:
                bound.append(record.with_epitope(epitope))
        return tuple(bound)
