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

"""Which MHC alleles allele-independent evidence is attributed to, and why.

Some predictions describe a peptide, not a peptide-MHC pair: antigen
processing and proteasomal cleavage depend on the peptide and its flanks
with no MHC involved. Vaxrank's scoring substrate is per-allele
(``CandidateEpitope.per_allele_scores``), so that peptide-level evidence has
to be attributed to one or more alleles before it can score anything.

That attribution is a *choice*, and it used to be made implicitly: the
evidence was projected onto whichever alleles the input file happened to
score that peptide against. Row incidence is neither the patient's genotype
nor a claim about which allele presents the peptide — on real LENS files
those sets are one or two alleles out of a five- or six-allele genotype, so
the implicit answer was an artifact of how the upstream tool chose to write
its rows.

This module makes the choice explicit and, crucially, **recordable**. Every
attributed allele carries an :class:`AlleleAttribution` saying where it came
from: it came from the input file, it was broadcast across the genotype, or
vaxrank selected it by ranking — and for a selection, on
which axis, by which predictor, at what value and what rank. A reader of a
finished run can therefore reconstruct exactly why a given allele was
credited with a given processing score, rather than having to re-derive it.

Note that this cannot be delegated to topiary's ``alleles=`` context option:
that declares one allele set for a whole frame, while every policy except
:data:`POLICY_ALL` needs a different set per peptide.
"""

from __future__ import annotations

from dataclasses import dataclass

from serializable import DataclassSerializable


# How the evidence reached this allele.
ALLELE_SOURCE_FROM_INPUT = "from_input"  # the source file paired these two
ALLELE_SOURCE_BROADCAST = "broadcast"    # nothing narrowed it; all alleles got it
ALLELE_SOURCE_SELECTED = "selected"      # vaxrank ranked the genotype and picked
ALLELE_SOURCES = frozenset({ALLELE_SOURCE_FROM_INPUT, ALLELE_SOURCE_BROADCAST, ALLELE_SOURCE_SELECTED})

# Policies for attributing allele-independent evidence.
POLICY_SELECTED = "selected"  # the best-ranked allele, else the genotype
POLICY_ALL = "all"              # every allele in the patient's genotype
POLICY_FROM_INPUT = "from_input"  # only pairings the source file carried
POLICY_TOP_PREFIX = "top:"      # ``top:N`` — the N best, else the genotype

# Ranking axes a selection may use, with the direction that means "better".
# Presentation is preferred over affinity when both are available: it is the
# model's own statement about presentation, which is the question being asked.
SELECTION_AXES = (
    ("pMHC_presentation", "score", True),   # higher score is better
    ("pMHC_affinity", "value", False),      # lower IC50 is better
)
SELECTION_AUTO = "auto"


class AllelePolicyError(ValueError):
    """Raised for an allele-attribution policy string that cannot be honored."""


@dataclass(frozen=True)
class AlleleAttribution(DataclassSerializable):
    """Why one allele was credited with a peptide's allele-free evidence.

    ``source`` is one of :data:`ALLELE_SOURCE_FROM_INPUT`, :data:`ALLELE_SOURCE_BROADCAST`, or
    :data:`ALLELE_SOURCE_SELECTED`. The remaining fields are populated only for a
    selection and are what makes it reproducible: the axis that ranked the
    alleles, the predictor whose value was read, the value itself, and the
    resulting 1-based position.
    """

    allele: str
    source: str
    rank_kind: str = ""
    rank_predictor: str = ""
    rank_value: object = None
    rank_position: object = None

    def __post_init__(self):
        if not self.allele:
            raise ValueError("An allele attribution requires an allele")
        if self.source not in ALLELE_SOURCES:
            raise ValueError(
                "Unknown allele attribution source %r (expected one of %s)"
                % (self.source, ", ".join(sorted(ALLELE_SOURCES))))
        if self.source == ALLELE_SOURCE_SELECTED and not self.rank_kind:
            raise ValueError(
                "A selected allele must record the axis that ranked it, or "
                "the selection cannot be reconstructed")
        if self.source != ALLELE_SOURCE_SELECTED and self.rank_kind:
            raise ValueError(
                "Only a selected allele carries ranking provenance")

    def describe(self) -> str:
        """One human-readable line, for reports and logs."""
        if self.source == ALLELE_SOURCE_FROM_INPUT:
            return "%s (from the input file)" % self.allele
        if self.source == ALLELE_SOURCE_BROADCAST:
            return "%s (broadcast: no allele-specific evidence)" % self.allele
        return "%s (selected: ranked #%s on %s%s = %s)" % (
            self.allele,
            self.rank_position if self.rank_position is not None else "?",
            self.rank_kind,
            " [%s]" % self.rank_predictor if self.rank_predictor else "",
            self.rank_value,
        )


@dataclass(frozen=True)
class AllelePolicy:
    """A parsed allele-attribution policy."""

    name: str = POLICY_SELECTED
    limit: object = None            # for ``top:N``
    axis: str = SELECTION_AUTO     # kind used to rank, or "auto"

    @classmethod
    def parse(cls, value, axis=SELECTION_AUTO) -> "AllelePolicy":
        """Parse a config string into a policy.

        Accepts ``selected``, ``all``, ``from_input``, and ``top:N``.
        """
        text = (value or POLICY_SELECTED).strip().lower()
        if text in (POLICY_SELECTED, POLICY_ALL, POLICY_FROM_INPUT):
            return cls(name=text, axis=axis)
        if text.startswith(POLICY_TOP_PREFIX):
            suffix = text[len(POLICY_TOP_PREFIX):]
            try:
                limit = int(suffix)
            except ValueError:
                raise AllelePolicyError(
                    "Allele policy %r must be 'top:N' with an integer N"
                    % value) from None
            if limit < 1:
                raise AllelePolicyError(
                    "Allele policy %r must keep at least one allele" % value)
            return cls(name=POLICY_TOP_PREFIX.rstrip(":"), limit=limit,
                       axis=axis)
        raise AllelePolicyError(
            "Unknown allele policy %r (expected 'selected', 'all', "
            "'from_input', or 'top:N')" % value)


def _ranked_alleles(epitope, alleles, axis):
    """Rank *alleles* by the peptide's own allele-scoped evidence.

    Returns ``(ordered, kind, predictor, values)``, or ``(None, ...)`` when
    the peptide carries no usable allele-scoped evidence — in which case the
    caller must fall back rather than invent a ranking.
    """
    candidates = (
        SELECTION_AXES if axis == SELECTION_AUTO
        else tuple(a for a in SELECTION_AXES if a[0] == axis))
    if axis != SELECTION_AUTO and not candidates:
        raise AllelePolicyError(
            "Unknown allele selection axis %r (expected one of %s or %r)"
            % (axis, ", ".join(a[0] for a in SELECTION_AXES),
               SELECTION_AUTO))

    allowed = set(alleles)
    for kind, attribute, higher_is_better in candidates:
        by_allele = {}
        predictor = ""
        for prediction in epitope.predictions_flat():
            if prediction.kind != kind or prediction.allele not in allowed:
                continue
            value = getattr(prediction, attribute, None)
            if value is None:
                continue
            predictor = predictor or prediction.predictor_name
            best = by_allele.get(prediction.allele)
            if best is None or (
                    value > best if higher_is_better else value < best):
                by_allele[prediction.allele] = value
        if by_allele:
            ordered = sorted(
                by_allele,
                key=lambda a: (-by_allele[a] if higher_is_better
                               else by_allele[a], a))
            return ordered, kind, predictor, by_allele
    return None, "", "", {}


def attribute_alleles(epitope, genotype, policy) -> tuple:
    """Return the :class:`AlleleAttribution` list for one candidate.

    ``genotype`` is the patient's allele set. ``policy`` is an
    :class:`AllelePolicy`. The result is ordered and deterministic, and every
    entry records why that allele was chosen.
    """
    from_input = tuple(a for a in epitope.patient_alleles if a)
    genotype = tuple(sorted({a for a in (genotype or ()) if a} | set(from_input)))
    if not genotype:
        return ()

    if policy.name == POLICY_FROM_INPUT:
        return tuple(
            AlleleAttribution(allele=a, source=ALLELE_SOURCE_FROM_INPUT)
            for a in sorted(from_input))

    if policy.name == POLICY_ALL:
        from_input_set = set(from_input)
        return tuple(
            AlleleAttribution(
                allele=a,
                source=ALLELE_SOURCE_FROM_INPUT if a in from_input_set else ALLELE_SOURCE_BROADCAST)
            for a in genotype)

    # selected / top:N
    ordered, kind, predictor, values = _ranked_alleles(
        epitope, genotype, policy.axis)
    if ordered is None:
        # No allele-scoped evidence to select from. "We don't know which
        # allele presents this" is not the same as "this allele does", so
        # keep the whole genotype rather than inventing a winner.
        from_input_set = set(from_input)
        return tuple(
            AlleleAttribution(
                allele=a,
                source=ALLELE_SOURCE_FROM_INPUT if a in from_input_set else ALLELE_SOURCE_BROADCAST)
            for a in genotype)

    limit = 1 if policy.name == POLICY_SELECTED else policy.limit
    return tuple(
        AlleleAttribution(
            allele=allele,
            source=ALLELE_SOURCE_SELECTED,
            rank_kind=kind,
            rank_predictor=predictor,
            rank_value=values[allele],
            rank_position=position,
        )
        for position, allele in enumerate(ordered[:limit], start=1))
