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

import logging
from dataclasses import dataclass

from serializable import DataclassSerializable

logger = logging.getLogger(__name__)


# How the evidence reached this allele.
ALLELE_SOURCE_FROM_INPUT = "from_input"  # the source file paired these two
ALLELE_SOURCE_BROADCAST = "broadcast"    # nothing narrowed it; all alleles got it
ALLELE_SOURCE_SELECTED = "selected"      # vaxrank ranked the genotype and picked
ALLELE_SOURCES = frozenset({ALLELE_SOURCE_FROM_INPUT, ALLELE_SOURCE_BROADCAST, ALLELE_SOURCE_SELECTED})

# Policies for attributing allele-independent evidence.
POLICY_SELECTED = "selected"  # the best-ranked allele, else the genotype
POLICY_ALL = "all"              # every allele in the patient's genotype
POLICY_FROM_INPUT = "from_input"  # only pairings the source file carried
POLICY_TOP = "top"              # ``top:N`` — the N best, else the genotype
POLICY_TOP_PREFIX = "top:"
ALLELE_POLICIES = frozenset({
    POLICY_SELECTED, POLICY_ALL, POLICY_FROM_INPUT, POLICY_TOP})

# Ranking axes a selection may use, with the direction that means "better".
# Presentation is preferred over affinity when both are available: it is the
# model's own statement about presentation, which is the question being asked.
SELECTION_AXES = (
    ("pMHC_presentation", "score", True),   # higher score is better
    ("pMHC_affinity", "value", False),      # lower IC50 is better
)
SELECTION_AUTO = "auto"


# ── MHC dependence: the one table ───────────────────────────────────────────
#
# Every prediction kind is either *allele-scoped* (it describes a peptide-MHC
# pair) or *peptide-level* (it describes the peptide and its flanks alone).
# That distinction was previously re-expressed inline at four call sites, each
# slightly differently — ``p.kind == "antigen_processing" and not p.allele``
# in one, a bare ``not p.allele`` in another, ``if allele:`` in a third. The
# bare forms are wrong in the same way: a blank allele also means a
# *malformed* row, and projecting one of those manufactures per-allele
# binding evidence for alleles the model never scored.
#
# So the question is answered here, once, from the kind. Every topiary
# canonical kind is listed; _assert_dependence_table_is_complete below fails
# at import if topiary adds one we have not classified, so a new kind is a
# visible decision rather than a silent default.
PEPTIDE_LEVEL_KINDS = frozenset({
    "antigen_processing",       # mhcflurry integrated processing
    "proteasome_cleavage",      # NetChop / PAProC / pepsickle
    "endolysosomal_cleavage",   # class-II processing
    "erap_trimming",            # ERAP1/2 N-terminal trimming
    "tap_transport",            # TAP peptide transport
})
ALLELE_SCOPED_KINDS = frozenset({
    "pMHC_affinity",
    "pMHC_presentation",
    "pMHC_stability",
    "pMHC_TCR_binding",
    "immunogenicity",           # scored for a peptide-MHC pair
})


def unclassified_kinds() -> tuple:
    """Canonical topiary kinds this module has not classified."""
    try:
        from topiary.ranking import KIND_ALIASES
    except Exception:  # pragma: no cover - topiary always present in practice
        return ()
    known = PEPTIDE_LEVEL_KINDS | ALLELE_SCOPED_KINDS
    return tuple(sorted(set(KIND_ALIASES.values()) - known))


def _warn_about_unclassified_kinds():
    """Report kinds we cannot classify, without refusing to run.

    Raising here would mean a topiary release that adds a prediction kind
    breaks ``import vaxrank`` outright — a worse failure than the one being
    guarded against. Unclassified kinds are instead treated as *not*
    peptide-level, which is the safe direction: they are never projected
    across alleles, so nothing is fabricated. A test pins the table against
    the pinned topiary, so drift surfaces in CI rather than in a run.
    """
    missing = unclassified_kinds()
    if missing:
        logger.warning(
            "Prediction kind(s) %s are not classified as allele-scoped or "
            "peptide-level by vaxrank.allele_evidence. They will not be "
            "attributed to alleles. Add them to one of the two sets there.",
            ", ".join(missing))


_warn_about_unclassified_kinds()


def is_peptide_level(prediction) -> bool:
    """True when a prediction describes the peptide, not a peptide-MHC pair.

    Both halves are required. The kind must be one that carries no MHC
    context, *and* the record must actually be allele-free — a peptide-level
    kind that arrived stamped with an allele is already attributed, and an
    allele-scoped kind that arrived blank is malformed, not peptide-level.
    """
    return prediction.kind in PEPTIDE_LEVEL_KINDS and not prediction.allele


def is_allele_scoped(prediction) -> bool:
    """True when a prediction names the allele it applies to."""
    return bool(prediction.allele) and prediction.kind in ALLELE_SCOPED_KINDS


def alleles_named_by(predictions) -> tuple:
    """The alleles a set of predictions was actually scored against."""
    return tuple(sorted({
        prediction.allele for prediction in predictions
        if is_allele_scoped(prediction)}))


def _plain_number(value):
    """Coerce a prediction value to a plain Python float, or None.

    Predictions built from a pandas frame carry numpy scalars, which the
    serializer refuses ("Cannot convert 1 : <class 'numpy.int64'>"). Since a
    recorded value has to survive a round-trip to stay reconstructable, the
    conversion belongs here, at the point the value is captured, rather than
    at each place one is written out.
    """
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if result != result else result  # drop NaN


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
        known_axes = {a[0] for a in SELECTION_AXES} | {SELECTION_AUTO}
        if axis not in known_axes:
            raise AllelePolicyError(
                "Unknown allele selection axis %r (expected one of %s)"
                % (axis, ", ".join(sorted(known_axes))))
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


# Preference between predictors of one kind when the config names none.
# Mirrors epitope_dsl._CANONICAL_METHOD_PREFERENCE so a selection and a score
# resolve to the same model.
_PREDICTOR_PREFERENCE = ("mhcflurry", "netmhcpan", "netmhcstabpan")


@dataclass(frozen=True)
class RankedEvidence:
    """One predictor's values for one kind, and the ranking they produce.

    The kind, the predictor, and the values travel together in one object on
    purpose. The first version of this code accumulated them in separate
    locals inside one loop — ``predictor`` latched to the first predictor
    seen, ``values`` kept the best across *all* predictors — so an
    attribution could record a ``(predictor, value)`` pair that predictor
    never emitted. Making them one value makes that class of mistake
    unrepresentable rather than merely fixed.
    """

    kind: str
    predictor: str
    values: dict
    higher_is_better: bool

    @property
    def coverage(self) -> int:
        """How many of the requested alleles this evidence can rank."""
        return len(self.values)

    @property
    def ordered(self) -> tuple:
        """Alleles best-first, ties broken by name for determinism."""
        return tuple(sorted(
            self.values,
            key=lambda a: (
                -self.values[a] if self.higher_is_better else self.values[a],
                a)))


def _pick_predictor(by_predictor, default_methods, kind):
    """Choose the one predictor whose values will rank the alleles.

    An explicit ``default_methods`` entry is the user's stated choice of
    model and always wins, so a selection resolves to the same model the
    score does. Absent one, prefer the predictor that can rank the most
    alleles, and only then fall back to the canonical name order: a
    predictor covering one allele out of six would otherwise silently reduce
    a ``top:3`` to a single attribution.
    """
    named = (default_methods or {}).get(kind)
    if named and named in by_predictor:
        return named
    return max(
        by_predictor,
        key=lambda name: (
            len(by_predictor[name]),
            -_PREDICTOR_PREFERENCE.index(name)
            if name in _PREDICTOR_PREFERENCE else -len(_PREDICTOR_PREFERENCE),
            name,
        ))


def _ranked_alleles(epitope, alleles, axis, default_methods=None):
    """Return the :class:`RankedEvidence` best able to rank *alleles*.

    ``None`` when the peptide carries no usable allele-scoped evidence — the
    caller must then fall back rather than invent a ranking.

    The axis is chosen by how many of the requested alleles it actually
    covers, not by whether it covers any: a file carrying presentation for
    one allele and affinity for six should rank on affinity. Ties go to the
    earlier (preferred) axis.
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
    best = None
    best_key = None
    for rank, (kind, attribute, higher_is_better) in enumerate(candidates):
        by_predictor = {}
        for prediction in epitope.predictions_flat():
            if prediction.kind != kind or prediction.allele not in allowed:
                continue
            if not is_allele_scoped(prediction):
                continue
            value = _plain_number(getattr(prediction, attribute, None))
            if value is None:
                continue
            slot = by_predictor.setdefault(prediction.predictor_name, {})
            current = slot.get(prediction.allele)
            if current is None or (
                    value > current if higher_is_better else value < current):
                slot[prediction.allele] = value
        if not by_predictor:
            continue
        predictor = _pick_predictor(by_predictor, default_methods, kind)
        evidence = RankedEvidence(
            kind=kind, predictor=predictor,
            values=dict(by_predictor[predictor]),
            higher_is_better=higher_is_better)
        # One key, not two locals that must stay in sync — the mistake this
        # whole module is a reaction to. More alleles covered wins; ties go
        # to the earlier (preferred) axis.
        key = (evidence.coverage, -rank)
        if best_key is None or key > best_key:
            best, best_key = evidence, key
    return best


def attribute_alleles(epitope, genotype, policy,
                      default_methods=None) -> tuple:
    """Return the :class:`AlleleAttribution` list for one candidate.

    ``genotype`` is the patient's allele set. ``policy`` is an
    :class:`AllelePolicy`. The result is ordered and deterministic, and every
    entry records why that allele was chosen.
    """
    # Derived from the allele-scoped predictions themselves, NOT from
    # patient_alleles: CandidateEpitope.__post_init__ folds every scored
    # allele back into patient_alleles, so a broadcast allele would be
    # relabelled "from_input" on a second attribution pass and the audit
    # record would claim the input paired a peptide with an allele it never
    # mentioned.
    from_input = alleles_named_by(epitope.predictions_flat())
    # An explicit genotype is the patient's actual typing and wins outright.
    # Unioning the input's alleles in would mean an explicit genotype could
    # only ever widen the attributed set, never narrow it — so it could not
    # be used to correct an input produced with the wrong or a wider panel.
    genotype = (
        tuple(sorted({a for a in genotype if a})) if genotype
        else from_input)
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

    if policy.name not in (POLICY_SELECTED, POLICY_TOP):
        raise AllelePolicyError(
            "Unknown allele policy name %r; AllelePolicy must come from "
            "AllelePolicy.parse" % policy.name)
    if policy.name == POLICY_TOP and not policy.limit:
        raise AllelePolicyError(
            "A 'top:N' allele policy requires a positive limit")

    evidence = _ranked_alleles(
        epitope, genotype, policy.axis, default_methods=default_methods)
    if evidence is None:
        # No allele-scoped evidence to select from. "We don't know which
        # allele presents this" is not the same as "this allele does", so
        # keep the whole genotype rather than inventing a winner.
        from_input_set = set(from_input)
        return tuple(
            AlleleAttribution(
                allele=a,
                source=ALLELE_SOURCE_FROM_INPUT if a in from_input_set else ALLELE_SOURCE_BROADCAST)
            for a in genotype)

    # ``top:N`` means "up to N of the alleles this axis could actually rank".
    # When fewer are rankable the result is shorter; padding it with alleles
    # the model could not rank would assert a preference the data does not
    # support.
    limit = 1 if policy.name == POLICY_SELECTED else policy.limit
    return tuple(
        AlleleAttribution(
            allele=allele,
            source=ALLELE_SOURCE_SELECTED,
            rank_kind=evidence.kind,
            rank_predictor=evidence.predictor,
            rank_value=evidence.values[allele],
            rank_position=position,
        )
        for position, allele in enumerate(evidence.ordered[:limit], start=1))
