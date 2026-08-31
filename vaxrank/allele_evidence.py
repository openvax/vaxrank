# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Which alleles are credited with a peptide's allele-free evidence.

Some predictions describe a peptide-MHC pair — affinity, presentation,
stability. Others describe the peptide alone: proteasome cleavage, TAP
transport, ERAP trimming, mhcflurry's integrated processing score. The
second group has no allele, and vaxrank's per-allele scores do, so
something has to decide which alleles that evidence counts for.

Broadcasting it to the whole genotype is one answer and the historical
one. It is also a claim: that every allele the patient carries presents
this peptide equally well, which the allele-scoped evidence usually
contradicts. This module makes the choice explicit and records it, so a
ranking can be read back and explained.

Two things are deliberately separate here, and confusing them is the
mistake this module exists to avoid:

*Broadcasting* is data shape. topiary's ``peptide_view()`` spreads an
allele-free value across the allele groups a peptide already has, and
:func:`vaxrank.epitope_dsl.genotype_lookup` declares which groups exist.
That is not policy — it is what makes an allele-free value readable at
all in a grouping keyed by allele.

*Attribution* is policy. It decides which of those alleles should be
credited, and on what evidence. That is what lives here.

The unit is a **topiary frame**, not a ``CandidateEpitope``. Both input
paths share a frame — the native predictor pipeline and the external
report readers — and only one of them ever builds CandidateEpitopes at
the point scoring happens. Attribution written against CandidateEpitope
applied to LENS runs and silently did nothing on the path most runs take
(openvax/vaxrank#349).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from serializable import DataclassSerializable

logger = logging.getLogger(__name__)


# ── Vocabulary ──────────────────────────────────────────────────────────────

# Why an allele ended up credited. Named for the mechanism that produced
# the attribution, not for the role the allele plays.
ALLELE_SOURCE_FROM_INPUT = "from_input"  # the source file paired these two
ALLELE_SOURCE_BROADCAST = "broadcast"    # nothing narrowed it; all alleles got it
ALLELE_SOURCE_SELECTED = "selected"      # vaxrank ranked the genotype and picked
ALLELE_SOURCES = frozenset({
    ALLELE_SOURCE_FROM_INPUT, ALLELE_SOURCE_BROADCAST, ALLELE_SOURCE_SELECTED})

# Policies for attributing allele-free evidence.
POLICY_ALL = "all"                # every allele in the patient's genotype
POLICY_SELECTED = "selected"      # the best-ranked allele, else the genotype
POLICY_FROM_INPUT = "from_input"  # only pairings the source file carried
POLICY_TOP = "top"                # ``top:N`` — the N best, else the genotype
POLICY_TOP_PREFIX = "top:"
ALLELE_POLICIES = frozenset({
    POLICY_ALL, POLICY_SELECTED, POLICY_FROM_INPUT, POLICY_TOP})

# Ranking axes a selection may use, with the direction that means "better"
# and the frame column carrying the value. Presentation is preferred over
# affinity when both are available: it is the model's own statement about
# presentation, which is the question being asked.
SELECTION_AXES = (
    ("pMHC_presentation", "score", True),   # higher score is better
    ("pMHC_affinity", "value", False),      # lower IC50 is better
)
SELECTION_AUTO = "auto"

# Preference between predictors of one kind when the config names none.
# Ordered to match topiary's CANONICAL_METHOD_PREFERENCE so a selection and
# a score resolve to the same model.
_PREDICTOR_PREFERENCE = ("mhcflurry", "netmhcpan", "netmhcstabpan")


class AllelePolicyError(ValueError):
    """A configured allele policy could not be understood."""


def is_peptide_level_kind(kind) -> bool:
    """True when a kind describes the peptide rather than a peptide-MHC pair.

    Asks ``topiary.KIND_MHC_DEPENDENCE``, which is the one table that knows
    (openvax/vaxrank#357). vaxrank kept its own copy of this partition while
    that table was private, and the copy named one kind where topiary names
    five.

    An unrecognized kind is *not* peptide-level. That is the safe direction:
    it is never attributed to an allele, so a topiary release adding a kind
    cannot make vaxrank credit an allele with evidence it has not classified.
    """
    from topiary import KIND_MHC_DEPENDENCE

    return KIND_MHC_DEPENDENCE.get(kind) == "none"


def is_allele_scoped_kind(kind) -> bool:
    """True when a kind describes a peptide-MHC pair."""
    from topiary import KIND_MHC_DEPENDENCE

    return KIND_MHC_DEPENDENCE.get(kind) == "single_allele"


@dataclass(frozen=True)
class AllelePolicy:
    """A parsed allele-attribution policy."""

    name: str = POLICY_ALL
    limit: object = None           # for ``top:N``
    axis: str = SELECTION_AUTO     # kind used to rank, or "auto"

    @classmethod
    def parse(cls, value, axis=SELECTION_AUTO) -> "AllelePolicy":
        """Parse a config string into a policy.

        Accepts ``all``, ``selected``, ``from_input``, and ``top:N``.
        """
        known_axes = {a[0] for a in SELECTION_AXES} | {SELECTION_AUTO}
        if axis not in known_axes:
            raise AllelePolicyError(
                "Unknown allele selection axis %r (expected one of %s)"
                % (axis, ", ".join(sorted(known_axes))))
        text = (value or POLICY_ALL).strip().lower()
        if text in (POLICY_ALL, POLICY_SELECTED, POLICY_FROM_INPUT):
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
            return cls(name=POLICY_TOP, limit=limit, axis=axis)
        raise AllelePolicyError(
            "Unknown allele policy %r (expected 'all', 'selected', "
            "'from_input', or 'top:N')" % value)

    @property
    def narrows(self) -> bool:
        """True when this policy can credit fewer alleles than the genotype.

        ``all`` cannot, so what is recorded and what is scored agree.
        Anything else can differ from the scores until openvax/topiary#232
        lands, which is what :func:`warn_if_scoring_ignores_policy` says.
        """
        return self.name != POLICY_ALL


@dataclass(frozen=True)
class AlleleAttribution(DataclassSerializable):
    """Why one allele was credited with a peptide's allele-free evidence.

    ``source`` is one of :data:`ALLELE_SOURCE_FROM_INPUT`,
    :data:`ALLELE_SOURCE_BROADCAST` or :data:`ALLELE_SOURCE_SELECTED`. The
    remaining fields are populated only for a selection and are what makes it
    reproducible: the axis that ranked the alleles, the predictor whose value
    was read, the value itself, and the resulting 1-based position.
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

    def reason(self) -> str:
        """Why this allele was credited, without naming it.

        Separate from :meth:`describe` because the report already has an
        Allele column, and the alternative was slicing the allele back off
        the front of ``describe()`` — string surgery that would break
        silently the first time the sentence changed.
        """
        if self.source == ALLELE_SOURCE_FROM_INPUT:
            return "from the input file"
        if self.source == ALLELE_SOURCE_BROADCAST:
            return "broadcast: no allele-specific evidence"
        return "selected: ranked #%s on %s%s = %s" % (
            self.rank_position if self.rank_position is not None else "?",
            self.rank_kind,
            " [%s]" % self.rank_predictor if self.rank_predictor else "",
            self.rank_value,
        )

    def describe(self) -> str:
        """One human-readable line, for reports and logs."""
        return "%s (%s)" % (self.allele, self.reason())


@dataclass(frozen=True)
class RankedEvidence:
    """One predictor's values for one kind, and the ranking they produce.

    The kind, the predictor and the values travel together in one object on
    purpose. An earlier version accumulated them in separate locals inside
    one loop — ``predictor`` latched to the first predictor seen while
    ``values`` kept the best across *all* predictors — so an attribution
    could record a ``(predictor, value)`` pair that predictor never emitted.
    Making them one value makes that class of mistake unrepresentable rather
    than merely fixed.
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


def _plain_number(value):
    """Coerce a frame cell to a plain float, or ``None`` if it is not one.

    numpy scalars survive arithmetic and comparison unchanged, so they reach
    the attribution record and then fail serialization with "Cannot convert
    1 : numpy.int64". Coercing here keeps the record made of builtins.
    """
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number:  # NaN
        return None
    return number


def _pick_predictor(by_predictor, default_methods, kind):
    """Choose the one predictor whose values will rank the alleles.

    An explicit ``default_methods`` entry is the user's stated choice of
    model and always wins, so a selection resolves to the same model the
    score does. Absent one, prefer the predictor that can rank the most
    alleles, and only then fall back to the canonical name order: a predictor
    covering one allele out of six would otherwise silently reduce a
    ``top:3`` to a single attribution.
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


def _ranked_alleles(rows, alleles, axis, default_methods=None):
    """Return the :class:`RankedEvidence` best able to rank *alleles*.

    ``rows`` is one peptide's slice of a topiary frame, as a list of dicts.
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
            % (axis, ", ".join(a[0] for a in SELECTION_AXES), SELECTION_AUTO))

    allowed = set(alleles)
    best = None
    best_key = None
    for rank, (kind, column, higher_is_better) in enumerate(candidates):
        by_predictor = {}
        for row in rows:
            allele = row.get("allele") or ""
            if row.get("kind") != kind or allele not in allowed:
                continue
            if not is_allele_scoped_kind(row.get("kind")):
                continue
            value = _plain_number(row.get(column))
            if value is None:
                continue
            slot = by_predictor.setdefault(
                row.get("prediction_method_name") or "", {})
            current = slot.get(allele)
            if current is None or (
                    value > current if higher_is_better else value < current):
                slot[allele] = value
        if not by_predictor:
            continue
        predictor = _pick_predictor(by_predictor, default_methods, kind)
        evidence = RankedEvidence(
            kind=kind, predictor=predictor,
            values=dict(by_predictor[predictor]),
            higher_is_better=higher_is_better)
        # One key, not two locals that must stay in sync. More alleles
        # covered wins; ties go to the earlier (preferred) axis.
        key = (evidence.coverage, -rank)
        if best_key is None or key > best_key:
            best, best_key = evidence, key
    return best


def attribute_peptide(rows, genotype, policy, default_methods=None) -> tuple:
    """Return the :class:`AlleleAttribution` list for one peptide's rows.

    ``rows`` is that peptide's slice of a topiary frame. ``genotype`` is the
    patient's allele set, or empty to use whatever the input paired.
    """
    # Derived from the allele-scoped rows themselves, NOT from the genotype:
    # a broadcast allele must not be relabelled "from_input" on a second
    # pass, or the audit record would claim the input paired a peptide with
    # an allele it never mentioned.
    from_input = tuple(sorted({
        row.get("allele") or ""
        for row in rows
        if (row.get("allele") or "") and is_allele_scoped_kind(row.get("kind"))
    }))
    # An explicit genotype is the patient's actual typing and wins outright.
    # Unioning the input's alleles in would mean an explicit genotype could
    # only ever widen the attributed set, never narrow it — so it could not
    # be used to correct an input produced with the wrong or a wider panel.
    genotype = (
        tuple(sorted({a for a in genotype if a})) if genotype else from_input)
    if not genotype:
        return ()

    if policy.name == POLICY_FROM_INPUT:
        return tuple(
            AlleleAttribution(allele=a, source=ALLELE_SOURCE_FROM_INPUT)
            for a in from_input)

    from_input_set = set(from_input)

    def unranked():
        return tuple(
            AlleleAttribution(
                allele=a,
                source=(ALLELE_SOURCE_FROM_INPUT if a in from_input_set
                        else ALLELE_SOURCE_BROADCAST))
            for a in genotype)

    if policy.name == POLICY_ALL:
        return unranked()

    if policy.name not in (POLICY_SELECTED, POLICY_TOP):
        raise AllelePolicyError(
            "Unknown allele policy name %r; AllelePolicy must come from "
            "AllelePolicy.parse" % policy.name)
    if policy.name == POLICY_TOP and not policy.limit:
        raise AllelePolicyError(
            "A 'top:N' allele policy requires a positive limit")

    evidence = _ranked_alleles(
        rows, genotype, policy.axis, default_methods=default_methods)
    if evidence is None:
        # No allele-scoped evidence to select from. "We don't know which
        # allele presents this" is not the same as "this allele does", so
        # keep the whole genotype rather than inventing a winner.
        return unranked()

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


def _identity_columns(group_columns):
    """The group columns that identify a peptide, i.e. all but the allele."""
    return [column for column in group_columns if column != "allele"]


def attribute_frame(df, policy, *, genotype_for=None, default_methods=None,
                    group_columns=None) -> dict:
    """Attribute every peptide in a topiary frame.

    Returns ``{identity_tuple: (AlleleAttribution, ...)}`` keyed by the group
    columns other than ``allele`` — the same identity both input paths group
    on, so the result joins back to either.

    ``genotype_for`` is called with an identity tuple and returns that
    peptide's genotype. Per peptide rather than per frame because the two
    paths differ: the native pipeline predicts every peptide against the
    patient's full typing, while LENS reports only the pairs that passed its
    own threshold, so candidates in one file arrive with different allele
    sets.

    Peptides with no allele-free evidence are skipped. Attribution is a
    statement about which alleles that evidence counts for, and a peptide
    that has none has nothing to attribute.
    """
    if df is None or len(df) == 0:
        return {}
    # Cheap vectorized pre-check before materializing records. Most frames
    # carry no allele-free evidence at all, and turning a large frame into
    # dicts to discover that costs real time on every scoring pass.
    if "kind" not in df.columns:
        return {}
    allele_free = df["allele"].map(lambda a: not (a or "")) if (
        "allele" in df.columns) else True
    peptide_level = df["kind"].map(is_peptide_level_kind)
    if not bool((peptide_level & allele_free).any()):
        return {}
    group_columns = list(group_columns or ())
    identity = _identity_columns(group_columns)
    if not identity:
        raise ValueError(
            "Allele attribution needs peptide-identity group columns; got %r"
            % (group_columns,))

    records = df.to_dict("records")
    rows_by_identity: dict = {}
    has_allele_free: set = set()
    for row in records:
        key = tuple(row.get(column) for column in identity)
        rows_by_identity.setdefault(key, []).append(row)
        if is_peptide_level_kind(row.get("kind")) and not (
                row.get("allele") or ""):
            has_allele_free.add(key)

    attributions = {}
    for key in has_allele_free:
        genotype = tuple(genotype_for(key)) if genotype_for else ()
        attributed = attribute_peptide(
            rows_by_identity[key], genotype, policy,
            default_methods=default_methods)
        if attributed:
            attributions[key] = attributed
    return attributions


_WARNED_ABOUT_SCORING = set()


def warn_if_scoring_ignores_policy(policy, attributions):
    """Say plainly that a narrowing policy does not yet change scores.

    Attribution records which alleles are credited with a peptide's
    allele-free evidence. It does not currently change the per-allele
    scores, because topiary projects a peptide-level value across every
    allele group a peptide has, keyed on the *kind* — a row that names an
    allele is projected onto the others anyway (openvax/topiary#232). So
    "credit this to the best allele" is recorded and reported, and every
    allele still receives the value when the score expression reads it.

    Writing the row onto the attributed alleles was tried and does nothing:
    the frame changes and the scores do not, which is worse than the gap
    itself. Better to state the limitation than to ship a frame that says
    something the evaluation discards.

    The default policy is ``all``, where recorded and scored agree and
    there is nothing to warn about.
    """
    if not policy.narrows or not attributions:
        return
    key = (policy.name, policy.limit)
    if key in _WARNED_ABOUT_SCORING:
        return
    _WARNED_ABOUT_SCORING.add(key)
    logger.warning(
        "allele_free_evidence=%s records which allele is credited with "
        "peptide-level evidence, and is reported per candidate, but does "
        "not yet change per-allele scores: topiary projects a peptide-level "
        "value across every allele group regardless of the allele on the "
        "row (openvax/topiary#232). Per-allele scores currently match "
        "allele_free_evidence=%s.",
        policy.name if policy.limit is None
        else "%s%s" % (POLICY_TOP_PREFIX, policy.limit),
        POLICY_ALL)
