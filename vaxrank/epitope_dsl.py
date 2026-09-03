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

"""Topiary DSL node construction for epitope filtering and scoring.

Two entry points:

- :func:`build_filter_node` — returns a boolean ``DSLNode`` (or ``None``) for
  the configured ``EpitopeConfig``; used with ``topiary.ranking.apply_filter``
  to drop rows wholesale before scoring.

- :func:`build_score_node` — returns a numeric ``DSLNode`` whose ``.eval(ctx)``
  produces a ``pd.Series`` indexed by exact prediction identity, peptide,
  peptide offset, and allele for current vaxrank frames.

When ``filter_expr`` / ``score_expr`` are set on the config, each is parsed
via :func:`topiary.ranking.parse`. When absent the nodes are built directly
from scalar fields (``binding_affinity_cutoff``, ``logistic_epitope_score_*``,
``scoring_mode``, ``percentile_rank_cutoff``). The default affinity-mode
score uses topiary 5.1's :class:`LogisticNormalizedExpr` so the output is in
``[0, 1]`` and matches the legacy
``logistic_epitope_score`` (pre-3.0)
byte-for-byte.

Multi-model inputs (e.g. LENS with both MHCflurry and netMHCpan producing
``pMHC_affinity``) are disambiguated via topiary 5.10's
``EvalContext(default_methods=...)`` and
``apply_filter(default_methods=...)``. Vaxrank computes the resolved
``default_methods`` dict by merging ``cfg.default_methods`` with an
auto-picked canonical (``mhcflurry`` > ``netmhcpan`` > alphabetical) for
Kinds the user hasn't specified, and forwards it to topiary.
"""

from __future__ import annotations

import functools
import logging

import pandas as pd


logger = logging.getLogger(__name__)


PREDICTION_ID_COLUMN = "prediction_id"

# Vaxrank works with two frame shapes, each with its own identity column.
# Both are named explicitly rather than left to topiary's inference: what
# inference keys on has changed twice (sample_name being prepended, then the
# group-index fixes in 5.17.1), and a frame whose grouping shifts underneath
# it merges candidates that are not the same candidate.
PREDICTION_GROUP_COLUMNS = (
    # Frames vaxrank builds from CandidateEpitope objects, and external
    # report readers. Identity is the complete provenance key.
    PREDICTION_ID_COLUMN,
    "peptide",
    "peptide_offset",
    "allele",
)
PIPELINE_GROUP_COLUMNS = (
    # Frames a TopiaryPredictor produces directly. There is one protein
    # fragment per variant, so the source sequence name *is* the identity;
    # these rows never carry a prediction_id.
    "source_sequence_name",
    "peptide",
    "peptide_offset",
    "allele",
)


# Method name -> topiary Kind for external-input predictions (LENS / pVACseq
# import). Listed methods come from the external report readers;
# anything not in this map is assumed to be a pMHC_affinity predictor.
_METHOD_KIND_MAP = {
    "mhcflurry": "pMHC_affinity",
    "netmhcpan": "pMHC_affinity",
    "netmhcstabpan": "pMHC_stability",
    # pVACseq aggregates across predictors without naming one; treat it as
    # a generic affinity predictor so default DSL nodes light up.
    "pvacseq": "pMHC_affinity",
}


@functools.lru_cache(maxsize=64)
def parse_epitope_expression(expr):
    """Parse a DSL string into a topiary node. Cached so the same config's
    ``filter_expr`` / ``score_expr`` aren't re-parsed multiple times within
    one :func:`write_neoepitope_report` call (validator + build_filter_node
    + build_score_node all need the AST). Parsed nodes are immutable, so
    sharing is safe.
    """
    from topiary.ranking import parse
    return parse(expr)


def prediction_kind_for_method(method_name):
    """Topiary Kind string for a prediction_method_name.

    Unknown methods default to pMHC_affinity so that a plain
    ``affinity.logistic(...)`` expression keeps working on external
    inputs even when the method name isn't in our registry.
    """
    if not method_name:
        return "pMHC_affinity"
    return _METHOD_KIND_MAP.get(str(method_name).lower(), "pMHC_affinity")



def epitopes_to_topiary_df(epitopes):
    """Convert a list of :class:`vaxrank.candidate_epitope.CandidateEpitope` into the
    topiary long-format DataFrame consumed by
    :class:`topiary.ranking.EvalContext` and
    :func:`topiary.ranking.apply_filter`.

    One row per leaf ``mhctools.Prediction`` in each CandidateEpitope's
    mutant context — the frame says what was predicted, and nothing else.
    Allele-free evidence stays one allele-free row; the patient's genotype
    reaches it through :func:`genotype_lookup`, which
    :func:`score_predictions` hands to topiary. When multiple predictions
    share a (peptide,
    allele) pair (e.g. MHCflurry and netMHCpan rows from a LENS file),
    DSL expressions can select between them via
    ``affinity['mhcflurry']`` / ``affinity['netmhcpan']``, and when
    the predictions carry ``predictor_version`` they can also
    disambiguate by version via ``affinity['mhcflurry', '2.1.1']``.
    Unqualified refs (``affinity.value``) resolve to the Kind's
    default method — see :func:`resolve_default_methods`.

    Comparator predictions (``epitope.wt``, ``nearest_self``, ...)
    are NOT emitted as rows. The DSL frame is mutant-only by
    convention — comparator data stays on the CandidateEpitope for display.
    """
    from .candidate_epitope import stated_or_blank

    rows = []
    for e in epitopes:
        ctx = e
        prediction_id = ctx.prediction_group_source
        predictions = ctx.predictions_flat()
        for p in predictions:
            rows.append({
                PREDICTION_ID_COLUMN: prediction_id,
                "source_sequence_name": ctx.source_sequence or ctx.sequence,
                "peptide": p.peptide,
                "peptide_offset": ctx.offset,
                "peptide_length": len(p.peptide),
                "allele": stated_or_blank(p.allele),
                "n_flank": ctx.n_flank,
                "c_flank": ctx.c_flank,
                "prediction_method_name": p.predictor_name,
                "predictor_version": stated_or_blank(
                    p.predictor_version),
                "kind": p.kind,
                "value": p.value,
                "affinity": p.value,
                "percentile_rank": p.percentile_rank,
                "score": p.score,
            })
    return pd.DataFrame(rows)


def prediction_group_columns(topiary_df):
    """Return the explicit score-group columns for a vaxrank frame.

    Picks by which identity column the frame carries: ``prediction_id`` for
    frames vaxrank built or an external reader produced, otherwise the
    predictor's own ``source_sequence_name``. A frame carrying neither has no
    identity to group on and is rejected rather than silently grouped by
    sequence, which would merge two biological sources whose peptides happen
    to be equal.
    """
    columns = set(topiary_df.columns)
    for candidate in (PREDICTION_GROUP_COLUMNS, PIPELINE_GROUP_COLUMNS):
        if set(candidate) <= columns:
            return list(candidate)
    raise ValueError(
        "Prediction frame carries no usable score-group identity. Expected "
        "%s (vaxrank-built or external-reader frames) or %s (frames from a "
        "TopiaryPredictor); got columns %s."
        % (", ".join(PREDICTION_GROUP_COLUMNS),
           ", ".join(PIPELINE_GROUP_COLUMNS),
           ", ".join(sorted(columns)) or "(none)"))



def resolve_default_methods(cfg, topiary_df):
    """Return the ``default_methods`` dict to pass to topiary's
    :class:`~topiary.ranking.EvalContext` / :func:`apply_filter`.

    The auto-pick — which model is "canonical" for a kind produced by
    several — is topiary's :func:`~topiary.ranking.resolve_default_methods`
    and its ``CANONICAL_METHOD_PREFERENCE``. vaxrank used to carry its own
    copy of that preference tuple, and the two had already drifted: ours
    ordered ``netmhcstabpan`` ahead of ``netmhcpan_el`` / ``netmhcpan_ba``
    where topiary orders it last, so the same frame could resolve to
    different models depending on which table was consulted.

    What stays here is the part topiary has no opinion about: an explicit
    ``cfg.default_methods`` entry is the user's stated choice and overrides
    the canonical pick. Kinds with a single model are omitted by topiary —
    the result only speaks where there is a real choice.
    """
    from topiary.ranking import resolve_default_methods as topiary_resolve

    if topiary_df.empty:
        return {}
    resolved = dict(topiary_resolve(topiary_df))
    configured = dict(cfg.default_methods or {})
    # Only announce a pick the user did not make. Logging it for a kind they
    # configured would tell them to set an entry they have already set.
    for kind, picked in sorted(resolved.items()):
        if kind in configured:
            continue
        logger.info(
            "Kind %s was produced by more than one model and no "
            "default_methods entry names one; topiary's canonical pick is "
            "%r. Set default_methods[%r] in the epitope config to override.",
            kind, picked, kind)
    # The user's explicit choice wins. An entry for a kind topiary did not
    # consider ambiguous is dropped: topiary only consults a default where
    # there is a real choice, so carrying it would change nothing, and
    # validate_default_methods has already rejected a name the data lacks.
    resolved.update({
        kind: method
        for kind, method in configured.items()
        if kind in resolved
    })
    return resolved


def validate_default_methods(cfg, topiary_df):
    """Error if ``cfg.default_methods`` names a kind or model absent from the
    data.

    Delegates to topiary's
    :func:`~topiary.ranking.validate_default_methods`, which checks the same
    thing for the same reason: an entry naming a model that never ran is
    silently inert, and stays inert until the day two models do produce that
    kind — at which point it starts deciding.
    """
    from topiary.ranking import validate_default_methods as topiary_validate

    if not cfg.default_methods or topiary_df.empty:
        return
    # Checked one entry at a time so the failure can name the exact config
    # key and value. topiary reports what is missing from the *data*, which
    # is the half it can know; which setting to edit is the half it cannot,
    # and is what the user actually needs in order to fix it.
    for kind, method in dict(cfg.default_methods).items():
        try:
            topiary_validate(topiary_df, {kind: method})
        except ValueError as error:
            raise ValueError(
                "default_methods[%r]=%r is not usable: %s Name a method "
                "present in the loaded predictions, or remove the entry."
                % (kind, method, error)) from error


def resolve_default_versions(cfg, topiary_df):
    """Return the ``default_versions`` mapping for an ambiguous frame.

    topiary raises on an unqualified reference when a model appears at more
    than one version. :func:`topiary.resolve_default_versions` picks one per
    ``(kind, model)`` — newest by PEP 440 — which is the same rule
    ``CandidateEpitope`` uses for unqualified access, so the DSL and the
    object API agree.

    vaxrank carried its own version of this while topiary lacked one
    (openvax/topiary#214). Two reasons the delegation is not merely
    de-duplication:

    * It resolves without dropping rows. The interim implementation narrowed
      the frame to the winning version, which silently removed the losing
      rows — so an expression pinning the older version in the same run
      found nothing.
    * It agrees with topiary about what counts as a version being named at
      all. The interim compared version strings, so a frame mixing one real
      version with rows recording none as the literal text ``"nan"`` treated
      that as a second version and dropped those rows.

    ``cfg`` is unused today — no setting pins a version, so there is nothing
    to layer over topiary's answer the way ``resolve_default_methods``
    layers ``epitopes.default_methods`` over its canonical pick. It is kept
    so the three resolvers read the same and so that setting, when it
    arrives, has an obvious home.
    """
    from topiary import describe_default_versions, resolve_default_versions

    if topiary_df.empty:
        return {}
    resolved = dict(resolve_default_versions(topiary_df))
    candidates = describe_default_versions(topiary_df)
    for (kind, method), version in sorted(resolved.items()):
        logger.warning(
            "%s reports %s at versions %s; scoring an unqualified reference "
            "with %s (newest). Name a version in the expression to pin one — "
            "every version stays available.",
            method, kind, ", ".join(candidates[(kind, method)]), version)
    return resolved


def genotype_map(epitopes):
    """``{prediction_group_key: genotype}`` for a list of CandidateEpitopes.

    The identity is the group columns other than ``allele``, in order, so
    this joins to a topiary frame built from the same epitopes. Shared by
    :func:`genotype_lookup`, which tells topiary which allele groups to
    create, and by allele attribution, which decides which of them get
    credited with allele-free evidence — one mapping so the two cannot
    disagree about a peptide's genotype.
    """
    return {
        e.prediction_group_key: tuple(e.patient_alleles or ())
        for e in epitopes or ()}


def genotype_lookup(epitopes, group_columns):
    """Return a per-peptide genotype callable for topiary's ``alleles=``.

    Evidence that describes the peptide rather than a peptide-MHC pair
    carries no allele, so it forms a group keyed on the empty string and
    contributes to no per-allele score. topiary closes that by taking the
    genotype each peptide should be evaluated against and adding the groups
    itself (openvax/topiary#182, #219).

    The genotype has to be *per peptide*, not per frame. LENS emits one row
    per (peptide, allele) that passed its own threshold, so each candidate
    arrives having been reported against only some of the patient's
    alleles — every LENS fixture in this repo holds between two and eight
    distinct allele sets. Declaring the union instead would add groups
    pairing a peptide with alleles it was never scored against, and any
    expression reading only peptide-level evidence would put a real number
    in every one of them.

    Returns ``None`` when no epitope declares a genotype, which leaves
    topiary to take the groups from the rows alone.
    """
    identity_column = group_columns[0]
    by_identity = genotype_map(epitopes)
    if not any(by_identity.values()):
        return None

    def lookup(keys):
        # An identity absent from the map declares nothing rather than
        # inheriting another candidate's genotype: the frame may have been
        # built elsewhere (pVACseq passes its own), and borrowing a genotype
        # is the failure this signature exists to prevent.
        return by_identity.get(
            (keys[identity_column], keys["peptide"], keys["peptide_offset"]),
            ())

    return lookup


def score_predictions(epitopes, cfg, *, topiary_df=None,
                      kind_support=None):
    """Score external-input epitopes using the configured Topiary DSL.

    Returns a ``pandas.Series`` of per-(peptide, allele) scores, indexed by
    ``(prediction_id, peptide, peptide_offset, allele)``.

    Multi-model data (e.g. both MHCflurry and netMHCpan for
    ``pMHC_affinity``) is handled via topiary's
    ``EvalContext(default_methods=...)`` and the same context for filtering:
    unqualified DSL references
    resolve to the per-Kind default (user-specified via
    ``cfg.default_methods`` or auto-picked canonical), while bracketed
    references like ``Affinity['netmhcpan']`` still pick up their
    specific rows.

    When ``cfg.score_expr`` / ``cfg.filter_expr`` are unset this applies
    the same default node that the main pipeline uses, which matches the
    legacy per-prediction logistic score byte-for-byte.

    Pass ``topiary_df`` to reuse an already-built frame (see
    :func:`epitopes_to_topiary_df`) rather than rebuilding from
    ``epitopes``.

    ``kind_support`` is topiary's per-(model, kind) MHC context, used to
    resolve ``mhc_dependence`` guards. Callers holding a
    ``TopiaryPredictor`` should forward its property; external report
    readers have no predictor and correctly leave it ``None``, in which
    case topiary resolves dependence from the rows.
    """
    from topiary.ranking import EvalContext, apply_filter

    df = (epitopes_to_topiary_df(epitopes)
          if topiary_df is None else topiary_df)
    if df.empty:
        return pd.Series(dtype=float)

    validate_default_methods(cfg, df)
    resolved = resolve_default_methods(cfg, df)
    resolved_versions = resolve_default_versions(cfg, df)
    group_columns = prediction_group_columns(df)
    alleles = genotype_lookup(epitopes, group_columns)

    filter_node = build_filter_node(cfg)
    if filter_node is not None:
        df = apply_filter(
            df, filter_node,
            group_keys=group_columns,
            default_methods=resolved,
            default_versions=resolved_versions,
            alleles=alleles,
            kind_support=kind_support)
        if df.empty:
            return pd.Series(dtype=float)

    score_node = build_score_node(cfg)
    ctx = EvalContext(
        df,
        group_keys=group_columns,
        default_methods=resolved,
        default_versions=resolved_versions,
        alleles=alleles,
        kind_support=kind_support,
    )
    return (
        score_node.eval(ctx)
        .reindex(ctx.group_index)
        .fillna(0.0)
    )


def attach_per_allele_scores(epitopes, cfg=None, *, topiary_df=None,
                             kind_support=None):
    """Score ``epitopes`` via the configured DSL and return a new list of
    :class:`~vaxrank.candidate_epitope.CandidateEpitope` instances with
    each one's ``per_allele_scores`` populated.

    External report readers (``read_pvacseq_report`` / ``read_lens_report``) build
    unscored CandidateEpitopes from a TSV. The upstream pipeline path
    populates ``per_allele_scores`` inside ``predict_epitopes`` from the
    DSL eval; this helper is the equivalent step for the loader path so
    downstream consumers (``VaccinePeptide.target_epitope_score`` etc.)
    see the same shape regardless of where the epitopes came from.

    ``cfg`` defaults to a fresh :class:`~vaxrank.epitope_config.EpitopeConfig`
    so the loader path's default scoring matches the legacy pre-3.1
    semantics exactly. Pass a custom config to override the threshold
    knobs or supply a user ``score_expr``.

    ``topiary_df`` can be supplied by loaders that already have a
    richer topiary frame than can be reconstructed from CandidateEpitope
    objects alone, e.g. pVACseq passthrough annotation columns.
    """
    from dataclasses import replace
    from .candidate_epitope import stated_or_blank
    from .epitope_config import EpitopeConfig

    from .allele_evidence import apply_attribution, attribute_frame

    if not epitopes:
        return epitopes
    if cfg is None:
        cfg = EpitopeConfig()
    # Build the frame once and share it: attribution and scoring must see
    # the same rows, and rebuilding is not free.
    if topiary_df is None:
        topiary_df = epitopes_to_topiary_df(epitopes)

    # Validate before scoring, not only in the report writer. topiary
    # raises on a bad reference too, with a good "did you mean" list, but
    # mid-run and without naming the setting — and the reader path scores
    # long before any report is written, so that is where a user meets the
    # error (#388).
    validate_dsl_against_predictions(cfg, epitopes, topiary_df=topiary_df)

    # Attribution runs on the unfiltered frame, so the policy ranks on the
    # evidence the input carried rather than on whatever survived a cutoff:
    # a filter on affinity would otherwise change which allele is credited
    # with a peptide's processing score.
    #
    # And it runs before scoring, because a narrowing policy rewrites the
    # frame the scores are computed from — that is what makes it a policy
    # about scoring rather than a note attached to the result.
    policy = cfg.allele_policy
    genotypes = genotype_map(epitopes)
    group_columns = prediction_group_columns(topiary_df)
    attributions = attribute_frame(
        topiary_df, policy,
        genotype_for=lambda key: genotypes.get(key, ()),
        default_methods=resolve_default_methods(cfg, topiary_df),
        group_columns=group_columns)
    scoring_df = apply_attribution(
        topiary_df, attributions, policy, group_columns=group_columns)

    score_series = score_predictions(
        epitopes, cfg, topiary_df=scoring_df, kind_support=kind_support)
    # score_series is keyed by prediction ID, peptide, offset, and allele.
    by_position: dict[tuple, dict[str, float]] = {}
    for idx, val in score_series.items():
        src_name, peptide, offset, allele = idx
        # Allele-independent predictions use a blank allele. They belong in
        # the nested prediction model, but never in a map whose public
        # contract is explicitly per patient allele.
        #
        # stated_or_blank rather than a truthiness test: this index comes
        # from the frame, which a caller may have built (pVACseq passes its
        # own), and topiary collapses a stringified null in a group key to a
        # real null — which is truthy. A bare ``if allele`` let NaN through
        # as a patient allele.
        allele = stated_or_blank(allele)
        if allele:
            by_position.setdefault(
                (src_name, peptide, offset), {})[allele] = float(val)
    return [
        replace(
            e,
            per_allele_scores=by_position.get(e.prediction_group_key, {}),
            allele_attributions=attributions.get(
                e.prediction_group_key, ()))
        for e in epitopes
    ]


def epitopes_for_ranking(epitopes, cfg=None):
    """Materialize external prediction groups eligible for ranking.

    :func:`attach_per_allele_scores` records one key for every retained
    peptide-allele group. This function applies the configured
    ``min_epitope_score`` to those retained groups and returns ranking-only
    copies containing the eligible allele-scoped mutant and comparator
    leaves. Groups below the minimum remain on the input objects for audit
    reporting but cannot become vaccine targets.

    The input objects remain unchanged and retain their complete
    ``patient_alleles`` membership. Callers must infer the patient genotype
    from that pre-filter collection; filtering target evidence must never
    redefine which alleles the patient has.
    """
    from dataclasses import replace
    from .epitope_config import EpitopeConfig

    if cfg is None:
        cfg = EpitopeConfig()

    retained_epitopes = []
    for epitope in epitopes:
        retained_alleles = {
            allele
            for allele, score in epitope.per_allele_scores.items()
            if score >= cfg.min_epitope_score
        }
        if not retained_alleles:
            continue
        retained_predictions = tuple(
            prediction
            for prediction in epitope.predictions_flat()
            if not prediction.allele
            or prediction.allele in retained_alleles
        )
        if not retained_predictions:
            continue
        retained_comparators = {}
        for name, comparator in epitope.comparators.items():
            retained_comparators[name] = replace(
                comparator,
                predictions=tuple(
                    prediction
                    for prediction in comparator.predictions_flat()
                    if not prediction.allele
                    or prediction.allele in retained_alleles
                ),
            )
        retained_epitopes.append(replace(
            epitope,
            predictions=retained_predictions,
            comparators=retained_comparators,
            patient_alleles=tuple(sorted(retained_alleles)),
            per_allele_scores={
                allele: epitope.per_allele_scores[allele]
                for allele in sorted(retained_alleles)
            },
        ))
    return retained_epitopes


def collect_dsl_references(node):
    """Walk a parsed DSL node and collect its external-data references.

    Returns a dict with two keys:

      ``columns`` - set of column names referenced via ``Column(...)``
      ``kinds``   - set of (kind, method, version) tuples from ``Field`` nodes

    ``method`` and ``version`` may be ``None`` when the formula used a
    kind without qualification (``affinity.value``). Column refs are
    already validated by topiary's ``apply_filter``; vaxrank only uses
    this to drive :func:`validate_dsl_against_predictions`.
    """
    from topiary.ranking.nodes import Field

    columns = set()
    kinds = set()
    stack = [node]
    while stack:
        n = stack.pop()
        if n is None:
            continue
        col_name = getattr(n, "col_name", None)
        if isinstance(col_name, str):
            columns.add(col_name)
        if isinstance(n, Field):
            kinds.add((n.kind, n.method, n.version))
        child_iter = getattr(n, "child_nodes", None)
        if callable(child_iter):
            stack.extend(child_iter())
    return {"columns": columns, "kinds": kinds}


def validate_dsl_against_predictions(cfg, epitopes, *, topiary_df=None):
    """Error early when ``filter_expr`` / ``score_expr`` references a
    predictor the loaded epitopes don't expose.

    Topiary's own ``apply_filter`` already validates ``Column(...)``
    references, but when a ``Field``'s (kind, method, version) isn't
    present in the data topiary can silently evaluate to NaN for that
    sub-expression. This function catches that up front and points the
    user at ``default_methods`` / ``--input-lens`` when appropriate.

    Pass ``topiary_df`` to reuse an already-built frame rather than
    rebuilding it from ``epitopes``.
    """
    # Paired with the setting each came from, so an error can name the key
    # to edit. topiary reports what the data lacks, which is the half it
    # can know; which config key produced the reference is the half it
    # cannot, and is the half the user needs.
    sources = []
    if cfg.filter_expr is not None:
        sources.append(("epitopes.filter_expr", cfg.filter_expr,
                        parse_epitope_expression(cfg.filter_expr)))
    if cfg.score_expr is not None:
        sources.append(("epitopes.score_expr", cfg.score_expr,
                        parse_epitope_expression(cfg.score_expr)))
    if not sources:
        return

    df = (epitopes_to_topiary_df(epitopes)
          if topiary_df is None else topiary_df)
    available_kinds = set(df["kind"].unique()) if not df.empty else set()
    available_methods = (
        set(df["prediction_method_name"].dropna().unique())
        if not df.empty else set())
    from topiary import is_named_version

    # Drop the "no version known" sentinels so a formula pinning a specific
    # version against predictions that don't carry any version errors out,
    # and so the error message's version list matches what we actually
    # checked against.
    available_versions = set()
    if not df.empty:
        available_versions = {
            v for v in df["predictor_version"].dropna().unique()
            if is_named_version(v)
        }

    available_columns = set(df.columns)

    for setting, expr, node in sources:
        refs = collect_dsl_references(node)
        # Column references. topiary raises on these too, at eval time and
        # with a good "did you mean" list — but mid-run rather than when
        # the config is read, and without naming the setting. Which columns
        # exist depends on the source: pVACseq's aggregated flavour carries
        # no transcript_expression where all_epitopes does, so one export can
        # make a working config and another a clear error (#388).
        for column in sorted(refs["columns"]):
            if column in available_columns:
                continue
            raise ValueError(
                f"{setting} references column {column!r}, which the loaded "
                f"input does not carry. Expression: {expr!r}. Columns "
                f"available on this input: {sorted(available_columns)}. "
                f"Which columns exist depends on the source and its "
                f"flavour, so an expression that works on one export may "
                f"not on another."
            )
        for kind, method, version in refs["kinds"]:
            if kind not in available_kinds:
                raise ValueError(
                    f"DSL expression references kind {kind!r} but loaded "
                    f"predictions only expose kinds {sorted(available_kinds)}. "
                    f"Check filter_expr / score_expr against the LENS "
                    f"file columns or input source."
                )
            if method is not None and method not in available_methods:
                raise ValueError(
                    f"DSL expression references predictor {method!r} but "
                    f"loaded predictions only expose methods "
                    f"{sorted(available_methods)}. Remove the bracket "
                    f"reference or load a data source that exposes this "
                    f"method."
                )
            if version is not None and version not in available_versions:
                raise ValueError(
                    f"DSL expression references predictor version "
                    f"{version!r} but loaded predictions only expose "
                    f"versions {sorted(available_versions)}. "
                    f"Drop the version from the bracket expression or "
                    f"update the LENS file / loader."
                )


# Default per-prediction score expressions. These string templates are
# the SINGLE source of truth for what the default scoring formula looks
# like — there is no parallel Python builder. When an EpitopeConfig has
# ``score_expr=None`` the relevant template is formatted with the
# config's scalar threshold knobs and then parsed by the topiary DSL,
# exactly as a user-supplied ``score_expr`` would be. Format placeholders
# correspond to EpitopeConfig fields; keep them in sync.
DEFAULT_AFFINITY_SCORE_EXPR_TEMPLATE = (
    "(affinity.value < {binding_affinity_cutoff}) * "
    "affinity.value.logistic_normalized("
    "{logistic_epitope_score_midpoint}, {logistic_epitope_score_width})"
)
DEFAULT_PERCENTILE_SCORE_EXPR_TEMPLATE = (
    "(percentile_rank < {percentile_rank_cutoff}) * "
    "(1.0 - percentile_rank / {percentile_rank_cutoff}).clip(0.0, none)"
)


def default_score_expr(cfg):
    """Return the default per-prediction score expression for ``cfg``.

    Formats :data:`DEFAULT_AFFINITY_SCORE_EXPR_TEMPLATE` or
    :data:`DEFAULT_PERCENTILE_SCORE_EXPR_TEMPLATE` against the
    config's scalar threshold fields. This is the string a user would
    write to reproduce the default — it is also exactly the string
    that gets parsed when ``cfg.score_expr`` is unset, so the default
    behavior travels through the same DSL pipeline as any override.
    """
    if cfg.scoring_mode == "percentile_rank":
        return DEFAULT_PERCENTILE_SCORE_EXPR_TEMPLATE.format(
            percentile_rank_cutoff=cfg.percentile_rank_cutoff)
    return DEFAULT_AFFINITY_SCORE_EXPR_TEMPLATE.format(
        binding_affinity_cutoff=cfg.binding_affinity_cutoff,
        logistic_epitope_score_midpoint=cfg.logistic_epitope_score_midpoint,
        logistic_epitope_score_width=cfg.logistic_epitope_score_width,
    )


def build_filter_node(cfg):
    """Return the ``DSLNode`` used to filter predictions, or ``None`` for no-op.

    The default is no filter — legacy vaxrank never hard-dropped rows; it
    scored above-cutoff rows as 0 and then applied ``min_epitope_score`` as
    a gate. The default score node preserves that by multiplying by the
    cutoff mask. Users who want a hard pre-filter must set ``filter_expr``.
    """
    if cfg.filter_expr is not None:
        return parse_epitope_expression(cfg.filter_expr)
    return None


def build_score_node(cfg):
    """Return the ``DSLNode`` that computes the per-(peptide, allele) score.

    Single code path: parse the resolved expression string. When
    ``cfg.score_expr`` is set the user's string is used; otherwise
    :func:`default_score_expr` returns the canonical default formatted
    against the config's scalar threshold knobs. The default is not a
    different mechanism — it is the same DSL machinery applied to a
    well-known formula.
    """
    return parse_epitope_expression(
        cfg.score_expr if cfg.score_expr is not None else default_score_expr(cfg)
    )
