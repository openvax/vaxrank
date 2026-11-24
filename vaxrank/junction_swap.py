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

"""Junction-aware linker selection for multi-antigen mRNA constructs.

Concatenating antigens with a fixed linker can create k-mers that span
the junction (end of antigen N + linker + start of antigen N+1) which
the patient's MHC can present even though the sequence was never
present in the wild proteome and is not the intended neoantigen.
These "chimeric" peptides are unintended candidate ligands; predicted
binding or presentation alone does not establish an immune response.

This module picks an inter-antigen linker per junction independently
to minimize predicted MHC presentation of chimeric k-mers. The
search space per junction is small (default 5 candidates) and the
junctions are independent under the assumption that linker choice
at junction j doesn't affect chimeric k-mers at junction k != j —
which is true given a length-cap on k that doesn't span two
junctions. So total work is N junctions × M candidates × K k-mer
lengths × |HLA alleles| predictions, all done up front.

Issue: openvax/vaxrank#247.
"""

import logging
import math
from dataclasses import dataclass, field
from numbers import Real

from topiary import is_stated, predict_peptide_occurrences

from .vaccine_library import JUNCTION_SWAP_CANDIDATES, get_linker

logger = logging.getLogger(__name__)


# Existing percentile-rank heuristics, lower = stronger prediction.
# Endpoint selection is recorded separately: affinity when available,
# otherwise presentation. These thresholds are selection rules, not
# probabilities or observations of presentation. The optimizer minimizes
# count of rank ≤ STRONG first, then rank ≤ MILD as tie-breaker, then
# maximizes the worst-case rank.
RANK_STRONG = 0.5
RANK_MILD = 2.0


@dataclass(frozen=True)
class MRNAJunctionOptimizationPolicy:
    """Explicit per-trial policy for joint mRNA antigen ordering, linker
    selection, and optional proteasomal cleavage optimization.

    All fields can be configured in a per-trial YAML policy under
    ``mrna.junction_optimization``. Default ``optimize_order=False``
    preserves fixed-order per-junction linker swapping for backward
    compatibility with saved designs.
    """
    optimize_order: bool = False
    exact_dp_max_antigens: int = 12
    include_signal_peptide_boundary: bool = True
    include_mitd_boundary: bool = True
    strong_burden_weight: float = 100.0
    mild_burden_weight: float = 10.0
    cleavage_penalty_weight: float = 5.0
    cleavage_scoring: str = 'none'
    pepsickle_human_only: bool = False
    pepsickle_threshold: float = 0.5
    linker_length_penalty: float = 0.01
    worst_rank_tiebreak_weight: float = 0.001
    max_2opt_passes: int = 50

    def __post_init__(self):
        for flag in (
                'optimize_order',
                'include_signal_peptide_boundary',
                'include_mitd_boundary',
                'pepsickle_human_only',
        ):
            if not isinstance(getattr(self, flag), bool):
                raise ValueError(f"junction_optimization.{flag} must be a boolean")
        if (isinstance(self.exact_dp_max_antigens, bool)
                or not isinstance(self.exact_dp_max_antigens, int)
                or self.exact_dp_max_antigens < 1
                or self.exact_dp_max_antigens > 20):
            raise ValueError(
                "junction_optimization.exact_dp_max_antigens must be an integer in [1, 20]")
        if (isinstance(self.max_2opt_passes, bool)
                or not isinstance(self.max_2opt_passes, int)
                or self.max_2opt_passes < 0):
            raise ValueError("junction_optimization.max_2opt_passes must be an integer >= 0")
        if self.cleavage_scoring not in ('none', 'pepsickle'):
            raise ValueError(
                "junction_optimization.cleavage_scoring must be 'none' or 'pepsickle'")
        if (isinstance(self.pepsickle_threshold, bool)
                or not isinstance(self.pepsickle_threshold, Real)
                or not math.isfinite(self.pepsickle_threshold)
                or not 0.0 <= float(self.pepsickle_threshold) <= 1.0):
            raise ValueError(
                "junction_optimization.pepsickle_threshold must be in [0.0, 1.0]")
        for name in (
                'strong_burden_weight',
                'mild_burden_weight',
                'cleavage_penalty_weight',
                'linker_length_penalty',
                'worst_rank_tiebreak_weight',
        ):
            val = getattr(self, name)
            if (isinstance(val, bool)
                    or not isinstance(val, Real)
                    or not math.isfinite(val)
                    or float(val) < 0.0):
                raise ValueError(
                    f"junction_optimization.{name} must be a finite number >= 0")

    @classmethod
    def resolve(cls, value=None):
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise TypeError(
                f"junction_optimization must be a mapping or MRNAJunctionOptimizationPolicy; "
                f"got {type(value).__name__}")
        defaults = cls().to_dict()
        unknown = set(value) - set(defaults)
        if unknown:
            raise ValueError(
                "Unknown mrna.junction_optimization fields: %s"
                % ", ".join(sorted(unknown)))
        merged = {k: (defaults[k] if value.get(k) is None else value[k]) for k in defaults}
        return cls(**merged)

    def to_dict(self):
        return {
            'optimize_order': self.optimize_order,
            'exact_dp_max_antigens': self.exact_dp_max_antigens,
            'include_signal_peptide_boundary': self.include_signal_peptide_boundary,
            'include_mitd_boundary': self.include_mitd_boundary,
            'strong_burden_weight': float(self.strong_burden_weight),
            'mild_burden_weight': float(self.mild_burden_weight),
            'cleavage_penalty_weight': float(self.cleavage_penalty_weight),
            'cleavage_scoring': self.cleavage_scoring,
            'pepsickle_human_only': self.pepsickle_human_only,
            'pepsickle_threshold': float(self.pepsickle_threshold),
            'linker_length_penalty': float(self.linker_length_penalty),
            'worst_rank_tiebreak_weight': float(self.worst_rank_tiebreak_weight),
            'max_2opt_passes': self.max_2opt_passes,
        }


@dataclass
class JunctionPrediction:
    """A single chimeric k-mer prediction at a junction."""
    junction_index: int
    linker_name: str
    kmer: str
    allele: str
    rank: float


@dataclass
class JunctionSwapResult:
    """Outcome of the per-junction linker swap or joint order+linker optimization."""
    chosen_linker_per_junction: list  # list[Linker]
    burden: int  # count of rank <= RANK_MILD across all junctions, post-swap
    strong_burden: int  # count of rank <= RANK_STRONG, post-swap
    chimeric_predictions: list = field(default_factory=list)  # all post-swap
    prediction_metadata: dict = field(default_factory=dict)
    antigen_order: list = field(default_factory=list)
    reordered: bool = False
    solver: str = 'fixed_order'
    fixed_order_burden: int | None = None
    fixed_order_strong_burden: int | None = None
    signal_peptide_burden: int = 0
    signal_peptide_strong_burden: int = 0
    cleavage_penalty: float = 0.0
    default_cleavage_penalty: float = 0.0
    cleavage_status: str = 'disabled'
    total_cost: float | None = None
    default_total_cost: float | None = None
    policy: dict | None = None

    def linker_names(self):
        return [link.name for link in self.chosen_linker_per_junction]


def boundary_kmers(left_aa, right_aa, k_lengths, reference_proteome=None):
    """Enumerate k-mers spanning a direct (zero-linker) boundary ``left_aa|right_aa``.

    A k-mer spans the boundary if it includes at least one residue from
    ``left_aa`` and at least one residue from ``right_aa``. Used for the
    ``signal_peptide -> first_antigen`` boundary.
    """
    if not left_aa or not right_aa:
        return []
    max_k = max(k_lengths)
    left_ctx = left_aa[-(max_k - 1):] if len(left_aa) >= max_k - 1 else left_aa
    right_ctx = right_aa[:max_k - 1] if len(right_aa) >= max_k - 1 else right_aa
    full = left_ctx + right_ctx
    split = len(left_ctx)

    out = []
    for k in k_lengths:
        for i in range(len(full) - k + 1):
            j = i + k
            if i < split < j:
                kmer = full[i:j]
                if reference_proteome is not None and kmer in reference_proteome:
                    continue
                out.append(kmer)
    return out


def junction_kmers(left_aa, linker_aa, right_aa, k_lengths,
                   reference_proteome=None):
    """Enumerate k-mers spanning a single junction.

    A k-mer "spans" the junction if it touches at least one residue of
    the linker. Pure k-mers entirely inside ``left_aa`` or ``right_aa``
    are excluded — they're real antigen content, not chimeric.

    ``reference_proteome`` is an optional callable or container that
    answers ``kmer in reference_proteome`` so already-tolerated
    sequences (e.g. self-peptides happening to match a junction) are
    filtered out.
    """
    max_k = max(k_lengths)
    # Limit context to (max_k - 1) on each side so the only k-mers we
    # produce are junction-spanning by construction.
    left_ctx = left_aa[-(max_k - 1):] if len(left_aa) >= max_k - 1 else left_aa
    right_ctx = right_aa[:max_k - 1] if len(right_aa) >= max_k - 1 else right_aa
    full = left_ctx + linker_aa + right_ctx
    linker_start = len(left_ctx)
    linker_end = linker_start + len(linker_aa)

    out = []
    for k in k_lengths:
        for i in range(len(full) - k + 1):
            j = i + k
            # Junction-spanning: must touch the linker, not entirely
            # inside left_ctx or right_ctx.
            if j <= linker_start or i >= linker_end:
                continue
            kmer = full[i:j]
            if reference_proteome is not None and kmer in reference_proteome:
                continue
            out.append(kmer)
    return out


def junction_rank_kind(predictor, rank_kind='auto'):
    """Resolve one declared allele-specific rank endpoint before inference."""
    if rank_kind not in ('auto', 'pMHC_affinity', 'pMHC_presentation'):
        raise ValueError('Junction rank kind must be auto, pMHC_affinity or pMHC_presentation')
    support = {str(getattr(kind, 'value', kind)): spec
               for kind, spec in predictor.kind_support().items()}
    selected = ('pMHC_affinity' if 'pMHC_affinity' in support else 'pMHC_presentation'
                ) if rank_kind == 'auto' else rank_kind
    if support.get(selected, {}).get('mhc_dependence') != 'single_allele':
        raise ValueError('Junction ranking requires a declared single-allele %s endpoint' % selected)
    return selected


def score_junction_kmers(kmers, alleles, predictor, models=None, *, rank_kind='auto'):
    """Predict exact peptide occurrences through Topiary and select one endpoint.

    ``auto`` explicitly preserves the historical affinity-first preference.
    The selected endpoint must have one model/version and a finite percentile
    rank for every requested occurrence and allele. Other endpoints remain in
    model metadata but never contribute to the rank burden. Missing coverage
    raises instead of making a linker appear to have zero burden.
    """
    if not kmers:
        return []
    from .input_scope import normalize_alleles
    selected = junction_rank_kind(predictor, rank_kind)
    requested = normalize_alleles(alleles, 'Junction prediction')
    if not requested or set(requested) != set(normalize_alleles(
            predictor.alleles, 'Junction predictor')):
        raise ValueError('Junction predictor allele coverage must match the requested genotype')
    frame = predict_peptide_occurrences([
        dict(prediction_id=str(i), peptide=peptide) for i, peptide in enumerate(kmers)
    ], predictor, use_flanks=False)
    if models is not None:
        for record in frame.to_dict('records'):
            model = {target: (record.get(source) if is_stated(record.get(source)) else None)
                     for target, source in (
                         ('predictor_name', 'prediction_method_name'),
                         ('predictor_version', 'predictor_version'), ('kind', 'kind'))}
            if model not in models:
                models.append(model)
    chosen = frame.loc[frame.kind.eq(selected)]
    if len(chosen[['prediction_method_name', 'predictor_version']].drop_duplicates()) != 1:
        raise ValueError('Junction ranks require one model/version for the selected endpoint')
    if not chosen.percentile_rank.map(
            lambda value: is_stated(value) and math.isfinite(float(value))
            and 0 <= float(value) <= 100).all():
        raise ValueError('Junction ranking requires a finite percentile_rank in [0, 100] '
                         'for every requested peptide and allele')
    return [(row.peptide, row.allele, float(row.percentile_rank))
            for row in chosen.itertuples()]


def _burden_key(rows, rank_strong=RANK_STRONG, rank_mild=RANK_MILD):
    """Tuple of (strong_count, mild_count, worst_rank_negated) used to
    pick the linker with lowest junction burden. Lower = better.

    Tie-breaking order:
      1. Fewer rank ≤ rank_strong predictions
      2. Fewer rank ≤ rank_mild hits
      3. Higher worst-case rank (push the worst hit's rank as high as
         possible — i.e. minimize negated worst-rank)

    The 100.0 default for empty-rows is the upper bound of the
    percentile-rank range (0..100). Empty rows are allowed only when no
    junction peptides were requested; missing predictions fail upstream.
    """
    strong = sum(1 for r in rows if r[2] <= rank_strong)
    mild = sum(1 for r in rows if r[2] <= rank_mild)
    worst = min((r[2] for r in rows), default=100.0)
    return (strong, mild, -worst)


def optimize_linkers(
        antigen_aas, alleles, predictor,
        candidate_names=JUNCTION_SWAP_CANDIDATES,
        k_lengths=(8, 9, 10, 11),
        reference_proteome=None,
        rank_strong=RANK_STRONG, rank_mild=RANK_MILD,
        rank_kind='auto',
        mitd_aa=None,
        default_linker_name=None):
    """Pick the candidate linker at each junction that minimizes
    predicted MHC presentation of chimeric k-mers.

    Parameters
    ----------
    antigen_aas : list[str]
        Antigen amino-acid sequences in the order they'll be
        concatenated.
    alleles : list[str]
        Patient HLA alleles, e.g. ``['HLA-A*02:01', 'HLA-B*07:02']``.
    predictor : mhctools.BasePredictor
        An MHC predictor with a ``predict(peptides)`` method returning
        MHCtools ``PeptideResult`` objects. The predictor's
        configured alleles must include those in ``alleles``.
    candidate_names : iterable of str
        Linker names parseable by ``vaccine_library.get_linker``.
        Default is ``JUNCTION_SWAP_CANDIDATES`` = G3S / G4S /
        (G3S)2 / (G4S)2 / AAA.
    k_lengths : iterable of int
        MHC-I k-mer lengths to enumerate (default 8–11).
    reference_proteome : optional, container-like
        ``in``-checkable set or set-of-kmers index. When provided,
        chimeric k-mers that already occur in the patient's reference
        proteome are dropped (already tolerated, not new presentation).
    mitd_aa : optional, str
        When provided, also optimize the last-antigen ↔ MITD junction.
        The result's ``chosen_linker_per_junction`` will have
        ``len(antigen_aas)`` entries (one per inter-antigen junction
        plus one for the pre-MITD junction). Without ``mitd_aa``, the
        result has ``len(antigen_aas) - 1`` entries.
    default_linker_name : optional, str
        If supplied, the optimizer also tracks the burden of using
        the default linker at every junction so the caller can
        decide whether to swap. The default's burden is returned in
        ``JunctionSwapResult.default_burden`` /
        ``default_strong_burden`` (added when this argument is set).
        Including the default in ``candidate_names`` is sufficient if
        it's already there — this argument lets callers track the
        baseline without duplicating predictor calls.

    Returns
    -------
    JunctionSwapResult
    """
    alleles = tuple(alleles or ())
    k_lengths = tuple(k_lengths)
    if (not all(isinstance(value, Real) and not isinstance(value, bool)
                and math.isfinite(value) for value in (rank_strong, rank_mild))
            or not 0 <= rank_strong <= rank_mild <= 100):
        raise ValueError('Junction thresholds require 0 <= strong <= mild <= 100')
    if len(antigen_aas) == 0:
        return JunctionSwapResult(
            chosen_linker_per_junction=[],
            burden=0, strong_burden=0, chimeric_predictions=[])

    # Build the "junctions" list: (left_aa, right_aa) tuples. For pure
    # inter-antigen junctions there are len(antigens)-1 entries. With
    # mitd_aa, append one more for the last-antigen → MITD junction.
    junctions = []
    for j in range(len(antigen_aas) - 1):
        junctions.append((antigen_aas[j], antigen_aas[j + 1]))
    if mitd_aa and antigen_aas:
        junctions.append((antigen_aas[-1], mitd_aa))

    if not junctions:
        return JunctionSwapResult(
            chosen_linker_per_junction=[],
            burden=0, strong_burden=0, chimeric_predictions=[])

    candidates = [get_linker(n) for n in candidate_names]
    default_idx_in_candidates = None
    if default_linker_name is not None:
        target = default_linker_name.upper()
        for idx, c in enumerate(candidates):
            if c.name.upper() == target:
                default_idx_in_candidates = idx
                break
        if default_idx_in_candidates is None:
            # Default isn't in the candidate set; add it so we can
            # measure its burden in the same sweep.
            candidates.append(get_linker(default_linker_name))
            default_idx_in_candidates = len(candidates) - 1

    chosen = []
    all_predictions = []
    total_burden = 0
    total_strong = 0
    default_total_burden = 0
    from topiary import TopiaryPredictor, describe_prediction_models
    configured = predictor if isinstance(predictor, TopiaryPredictor) else TopiaryPredictor(models=[predictor])
    prediction_metadata = {
        'configured_models': describe_prediction_models(configured),
        'predictor_class': type(predictor).__module__ + '.' + type(predictor).__name__,
        'alleles': list(alleles),
        'peptide_lengths': list(k_lengths),
        'rank_kind_policy': ('affinity_if_available_else_presentation'
                             if rank_kind == 'auto' else 'explicit'),
        'rank_kind': junction_rank_kind(predictor, rank_kind),
        'rank_metric': 'percentile_rank',
        'rank_strong': rank_strong,
        'rank_mild': rank_mild,
        'context': 'peptide_only',
        'coverage': 'complete',
        'models': [],
        'requests': [],
        'scored_prediction_count': 0,
        'junction_decisions': [],
    }
    default_total_strong = 0

    for j, (left_aa, right_aa) in enumerate(junctions):
        best = None
        per_cand_keys = []
        for cand in candidates:
            kmers = junction_kmers(
                left_aa, cand.amino_acids, right_aa, k_lengths,
                reference_proteome=reference_proteome)
            if kmers:
                prediction_metadata['requests'].append({
                    'junction_index': j, 'linker_name': cand.name, 'peptides': list(kmers)})
            rows = score_junction_kmers(
                kmers, alleles, predictor, models=prediction_metadata['models'], rank_kind=rank_kind)
            prediction_metadata['scored_prediction_count'] += len(rows)
            key = _burden_key(rows, rank_strong, rank_mild)
            per_cand_keys.append(key)
            if best is None or key < best[0]:
                best = (key, cand, rows)
        key, cand, rows = best
        prediction_metadata['junction_decisions'].append(dict(
            junction_index=j, preferred_linker=cand.name,
            ordering='strong_count, mild_count, negative_minimum_rank, candidate_order',
            candidates=[dict(linker=candidate.name, key=list(key), candidate_order=index)
                        for index, (candidate, key) in enumerate(zip(candidates, per_cand_keys))]))
        strong, mild, _ = key
        chosen.append(cand)
        total_burden += mild
        total_strong += strong
        for kmer, allele, rank in rows:
            all_predictions.append(JunctionPrediction(
                junction_index=j, linker_name=cand.name,
                kmer=kmer, allele=allele, rank=rank))
        if default_idx_in_candidates is not None:
            d_key = per_cand_keys[default_idx_in_candidates]
            default_total_strong += d_key[0]
            default_total_burden += d_key[1]
        logger.info(
            "Junction %d: chose linker %s (strong=%d, mild=%d).",
            j, cand.name, strong, mild)

    result = JunctionSwapResult(
        chosen_linker_per_junction=chosen,
        burden=total_burden,
        strong_burden=total_strong,
        chimeric_predictions=all_predictions,
        prediction_metadata=prediction_metadata,
    )
    if default_linker_name is not None:
        # Attach the default's burden so the caller can decide whether
        # to swap without a second predictor sweep.
        result.default_burden = default_total_burden
        result.default_strong_burden = default_total_strong
    return result


def _edge_composite_cost(rows, linker_aa, cleavage_penalty, policy,
                         rank_strong=RANK_STRONG, rank_mild=RANK_MILD):
    """Compute scalar objective cost and burden tuple for a directed junction."""
    strong, mild, neg_worst = _burden_key(rows, rank_strong, rank_mild)
    worst = -neg_worst
    cost = (
        strong * float(policy.strong_burden_weight)
        + mild * float(policy.mild_burden_weight)
        + float(cleavage_penalty) * float(policy.cleavage_penalty_weight)
        + len(linker_aa) * float(policy.linker_length_penalty)
        + (100.0 - worst) * float(policy.worst_rank_tiebreak_weight)
    )
    return cost, strong, mild, worst


def _compute_pepsickle_boundary_penalties(junction_contexts, policy):
    """Batch-score C-terminal proteasomal cleavage recovery via Pepsickle.

    Each item in ``junction_contexts`` is ``(key, left_aa, linker_aa, right_aa)``.
    For each junction, we construct a local sequence window around the C-terminal
    bond of ``left_aa`` (up to 16 aa upstream and 16 aa downstream so Pepsickle's
    8-aa context window is unpadded whenever sequence length allows) and measure
    cleavage probability at bond ``len(left_ctx)``. Penalty is ``1.0 - prob``
    (lower penalty = better proteasomal release of ``left_aa``'s C-terminus).
    """
    if not junction_contexts:
        return {}, 'disabled'
    try:
        from mhctools import PeptideInput
        from .cleavage_inference import audit_pepsickle_inputs
    except Exception as exc:
        logger.warning("Pepsickle cleavage scoring unavailable: %s", exc)
        return {key: 0.0 for key, *_ in junction_contexts}, 'pepsickle_unassessed'

    seq_to_input = {}
    key_to_seq_and_bond = {}
    for key, left_aa, linker_aa, right_aa in junction_contexts:
        left_ctx = left_aa[-16:] if len(left_aa) > 16 else left_aa
        right_tail = (linker_aa + right_aa)[:16]
        if not left_ctx or not right_tail:
            key_to_seq_and_bond[key] = (None, None)
            continue
        seq = left_ctx + right_tail
        bond = len(left_ctx)
        key_to_seq_and_bond[key] = (seq, bond)
        if seq not in seq_to_input:
            seq_to_input[seq] = PeptideInput(seq)

    if not seq_to_input:
        return {key: 0.0 for key, *_ in junction_contexts}, 'pepsickle_unassessed'

    unique_inputs = list(seq_to_input.values())
    profiles = audit_pepsickle_inputs(
        unique_inputs,
        human_only=policy.pepsickle_human_only,
        threshold=policy.pepsickle_threshold,
    )
    by_seq = {p.peptide.sequence: p for p in profiles}
    any_observed = any(p.status == 'observations_returned' for p in profiles)
    if not any_observed:
        reasons = sorted({r for p in profiles for r in p.reason_codes})
        logger.warning(
            "Pepsickle cleavage scoring produced no observations (%s); "
            "continuing with zero cleavage penalty.",
            ', '.join(reasons) or 'unassessed',
        )
        return {key: 0.0 for key, *_ in junction_contexts}, 'pepsickle_unassessed'

    penalties = {}
    for key, (seq, bond) in key_to_seq_and_bond.items():
        if seq is None:
            penalties[key] = 0.0
            continue
        profile = by_seq.get(seq)
        if profile is None or profile.status != 'observations_returned':
            penalties[key] = 0.0
            continue
        site = next((s for s in profile.sites if s.bond == bond), None)
        if site is None or site.score is None:
            penalties[key] = 0.0
        else:
            penalties[key] = max(0.0, min(1.0, 1.0 - float(site.score)))
    return penalties, 'pepsickle_applied'


def _path_cost(order, start_cost, edge_cost, end_cost):
    if not order:
        return 0.0
    total = start_cost[order[0]] + end_cost[order[-1]]
    for idx in range(len(order) - 1):
        total += edge_cost[order[idx]][order[idx + 1]]
    return total


def _solve_held_karp(n, start_cost, edge_cost, end_cost):
    """Exact Hamiltonian path optimization via Held-Karp bitmask DP.

    Minimizes ``start_cost[p_0] + sum(edge_cost[p_k][p_{k+1}]) + end_cost[p_{n-1}]``
    across all permutations of ``range(n)``. Exact ties prefer the original
    order ``[0, 1, ..., n-1]`` (and lexicographically smaller permutations).
    """
    if n <= 1:
        return list(range(n))

    eps = 1e-12
    full_mask = (1 << n) - 1
    inf = float('inf')

    # dp[mask][last] = minimal cost to visit `mask` ending at `last`
    dp = [[inf] * n for _ in range(1 << n)]
    parent = [[-1] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = start_cost[i]

    for mask in range(1, 1 << n):
        # Skip single-bit masks (already initialized)
        if mask & (mask - 1) == 0:
            continue
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            prev_mask = mask ^ (1 << last)
            best_val = inf
            best_prev = -1
            for prev in range(n):
                if not (prev_mask & (1 << prev)):
                    continue
                cand = dp[prev_mask][prev] + edge_cost[prev][last]
                if cand < best_val - eps:
                    best_val = cand
                    best_prev = prev
            dp[mask][last] = best_val
            parent[mask][last] = best_prev

    # Prefer the identity order [0, 1, ..., n-1] on exact ties.
    identity_order = list(range(n))
    identity_cost = _path_cost(identity_order, start_cost, edge_cost, end_cost)

    best_total = inf
    best_last = 0
    for last in range(n):
        total = dp[full_mask][last] + end_cost[last]
        if total < best_total - eps:
            best_total = total
            best_last = last

    if identity_cost <= best_total + eps:
        return identity_order

    order = []
    mask = full_mask
    curr = best_last
    while curr != -1:
        order.append(curr)
        prev = parent[mask][curr]
        mask ^= (1 << curr)
        curr = prev
    order.reverse()
    return order


def _solve_two_opt(n, start_cost, edge_cost, end_cost, max_passes=50):
    """Deterministic local search (greedy multi-start + 2-opt + Or-opt-1) for
    directed Hamiltonian paths on large constructs (``n > exact_dp_max_antigens``).
    """
    if n <= 1:
        return list(range(n))

    eps = 1e-12
    best_order = list(range(n))
    best_cost = _path_cost(best_order, start_cost, edge_cost, end_cost)

    # Multi-start greedy nearest-neighbor initialization
    for start in range(n):
        remaining = set(range(n))
        remaining.remove(start)
        tour = [start]
        while remaining:
            last = tour[-1]
            nxt = min(remaining, key=lambda node: (edge_cost[last][node], node))
            tour.append(nxt)
            remaining.remove(nxt)
        c = _path_cost(tour, start_cost, edge_cost, end_cost)
        if c < best_cost - eps:
            best_cost = c
            best_order = tour

    # Steepest-descent 2-opt reversals + Or-opt single-node relocations
    order = list(best_order)
    current_cost = best_cost
    for _ in range(max_passes):
        improved = False
        best_cand_order = None
        best_cand_cost = current_cost - eps

        # 2-opt segment reversals
        for i in range(n - 1):
            for j in range(i + 1, n):
                cand = order[:i] + order[i:j + 1][::-1] + order[j + 1:]
                c = _path_cost(cand, start_cost, edge_cost, end_cost)
                if c < best_cand_cost:
                    best_cand_cost = c
                    best_cand_order = cand
                    improved = True

        # Or-opt-1 single antigen relocation (effective on asymmetric directed costs)
        for i in range(n):
            node = order[i]
            without = order[:i] + order[i + 1:]
            for j in range(n):
                if j == i:
                    continue
                cand = without[:j] + [node] + without[j:]
                c = _path_cost(cand, start_cost, edge_cost, end_cost)
                if c < best_cand_cost:
                    best_cand_cost = c
                    best_cand_order = cand
                    improved = True

        if not improved or best_cand_order is None:
            break
        order = best_cand_order
        current_cost = best_cand_cost

    identity_cost = _path_cost(list(range(n)), start_cost, edge_cost, end_cost)
    if identity_cost <= current_cost + eps:
        return list(range(n))
    return order


def optimize_antigen_order_and_linkers(
        antigen_aas, alleles, predictor,
        candidate_names=JUNCTION_SWAP_CANDIDATES,
        k_lengths=(8, 9, 10, 11),
        reference_proteome=None,
        rank_strong=RANK_STRONG, rank_mild=RANK_MILD,
        rank_kind='auto',
        signal_peptide_aa=None,
        mitd_aa=None,
        default_linker_name=None,
        policy=None,
        cleavage_scorer=None):
    """Jointly optimize mRNA antigen order and per-junction linkers.

    Precomputes the directed ``N x N x M`` inter-antigen junction cost matrix
    (plus ``signal_peptide -> first_antigen`` and ``last_antigen -> linker -> MITD``
    boundaries when present and enabled), batching all unique chimeric k-mers
    into a single MHC predictor call. Uses exact Held-Karp dynamic programming
    for ``N <= policy.exact_dp_max_antigens`` and deterministic 2-opt/Or-opt
    local search for larger poly-epitope strings.

    Parameters
    ----------
    antigen_aas : list[str]
        Candidate antigen amino-acid sequences in initial ranked order.
    alleles : list[str]
        Patient MHC alleles (any species supported by mhcgnomes/predictor).
    predictor : mhctools.BasePredictor or topiary.TopiaryPredictor
        MHC binding/presentation predictor covering ``alleles``.
    candidate_names : iterable of str
        Candidate linker names resolvable by ``vaccine_library.get_linker``.
    k_lengths : iterable of int
        Chimeric k-mer lengths to evaluate (default 8–11).
    reference_proteome : optional container
        Self k-mer set excluded from chimeric burden.
    rank_strong, rank_mild : float
        Percentile-rank thresholds (0 <= strong <= mild <= 100).
    rank_kind : str
        'auto', 'pMHC_affinity', or 'pMHC_presentation'.
    signal_peptide_aa : optional str
        N-terminal signal peptide AA sequence (boundary scored without linker).
    mitd_aa : optional str
        C-terminal MITD AA sequence (boundary scored across candidate linkers).
    default_linker_name : optional str
        Shared default linker name used to compute the fixed-order baseline.
    policy : dict or MRNAJunctionOptimizationPolicy or None
        Per-trial weights and solver configuration.
    cleavage_scorer : optional callable
        Optional ``(left_aa, linker_aa, right_aa) -> float`` returning a
        cleavage penalty in ``[0.0, 1.0]`` (0 = ideal cleavage, 1 = poor).
        Overrides ``policy.cleavage_scoring='pepsickle'`` when provided.

    Returns
    -------
    JunctionSwapResult
    """
    alleles = tuple(alleles or ())
    k_lengths = tuple(k_lengths)
    policy = MRNAJunctionOptimizationPolicy.resolve(policy)
    if (not all(isinstance(value, Real) and not isinstance(value, bool)
                and math.isfinite(value) for value in (rank_strong, rank_mild))
            or not 0 <= rank_strong <= rank_mild <= 100):
        raise ValueError('Junction thresholds require 0 <= strong <= mild <= 100')

    n = len(antigen_aas)
    if n == 0:
        return JunctionSwapResult(
            chosen_linker_per_junction=[],
            burden=0, strong_burden=0, chimeric_predictions=[],
            antigen_order=[], reordered=False, solver='held_karp_dp',
            policy=policy.to_dict(),
        )

    candidates = [get_linker(name) for name in candidate_names]
    default_idx = None
    if default_linker_name is not None:
        target = default_linker_name.upper()
        for idx, c in enumerate(candidates):
            if c.name.upper() == target:
                default_idx = idx
                break
        if default_idx is None:
            candidates.append(get_linker(default_linker_name))
            default_idx = len(candidates) - 1

    use_sp = bool(signal_peptide_aa and policy.include_signal_peptide_boundary)
    use_mitd = bool(mitd_aa and policy.include_mitd_boundary)

    # 1. Enumerate all chimeric k-mers across the entire N x N x M matrix + boundaries
    sp_kmers = {}
    if use_sp:
        for i in range(n):
            sp_kmers[i] = boundary_kmers(
                signal_peptide_aa, antigen_aas[i], k_lengths,
                reference_proteome=reference_proteome)

    pair_kmers = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for c_idx, cand in enumerate(candidates):
                pair_kmers[(i, j, c_idx)] = junction_kmers(
                    antigen_aas[i], cand.amino_acids, antigen_aas[j],
                    k_lengths, reference_proteome=reference_proteome)

    mitd_kmers = {}
    if use_mitd:
        for i in range(n):
            for c_idx, cand in enumerate(candidates):
                mitd_kmers[(i, c_idx)] = junction_kmers(
                    antigen_aas[i], cand.amino_acids, mitd_aa,
                    k_lengths, reference_proteome=reference_proteome)

    # Deduplicated single batch MHC prediction call across all candidate junctions
    unique_kmers = []
    seen_kmers = set()
    for kmer_list in (*sp_kmers.values(), *pair_kmers.values(), *mitd_kmers.values()):
        for kmer in kmer_list:
            if kmer not in seen_kmers:
                seen_kmers.add(kmer)
                unique_kmers.append(kmer)

    from topiary import TopiaryPredictor, describe_prediction_models
    configured = (predictor if isinstance(predictor, TopiaryPredictor)
                  else TopiaryPredictor(models=[predictor]))
    models_meta = []
    all_rows = score_junction_kmers(
        unique_kmers, alleles, predictor, models=models_meta, rank_kind=rank_kind)
    rows_by_kmer = {}
    for row in all_rows:
        rows_by_kmer.setdefault(row[0], []).append(row)

    def _rows_for(kmer_list):
        out = []
        for km in kmer_list:
            out.extend(rows_by_kmer.get(km, ()))
        return out

    # 2. Optional proteasomal cleavage scoring
    cleavage_penalties = {}
    cleavage_status = 'disabled'
    if cleavage_scorer is not None:
        cleavage_status = 'custom_scorer'
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                for c_idx, cand in enumerate(candidates):
                    val = float(cleavage_scorer(
                        antigen_aas[i], cand.amino_acids, antigen_aas[j]))
                    if not math.isfinite(val) or val < 0.0:
                        raise ValueError(
                            f"cleavage_scorer returned invalid penalty {val!r}")
                    cleavage_penalties[('pair', i, j, c_idx)] = val
            if use_mitd:
                for c_idx, cand in enumerate(candidates):
                    val = float(cleavage_scorer(
                        antigen_aas[i], cand.amino_acids, mitd_aa))
                    if not math.isfinite(val) or val < 0.0:
                        raise ValueError(
                            f"cleavage_scorer returned invalid penalty {val!r}")
                    cleavage_penalties[('mitd', i, c_idx)] = val
    elif policy.cleavage_scoring == 'pepsickle' and policy.cleavage_penalty_weight > 0:
        contexts = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                for c_idx, cand in enumerate(candidates):
                    contexts.append((
                        ('pair', i, j, c_idx),
                        antigen_aas[i], cand.amino_acids, antigen_aas[j]))
            if use_mitd:
                for c_idx, cand in enumerate(candidates):
                    contexts.append((
                        ('mitd', i, c_idx),
                        antigen_aas[i], cand.amino_acids, mitd_aa))
        cleavage_penalties, cleavage_status = _compute_pepsickle_boundary_penalties(
            contexts, policy)

    # 3. Precompute boundary costs and best linker per directed edge
    start_cost = [0.0] * n
    start_stats = [(0, 0, 100.0, []) for _ in range(n)]
    if use_sp:
        for i in range(n):
            rows = _rows_for(sp_kmers[i])
            cost, strong, mild, worst = _edge_composite_cost(
                rows, '', 0.0, policy, rank_strong, rank_mild)
            start_cost[i] = cost
            start_stats[i] = (strong, mild, worst, rows)

    edge_cost = [[float('inf')] * n for _ in range(n)]
    best_pair_choice = {}  # (i, j) -> (cost, c_idx, strong, mild, worst, cleavage_pen, rows)
    pair_cand_stats = {}   # (i, j, c_idx) -> (cost, strong, mild, worst, cleavage_pen, rows)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            best = None
            for c_idx, cand in enumerate(candidates):
                rows = _rows_for(pair_kmers[(i, j, c_idx)])
                c_pen = cleavage_penalties.get(('pair', i, j, c_idx), 0.0)
                cost, strong, mild, worst = _edge_composite_cost(
                    rows, cand.amino_acids, c_pen, policy, rank_strong, rank_mild)
                pair_cand_stats[(i, j, c_idx)] = (cost, strong, mild, worst, c_pen, rows)
                tie_key = (cost, strong, mild, -worst, c_idx)
                if best is None or tie_key < best[0]:
                    best = (tie_key, (cost, c_idx, strong, mild, worst, c_pen, rows))
            edge_cost[i][j] = best[1][0]
            best_pair_choice[(i, j)] = best[1]

    end_cost = [0.0] * n
    best_mitd_choice = {}  # i -> (cost, c_idx, strong, mild, worst, cleavage_pen, rows)
    mitd_cand_stats = {}   # (i, c_idx) -> (cost, strong, mild, worst, cleavage_pen, rows)
    if use_mitd:
        for i in range(n):
            best = None
            for c_idx, cand in enumerate(candidates):
                rows = _rows_for(mitd_kmers[(i, c_idx)])
                c_pen = cleavage_penalties.get(('mitd', i, c_idx), 0.0)
                cost, strong, mild, worst = _edge_composite_cost(
                    rows, cand.amino_acids, c_pen, policy, rank_strong, rank_mild)
                mitd_cand_stats[(i, c_idx)] = (cost, strong, mild, worst, c_pen, rows)
                tie_key = (cost, strong, mild, -worst, c_idx)
                if best is None or tie_key < best[0]:
                    best = (tie_key, (cost, c_idx, strong, mild, worst, c_pen, rows))
            end_cost[i] = best[1][0]
            best_mitd_choice[i] = best[1]

    # 4. Solve for optimal antigen order
    if n <= policy.exact_dp_max_antigens:
        solver = 'held_karp_dp'
        order = _solve_held_karp(n, start_cost, edge_cost, end_cost)
    else:
        solver = 'two_opt'
        order = _solve_two_opt(
            n, start_cost, edge_cost, end_cost, max_passes=policy.max_2opt_passes)

    # 5. Evaluate chosen solution, fixed-order optimal linkers, and fixed-order default linker
    chosen_linkers = []
    all_predictions = []
    junction_decisions = []
    requests_meta = []

    sp_strong, sp_mild, _, sp_rows = start_stats[order[0]]
    total_strong = sp_strong
    total_mild = sp_mild
    total_cleavage_pen = 0.0
    total_cost = start_cost[order[0]]
    scored_count = len(all_rows)

    for kmer, allele, rank in sp_rows:
        all_predictions.append(JunctionPrediction(
            junction_index=-1, linker_name='signal_peptide',
            kmer=kmer, allele=allele, rank=rank))

    for pos in range(n - 1):
        left_i, right_j = order[pos], order[pos + 1]
        cost, c_idx, strong, mild, worst, c_pen, rows = best_pair_choice[(left_i, right_j)]
        cand = candidates[c_idx]
        chosen_linkers.append(cand)
        total_strong += strong
        total_mild += mild
        total_cleavage_pen += c_pen
        total_cost += cost
        for kmer, allele, rank in rows:
            all_predictions.append(JunctionPrediction(
                junction_index=pos, linker_name=cand.name,
                kmer=kmer, allele=allele, rank=rank))
        per_cand_entries = []
        for idx, c in enumerate(candidates):
            c_cost, c_strong, c_mild, c_worst, _, _ = pair_cand_stats[(left_i, right_j, idx)]
            per_cand_entries.append(dict(
                linker=c.name,
                key=[c_strong, c_mild, -c_worst],
                composite_cost=round(c_cost, 6),
                candidate_order=idx,
            ))
            km = pair_kmers[(left_i, right_j, idx)]
            if km:
                requests_meta.append({
                    'junction_index': pos, 'linker_name': c.name, 'peptides': list(km)})
        junction_decisions.append(dict(
            junction_index=pos, preferred_linker=cand.name,
            ordering='composite_cost (strong, mild, cleavage, length, worst_rank), candidate_order',
            candidates=per_cand_entries,
        ))

    if use_mitd:
        last_i = order[-1]
        cost, c_idx, strong, mild, worst, c_pen, rows = best_mitd_choice[last_i]
        cand = candidates[c_idx]
        chosen_linkers.append(cand)
        total_strong += strong
        total_mild += mild
        total_cleavage_pen += c_pen
        total_cost += cost
        mitd_pos = n - 1
        for kmer, allele, rank in rows:
            all_predictions.append(JunctionPrediction(
                junction_index=mitd_pos, linker_name=cand.name,
                kmer=kmer, allele=allele, rank=rank))
        per_cand_entries = []
        for idx, c in enumerate(candidates):
            c_cost, c_strong, c_mild, c_worst, _, _ = mitd_cand_stats[(last_i, idx)]
            per_cand_entries.append(dict(
                linker=c.name,
                key=[c_strong, c_mild, -c_worst],
                composite_cost=round(c_cost, 6),
                candidate_order=idx,
            ))
            km = mitd_kmers[(last_i, idx)]
            if km:
                requests_meta.append({
                    'junction_index': mitd_pos, 'linker_name': c.name, 'peptides': list(km)})
        junction_decisions.append(dict(
            junction_index=mitd_pos, preferred_linker=cand.name,
            ordering='composite_cost (strong, mild, cleavage, length, worst_rank), candidate_order',
            candidates=per_cand_entries,
        ))

    # Fixed-order optimal linkers baseline (original order [0..n-1])
    fixed_strong = start_stats[0][0]
    fixed_mild = start_stats[0][1]
    for pos in range(n - 1):
        _, _, s, m, _, _, _ = best_pair_choice[(pos, pos + 1)]
        fixed_strong += s
        fixed_mild += m
    if use_mitd:
        _, _, s, m, _, _, _ = best_mitd_choice[n - 1]
        fixed_strong += s
        fixed_mild += m

    # Fixed-order shared default linker baseline
    default_strong = start_stats[0][0]
    default_mild = start_stats[0][1]
    default_cleavage_pen = 0.0
    default_cost = start_cost[0]
    if default_idx is not None:
        for pos in range(n - 1):
            c_cost, s, m, _, c_pen, _ = pair_cand_stats[(pos, pos + 1, default_idx)]
            default_strong += s
            default_mild += m
            default_cleavage_pen += c_pen
            default_cost += c_cost
        if use_mitd:
            c_cost, s, m, _, c_pen, _ = mitd_cand_stats[(n - 1, default_idx)]
            default_strong += s
            default_mild += m
            default_cleavage_pen += c_pen
            default_cost += c_cost

    prediction_metadata = {
        'configured_models': describe_prediction_models(configured),
        'predictor_class': type(predictor).__module__ + '.' + type(predictor).__name__,
        'alleles': list(alleles),
        'peptide_lengths': list(k_lengths),
        'rank_kind_policy': ('affinity_if_available_else_presentation'
                             if rank_kind == 'auto' else 'explicit'),
        'rank_kind': junction_rank_kind(predictor, rank_kind),
        'rank_metric': 'percentile_rank',
        'rank_strong': rank_strong,
        'rank_mild': rank_mild,
        'context': 'peptide_only',
        'coverage': 'complete',
        'models': models_meta,
        'requests': requests_meta,
        'scored_prediction_count': scored_count,
        'unique_kmers_evaluated': len(unique_kmers),
        'solver': solver,
        'cleavage_status': cleavage_status,
        'junction_decisions': junction_decisions,
    }

    result = JunctionSwapResult(
        chosen_linker_per_junction=chosen_linkers,
        burden=total_mild,
        strong_burden=total_strong,
        chimeric_predictions=all_predictions,
        prediction_metadata=prediction_metadata,
        antigen_order=order,
        reordered=(order != list(range(n))),
        solver=solver,
        fixed_order_burden=fixed_mild,
        fixed_order_strong_burden=fixed_strong,
        signal_peptide_burden=sp_mild,
        signal_peptide_strong_burden=sp_strong,
        cleavage_penalty=round(total_cleavage_pen, 6),
        default_cleavage_penalty=round(default_cleavage_pen, 6),
        cleavage_status=cleavage_status,
        total_cost=round(total_cost, 6),
        default_total_cost=(round(default_cost, 6) if default_idx is not None else None),
        policy=policy.to_dict(),
    )
    if default_idx is not None:
        result.default_burden = default_mild
        result.default_strong_burden = default_strong
    return result
