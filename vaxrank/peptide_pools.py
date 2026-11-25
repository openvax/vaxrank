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

"""Multi-pool synthetic long peptide (SLP) formulation partitioning.

Clinical personalized peptide vaccines (e.g. PGV-001, NEO-PV-01, NeoVax)
commonly split 16–20 emitted peptides across 2–4 injection pools delivered at
separate anatomical sites so:

1. Peptides competing for the **same patient MHC allele** are placed in
   different pools, reducing MHC-loading competition and immunodominance in
   the same draining lymph node.
2. **Cysteine-containing peptides** are spread across pools to reduce
   inter-peptide disulfide cross-linking risk during formulation.
3. **Hydrophobicity (max 7-mer GRAVY)** and **net charge at formulation pH**
   are balanced across pools to reduce co-precipitation and solubility
   failures.

Every tier cutoff and objective weight lives on :class:`PeptidePoolPolicy`
so clinical trials and non-human protocols can persist and version their
formulation policy in YAML (``peptide.pool_partition:``). The optimizer is
deterministic, species-agnostic (``mhcgnomes``-backed canonical MHC alleles),
enforces balanced pool cardinalities (``floor(N/K)..ceil(N/K)``), and leaves
single-pool (``num_pools == 1``) designs untouched.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, fields
import math
from typing import Iterable, Optional

from .identifiers import normalize_mhc_allele


# Lehninger Principles of Biochemistry standard pKa values for
# Henderson-Hasselbalch peptide net-charge calculation.
_PKA_N_TERM = 9.69
_PKA_C_TERM = 2.34
_PKA_POSITIVE = {'K': 10.53, 'R': 12.48, 'H': 6.00}
_PKA_NEGATIVE = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10.07}


def _check_finite_nonnegative(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f"PeptidePoolPolicy.{name} must be finite and >= 0; got {value!r}")
    return float(value)


@dataclass(frozen=True)
class PeptidePoolPolicy:
    """Explicit, serializable per-trial policy for multi-pool SLP partitioning.

    All cutoffs and objective weights can be saved in a trial YAML config
    under ``peptide.pool_partition:`` and are recorded in ``selection_policy.json``,
    ``manifest.json``, and reports.
    """

    strong_rank_threshold: float = 0.5
    medium_rank_threshold: float = 2.0
    strong_score_threshold: float = 0.5
    medium_score_threshold: float = 0.1
    strong_allele_weight: float = 3.0
    medium_allele_weight: float = 1.0
    allele_competition_weight: float = 100.0
    cysteine_pair_weight: float = 10.0
    gravy_balance_weight: float = 2.0
    charge_balance_weight: float = 1.0
    charge_ph: float = 7.0
    max_swap_passes: int = 50

    def __post_init__(self):
        sr = _check_finite_nonnegative("strong_rank_threshold", self.strong_rank_threshold)
        mr = _check_finite_nonnegative("medium_rank_threshold", self.medium_rank_threshold)
        if sr <= 0:
            raise ValueError(f"PeptidePoolPolicy.strong_rank_threshold must be > 0; got {sr!r}")
        if mr < sr:
            raise ValueError(
                f"PeptidePoolPolicy.medium_rank_threshold ({mr}) must be >= strong_rank_threshold ({sr})"
            )
        ss = _check_finite_nonnegative("strong_score_threshold", self.strong_score_threshold)
        ms = _check_finite_nonnegative("medium_score_threshold", self.medium_score_threshold)
        if ss < ms:
            raise ValueError(
                f"PeptidePoolPolicy.strong_score_threshold ({ss}) must be >= medium_score_threshold ({ms})"
            )
        for name in (
            "strong_allele_weight",
            "medium_allele_weight",
            "allele_competition_weight",
            "cysteine_pair_weight",
            "gravy_balance_weight",
            "charge_balance_weight",
        ):
            _check_finite_nonnegative(name, getattr(self, name))
        ph = _check_finite_nonnegative("charge_ph", self.charge_ph)
        if ph > 14.0:
            raise ValueError(f"PeptidePoolPolicy.charge_ph must be in [0, 14]; got {ph!r}")
        if isinstance(self.max_swap_passes, bool) or not isinstance(self.max_swap_passes, int) or self.max_swap_passes < 0:
            raise ValueError(
                f"PeptidePoolPolicy.max_swap_passes must be a nonnegative integer; got {self.max_swap_passes!r}"
            )

    @classmethod
    def resolve(cls, value=None) -> "PeptidePoolPolicy":
        """Resolve ``None``, a ``PeptidePoolPolicy``, or a dict to a validated policy."""
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise TypeError(
                f"peptide.pool_partition must be a mapping or PeptidePoolPolicy; got {type(value).__name__}"
            )
        valid = {f.name for f in fields(cls)}
        unknown = sorted(set(value) - valid)
        if unknown:
            raise ValueError(
                f"Unknown peptide.pool_partition keys: {', '.join(unknown)}"
            )
        cleaned = {k: v for k, v in value.items() if v is not None}
        return cls(**cleaned)

    def to_dict(self) -> dict:
        return asdict(self)


def net_charge_at_ph(
    sequence: str,
    *,
    ph: float = 7.0,
    n_terminal_acetylation: bool = False,
    c_terminal_amidation: bool = False,
) -> float:
    """Estimate peptide net charge at ``ph`` via Henderson-Hasselbalch pKa."""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    charge = 0.0
    for aa in seq:
        pka_pos = _PKA_POSITIVE.get(aa)
        if pka_pos is not None:
            charge += 1.0 / (1.0 + 10.0 ** (ph - pka_pos))
        else:
            pka_neg = _PKA_NEGATIVE.get(aa)
            if pka_neg is not None:
                charge -= 1.0 / (1.0 + 10.0 ** (pka_neg - ph))
    if not n_terminal_acetylation:
        charge += 1.0 / (1.0 + 10.0 ** (ph - _PKA_N_TERM))
    if not c_terminal_amidation:
        charge -= 1.0 / (1.0 + 10.0 ** (_PKA_C_TERM - ph))
    rounded = round(charge, 2)
    return 0.0 if rounded == 0.0 else rounded


def net_charge_ph7(
    sequence: str,
    *,
    n_terminal_acetylation: bool = False,
    c_terminal_amidation: bool = False,
) -> float:
    """Estimate peptide net charge at physiological pH (7.0)."""
    return net_charge_at_ph(
        sequence,
        ph=7.0,
        n_terminal_acetylation=n_terminal_acetylation,
        c_terminal_amidation=c_terminal_amidation,
    )


@dataclass(frozen=True)
class PeptidePoolProfile:
    """Immunological and physicochemical profile of one emitted peptide."""

    index: int
    name: str
    allele_weights: tuple[tuple[str, float], ...]
    strong_alleles: tuple[str, ...]
    medium_alleles: tuple[str, ...]
    max_gravy: float
    cterm_gravy: float
    has_cysteine: bool
    net_charge: float

    @property
    def weight_map(self) -> dict[str, float]:
        return dict(self.allele_weights)

    @property
    def all_alleles(self) -> tuple[str, ...]:
        return tuple(sorted(set(self.strong_alleles) | set(self.medium_alleles)))


def _canonical_allele(allele: str) -> str:
    try:
        return normalize_mhc_allele(allele)
    except Exception:
        return str(allele).strip()


def _extract_construct_epitopes(construct, vp_by_sequence: dict[str, list]):
    """Collect recorded target epitopes belonging to an emitted peptide construct."""
    matched = []
    seen_ids = set()

    def _add(vp):
        if id(vp) not in seen_ids:
            seen_ids.add(id(vp))
            matched.append(vp)

    for vp in vp_by_sequence.get(construct.sequence, ()):
        _add(vp)
    for seg in construct.antigens or ():
        for vp in vp_by_sequence.get(seg.get('aa', ''), ()):
            _add(vp)
        for src in seg.get('sources') or ():
            if isinstance(src, dict):
                for vp in vp_by_sequence.get(src.get('candidate_sequence', ''), ()):
                    _add(vp)

    epitopes = []
    for vp in matched:
        for ep in getattr(vp, 'epitopes', ()) or ():
            if getattr(ep, 'is_targetable', True):
                epitopes.append(ep)
    return epitopes


def _allele_weights_for_epitopes(
    epitopes,
    target_alleles: Optional[Iterable[str]] = None,
    *,
    policy: Optional[PeptidePoolPolicy] = None,
) -> tuple[dict[str, float], tuple[str, ...], tuple[str, ...]]:
    """Map each bound MHC allele to a competition weight under ``policy``.

    Evaluates best percentile rank (across ``pMHC_presentation`` and
    ``pMHC_affinity`` predictions) against ``policy.strong_rank_threshold`` and
    ``policy.medium_rank_threshold``. Falls back to ``per_allele_scores``
    against ``policy.strong_score_threshold`` / ``policy.medium_score_threshold``
    when percentile ranks are absent.
    """
    policy = PeptidePoolPolicy.resolve(policy)
    allowed = None
    if target_alleles:
        allowed = {_canonical_allele(a) for a in target_alleles if a}
        allowed.discard('')

    if not epitopes:
        return {}, (), ()

    best_pct: dict[str, float] = {}
    for ep in epitopes:
        preds = ep.predictions_flat() if hasattr(ep, 'predictions_flat') else ()
        for pred in preds:
            if pred.kind not in ('pMHC_affinity', 'pMHC_presentation'):
                continue
            pct = getattr(pred, 'percentile_rank', None)
            if pct is None or not pred.allele:
                continue
            try:
                pct_f = float(pct)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(pct_f):
                continue
            canon = _canonical_allele(pred.allele)
            if not canon or (allowed is not None and canon not in allowed):
                continue
            prev = best_pct.get(canon)
            if prev is None or pct_f < prev:
                best_pct[canon] = pct_f

    strong = set()
    medium = set()
    if best_pct:
        for allele, pct in best_pct.items():
            if pct <= policy.strong_rank_threshold:
                strong.add(allele)
            elif pct <= policy.medium_rank_threshold:
                medium.add(allele)
    else:
        best_score: dict[str, float] = {}
        for ep in epitopes:
            for allele, score in (getattr(ep, 'per_allele_scores', {}) or {}).items():
                if score is None:
                    continue
                canon = _canonical_allele(allele)
                if not canon or (allowed is not None and canon not in allowed):
                    continue
                score_f = float(score)
                prev = best_score.get(canon, 0.0)
                if score_f > prev:
                    best_score[canon] = score_f
        for allele, score in best_score.items():
            if score >= policy.strong_score_threshold:
                strong.add(allele)
            elif score >= policy.medium_score_threshold:
                medium.add(allele)

    medium -= strong
    weights = {a: policy.strong_allele_weight for a in strong}
    weights.update({a: policy.medium_allele_weight for a in medium})
    return weights, tuple(sorted(strong)), tuple(sorted(medium))


def _build_profiles(
    constructs,
    ranked_vaccine_peptides=(),
    target_alleles=None,
    *,
    policy: Optional[PeptidePoolPolicy] = None,
) -> list[PeptidePoolProfile]:
    policy = PeptidePoolPolicy.resolve(policy)
    vp_by_sequence: dict[str, list] = {}
    for _, vps in ranked_vaccine_peptides or ():
        for vp in vps:
            vp_by_sequence.setdefault(vp.amino_acids, []).append(vp)

    profiles = []
    for idx, construct in enumerate(constructs):
        epitopes = _extract_construct_epitopes(construct, vp_by_sequence)
        weights, strong, medium = _allele_weights_for_epitopes(
            epitopes, target_alleles, policy=policy
        )
        mfg = construct.manufacturability or {}
        max_gravy = float(mfg.get('max_7mer_gravy_score', 0.0) or 0.0)
        cterm_gravy = float(mfg.get('cterm_7mer_gravy_score', 0.0) or 0.0)
        cys_count = int(mfg.get('cysteine_count', 0) or 0)
        charge = net_charge_at_ph(
            construct.sequence,
            ph=policy.charge_ph,
            n_terminal_acetylation=bool(construct.components.get('n_terminal_acetylation', False)),
            c_terminal_amidation=bool(construct.components.get('c_terminal_amidation', False)),
        )
        profiles.append(PeptidePoolProfile(
            index=idx,
            name=construct.name,
            allele_weights=tuple(sorted(weights.items())),
            strong_alleles=strong,
            medium_alleles=medium,
            max_gravy=max_gravy,
            cterm_gravy=cterm_gravy,
            has_cysteine=cys_count > 0,
            net_charge=charge,
        ))
    return profiles


def _pool_cost(
    pool_profiles: list[PeptidePoolProfile],
    mean_gravy: float,
    mean_charge: float,
    gravy_scale: float = 1.0,
    charge_scale: float = 1.0,
    *,
    policy: Optional[PeptidePoolPolicy] = None,
) -> float:
    """Compute objective penalty for one pool under ``policy``."""
    if not pool_profiles:
        return 0.0
    policy = policy or PeptidePoolPolicy()

    # Pairwise same-allele weighted competition: sum_{i<j} w_{i,a} * w_{j,a}
    # = 0.5 * ((sum w)^2 - sum(w^2)) per allele.
    sum_w: dict[str, float] = {}
    sum_w2: dict[str, float] = {}
    cys = 0
    g_sum = 0.0
    q_sum = 0.0
    for prof in pool_profiles:
        for allele, w in prof.allele_weights:
            sum_w[allele] = sum_w.get(allele, 0.0) + w
            sum_w2[allele] = sum_w2.get(allele, 0.0) + w * w
        if prof.has_cysteine:
            cys += 1
        g_sum += prof.max_gravy
        q_sum += prof.net_charge

    allele_pen = 0.5 * sum(
        sum_w[a] * sum_w[a] - sum_w2[a]
        for a in sum_w
    )
    cys_pen = cys * (cys - 1) / 2.0
    n = len(pool_profiles)
    gravy_dev = ((g_sum / n - mean_gravy) / gravy_scale) ** 2
    charge_dev = ((q_sum / n - mean_charge) / charge_scale) ** 2
    return (
        policy.allele_competition_weight * allele_pen
        + policy.cysteine_pair_weight * cys_pen
        + policy.gravy_balance_weight * gravy_dev
        + policy.charge_balance_weight * charge_dev
    )


def partition_peptide_constructs(
    constructs,
    num_pools: int,
    *,
    ranked_vaccine_peptides=(),
    target_alleles=None,
    policy=None,
) -> Optional[dict]:
    """Assign emitted peptide constructs across ``1..num_pools`` formulation pools.

    Mutates each construct's ``components`` dict in place to add ``'pool'``
    (1-based integer) and ``'pool_profile'`` (bound MHC alleles, net charge,
    formulation pH, and policy) when ``num_pools > 1``. Returns a
    JSON-serializable pool summary dict, or ``None`` when ``num_pools <= 1``
    or ``constructs`` is empty.
    """
    if isinstance(num_pools, bool) or not isinstance(num_pools, int) or num_pools < 1:
        raise ValueError("num_pools must be a positive integer; got %r" % (num_pools,))
    resolved_policy = PeptidePoolPolicy.resolve(policy)
    if num_pools == 1 or not constructs:
        return None

    k = min(num_pools, len(constructs))
    profiles = _build_profiles(
        constructs,
        ranked_vaccine_peptides=ranked_vaccine_peptides,
        target_alleles=target_alleles,
        policy=resolved_policy,
    )
    n = len(profiles)
    mean_gravy = sum(p.max_gravy for p in profiles) / n
    mean_charge = sum(p.net_charge for p in profiles) / n
    gravy_scale = max(
        math.sqrt(sum((p.max_gravy - mean_gravy) ** 2 for p in profiles) / n),
        1.0,
    )
    charge_scale = max(
        math.sqrt(sum((p.net_charge - mean_charge) ** 2 for p in profiles) / n),
        1.0,
    )

    # Balanced pool cardinality caps: first (n % k) pools hold ceil(n / k),
    # remaining pools hold floor(n / k).
    base_size, extra = divmod(n, k)
    caps = [base_size + (1 if i < extra else 0) for i in range(k)]

    # Greedy initialization: place most-constrained peptides first
    # (most strong alleles -> most total alleles -> cysteine -> high GRAVY -> rank).
    order = sorted(
        profiles,
        key=lambda p: (
            -len(p.strong_alleles),
            -len(p.all_alleles),
            -int(p.has_cysteine),
            -abs(p.max_gravy),
            p.index,
        ),
    )
    pools: list[list[PeptidePoolProfile]] = [[] for _ in range(k)]
    for prof in order:
        best_p = None
        best_delta = math.inf
        best_size = math.inf
        for pi in range(k):
            if len(pools[pi]) >= caps[pi]:
                continue
            before = _pool_cost(
                pools[pi], mean_gravy, mean_charge, gravy_scale, charge_scale,
                policy=resolved_policy,
            )
            after = _pool_cost(
                pools[pi] + [prof], mean_gravy, mean_charge, gravy_scale, charge_scale,
                policy=resolved_policy,
            )
            delta = after - before
            if (
                delta < best_delta - 1e-12
                or (abs(delta - best_delta) <= 1e-12 and len(pools[pi]) < best_size)
            ):
                best_delta = delta
                best_size = len(pools[pi])
                best_p = pi
        pools[best_p].append(prof)

    # Deterministic steepest-descent 2-swap local search preserving exact
    # balanced pool sizes.
    for _ in range(resolved_policy.max_swap_passes):
        improved = False
        best_swap = None
        best_gain = 1e-12
        for p1 in range(k):
            for p2 in range(p1 + 1, k):
                c12_before = (
                    _pool_cost(
                        pools[p1], mean_gravy, mean_charge, gravy_scale, charge_scale,
                        policy=resolved_policy,
                    )
                    + _pool_cost(
                        pools[p2], mean_gravy, mean_charge, gravy_scale, charge_scale,
                        policy=resolved_policy,
                    )
                )
                for i1, item1 in enumerate(pools[p1]):
                    for i2, item2 in enumerate(pools[p2]):
                        cand1 = pools[p1][:i1] + [item2] + pools[p1][i1 + 1:]
                        cand2 = pools[p2][:i2] + [item1] + pools[p2][i2 + 1:]
                        c12_after = (
                            _pool_cost(
                                cand1, mean_gravy, mean_charge, gravy_scale, charge_scale,
                                policy=resolved_policy,
                            )
                            + _pool_cost(
                                cand2, mean_gravy, mean_charge, gravy_scale, charge_scale,
                                policy=resolved_policy,
                            )
                        )
                        gain = c12_before - c12_after
                        if gain > best_gain:
                            best_gain = gain
                            best_swap = (p1, i1, p2, i2)
        if best_swap is not None:
            p1, i1, p2, i2 = best_swap
            pools[p1][i1], pools[p2][i2] = pools[p2][i2], pools[p1][i1]
            improved = True
        if not improved:
            break

    # Canonical deterministic numbering: sort pools by lowest original construct
    # rank index so Pool 1 always contains the top-ranked emitted construct.
    for pool in pools:
        pool.sort(key=lambda p: p.index)
    pools.sort(key=lambda pool: pool[0].index)

    policy_dict = resolved_policy.to_dict()
    for pool_num, pool in enumerate(pools, start=1):
        for prof in pool:
            construct = constructs[prof.index]
            construct.components['pool'] = pool_num
            construct.components['pool_profile'] = {
                'strong_alleles': list(prof.strong_alleles),
                'medium_alleles': list(prof.medium_alleles),
                'net_charge': prof.net_charge,
                'net_charge_ph7': prof.net_charge if resolved_policy.charge_ph == 7.0 else net_charge_ph7(
                    construct.sequence,
                    n_terminal_acetylation=bool(construct.components.get('n_terminal_acetylation', False)),
                    c_terminal_amidation=bool(construct.components.get('c_terminal_amidation', False)),
                ),
                'charge_ph': resolved_policy.charge_ph,
            }
            construct.components['pool_policy'] = policy_dict

    return summarize_peptide_pools(
        constructs,
        num_pools_requested=num_pools,
        policy=resolved_policy,
    )


def summarize_peptide_pools(
    constructs,
    *,
    num_pools_requested: Optional[int] = None,
    policy=None,
) -> Optional[dict]:
    """Summarize pool assignments, MHC competition reduction, balance, and policy."""
    pooled = [c for c in constructs if c.components.get('pool') is not None]
    if not pooled:
        return None

    recorded_policy = None
    if policy is not None:
        recorded_policy = PeptidePoolPolicy.resolve(policy).to_dict()
    else:
        for c in pooled:
            if c.components.get('pool_policy'):
                recorded_policy = PeptidePoolPolicy.resolve(c.components['pool_policy']).to_dict()
                break
    if recorded_policy is None:
        recorded_policy = PeptidePoolPolicy().to_dict()
    charge_ph = float(recorded_policy.get('charge_ph', 7.0))

    by_pool: dict[int, list] = {}
    all_allele_counts = Counter()
    for c in pooled:
        p = int(c.components['pool'])
        by_pool.setdefault(p, []).append(c)
        profile = c.components.get('pool_profile') or {}
        alleles = set(profile.get('strong_alleles') or ()) | set(profile.get('medium_alleles') or ())
        for a in alleles:
            all_allele_counts[a] += 1

    single_pool_collisions = sum(cnt * (cnt - 1) // 2 for cnt in all_allele_counts.values())
    partitioned_collisions = 0
    pool_rows = []
    for pool_num in sorted(by_pool):
        members = by_pool[pool_num]
        allele_counts = Counter()
        strong_counts = Counter()
        cys_peptides = 0
        gravys = []
        charges = []
        for c in members:
            profile = c.components.get('pool_profile') or {}
            strong = set(profile.get('strong_alleles') or ())
            medium = set(profile.get('medium_alleles') or ())
            for a in strong | medium:
                allele_counts[a] += 1
            for a in strong:
                strong_counts[a] += 1
            mfg = c.manufacturability or {}
            if int(mfg.get('cysteine_count', 0) or 0) > 0:
                cys_peptides += 1
            if 'max_7mer_gravy_score' in mfg and mfg['max_7mer_gravy_score'] is not None:
                gravys.append(float(mfg['max_7mer_gravy_score']))
            if 'net_charge' in profile:
                charges.append(float(profile['net_charge']))
            elif 'net_charge_ph7' in profile and charge_ph == 7.0:
                charges.append(float(profile['net_charge_ph7']))
            else:
                charges.append(net_charge_at_ph(
                    c.sequence,
                    ph=charge_ph,
                    n_terminal_acetylation=bool(c.components.get('n_terminal_acetylation', False)),
                    c_terminal_amidation=bool(c.components.get('c_terminal_amidation', False)),
                ))
        collisions = sum(cnt * (cnt - 1) // 2 for cnt in allele_counts.values())
        partitioned_collisions += collisions
        pool_rows.append({
            'pool': pool_num,
            'count': len(members),
            'construct_names': [c.name for c in members],
            'allele_counts': dict(sorted(allele_counts.items())),
            'strong_allele_counts': dict(sorted(strong_counts.items())),
            'same_allele_pairs': collisions,
            'cys_peptides': cys_peptides,
            'mean_max_gravy': round(sum(gravys) / len(gravys), 2) if gravys else None,
            'mean_net_charge': round(sum(charges) / len(charges), 2) if charges else None,
        })

    return {
        'num_pools_requested': int(num_pools_requested or len(pool_rows)),
        'num_pools_assigned': len(pool_rows),
        'single_pool_same_allele_pairs': single_pool_collisions,
        'partitioned_same_allele_pairs': partitioned_collisions,
        'charge_ph': charge_ph,
        'policy': recorded_policy,
        'pools': pool_rows,
    }
