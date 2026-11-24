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
3. **Hydrophobicity (max 7-mer GRAVY)** and **net charge at pH 7.0** are
   balanced across pools to reduce co-precipitation and solubility failures.

The optimizer is deterministic, species-agnostic (uses ``mhcgnomes``-backed
canonical allele names), enforces balanced pool cardinalities
(``floor(N/K)..ceil(N/K)``), and leaves single-pool (``num_pools == 1``)
designs untouched.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from typing import Iterable, Optional

from .coverage import compute_coverage
from .identifiers import normalize_mhc_allele


# Tier weights for same-allele competition within a pool.
_STRONG_WEIGHT = 3.0
_MEDIUM_WEIGHT = 1.0

# Relative weights between the primary (MHC competition), secondary (cysteine
# cross-linking), and tertiary (scale-normalized GRAVY & net-charge balance)
# objectives.
_W_ALLELE = 100.0
_W_CYSTEINE = 10.0
_W_GRAVY = 2.0
_W_CHARGE = 1.0


def net_charge_ph7(
    sequence: str,
    *,
    n_terminal_acetylation: bool = False,
    c_terminal_amidation: bool = False,
) -> float:
    """Estimate peptide net charge at physiological pH (7.0).

    Uses standard side-chain contributions (K/R +1.0, H +0.1, D/E -1.0) plus
    free-terminus contributions (+1.0 for unacetylated N-terminus, -1.0 for
    unamidated C-terminus).
    """
    if not sequence:
        return 0.0
    seq = sequence.upper()
    charge = 0.0
    for aa in seq:
        if aa in ('K', 'R'):
            charge += 1.0
        elif aa == 'H':
            charge += 0.1
        elif aa in ('D', 'E'):
            charge -= 1.0
    if not n_terminal_acetylation:
        charge += 1.0
    if not c_terminal_amidation:
        charge -= 1.0
    return round(charge, 2)


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

    epitopes = []
    for vp in matched:
        for ep in getattr(vp, 'epitopes', ()) or ():
            if getattr(ep, 'is_targetable', True):
                epitopes.append(ep)
    return epitopes


def _allele_weights_for_epitopes(
    epitopes,
    target_alleles: Optional[Iterable[str]] = None,
) -> tuple[dict[str, float], tuple[str, ...], tuple[str, ...]]:
    """Map each bound allele to a competition weight (strong=3.0, medium=1.0).

    Prefers percentile-rank tiers via :func:`vaxrank.coverage.compute_coverage`
    (presentation or affinity <=0.5% strong, <=2.0% medium/low). Falls back
    to ``per_allele_scores`` (>=0.5 strong, >=0.1 medium) when percentile ranks
    are absent (e.g. IC50-only or synthetic runs).
    """
    alleles = []
    seen = set()
    for a in target_alleles or ():
        canon = _canonical_allele(a)
        if canon and canon not in seen:
            seen.add(canon)
            alleles.append(canon)
    for ep in epitopes:
        for a in getattr(ep, 'per_allele_scores', {}) or {}:
            canon = _canonical_allele(a)
            if canon and canon not in seen:
                seen.add(canon)
                alleles.append(canon)
        for pred in ep.predictions_flat() if hasattr(ep, 'predictions_flat') else ():
            if pred.allele:
                canon = _canonical_allele(pred.allele)
                if canon and canon not in seen:
                    seen.add(canon)
                    alleles.append(canon)

    if not alleles or not epitopes:
        return {}, (), ()

    coverage = compute_coverage(epitopes, alleles)
    strong = set()
    medium = set()
    any_pct_tier = any(cov.tier is not None for cov in coverage.values())
    if any_pct_tier:
        for allele, cov in coverage.items():
            if cov.tier == 'strong':
                strong.add(allele)
            elif cov.tier in ('medium', 'low'):
                medium.add(allele)
    else:
        best_score = {}
        for ep in epitopes:
            for allele, score in (getattr(ep, 'per_allele_scores', {}) or {}).items():
                if score is None:
                    continue
                canon = _canonical_allele(allele)
                prev = best_score.get(canon, 0.0)
                if score > prev:
                    best_score[canon] = float(score)
        for allele, score in best_score.items():
            if score >= 0.5:
                strong.add(allele)
            elif score >= 0.1:
                medium.add(allele)

    medium -= strong
    weights = {a: _STRONG_WEIGHT for a in strong}
    weights.update({a: _MEDIUM_WEIGHT for a in medium})
    return weights, tuple(sorted(strong)), tuple(sorted(medium))


def _build_profiles(
    constructs,
    ranked_vaccine_peptides=(),
    target_alleles=None,
) -> list[PeptidePoolProfile]:
    vp_by_sequence: dict[str, list] = {}
    for _, vps in ranked_vaccine_peptides or ():
        for vp in vps:
            vp_by_sequence.setdefault(vp.amino_acids, []).append(vp)

    profiles = []
    for idx, construct in enumerate(constructs):
        epitopes = _extract_construct_epitopes(construct, vp_by_sequence)
        weights, strong, medium = _allele_weights_for_epitopes(epitopes, target_alleles)
        mfg = construct.manufacturability or {}
        max_gravy = float(mfg.get('max_7mer_gravy_score', 0.0) or 0.0)
        cterm_gravy = float(mfg.get('cterm_7mer_gravy_score', 0.0) or 0.0)
        cys_count = int(mfg.get('cysteine_count', 0) or 0)
        charge = net_charge_ph7(
            construct.sequence,
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
) -> float:
    """Compute objective penalty for one pool."""
    if not pool_profiles:
        return 0.0

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
        _W_ALLELE * allele_pen
        + _W_CYSTEINE * cys_pen
        + _W_GRAVY * gravy_dev
        + _W_CHARGE * charge_dev
    )


def _total_cost(
    pools: list[list[PeptidePoolProfile]],
    mean_gravy: float,
    mean_charge: float,
    gravy_scale: float = 1.0,
    charge_scale: float = 1.0,
) -> float:
    return sum(
        _pool_cost(p, mean_gravy, mean_charge, gravy_scale, charge_scale)
        for p in pools
    )


def _same_allele_pairs(pool_profiles: list[PeptidePoolProfile]) -> int:
    counts = Counter()
    for prof in pool_profiles:
        for allele in prof.all_alleles:
            counts[allele] += 1
    return sum(c * (c - 1) // 2 for c in counts.values())


def partition_peptide_constructs(
    constructs,
    num_pools: int,
    *,
    ranked_vaccine_peptides=(),
    target_alleles=None,
) -> Optional[dict]:
    """Assign emitted peptide constructs across ``1..num_pools`` formulation pools.

    Mutates each construct's ``components`` dict in place to add ``'pool'``
    (1-based integer) and ``'pool_profile'`` (bound alleles + net charge) when
    ``num_pools > 1``. Returns a JSON-serializable pool summary dict, or
    ``None`` when ``num_pools <= 1`` or ``constructs`` is empty.
    """
    if isinstance(num_pools, bool) or not isinstance(num_pools, int) or num_pools < 1:
        raise ValueError("num_pools must be a positive integer; got %r" % (num_pools,))
    if num_pools == 1 or not constructs:
        return None

    k = min(num_pools, len(constructs))
    profiles = _build_profiles(
        constructs,
        ranked_vaccine_peptides=ranked_vaccine_peptides,
        target_alleles=target_alleles,
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
            before = _pool_cost(pools[pi], mean_gravy, mean_charge, gravy_scale, charge_scale)
            after = _pool_cost(pools[pi] + [prof], mean_gravy, mean_charge, gravy_scale, charge_scale)
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
    max_passes = 50
    for _ in range(max_passes):
        improved = False
        best_swap = None
        best_gain = 1e-12
        for p1 in range(k):
            for p2 in range(p1 + 1, k):
                c12_before = (
                    _pool_cost(pools[p1], mean_gravy, mean_charge, gravy_scale, charge_scale)
                    + _pool_cost(pools[p2], mean_gravy, mean_charge, gravy_scale, charge_scale)
                )
                for i1, item1 in enumerate(pools[p1]):
                    for i2, item2 in enumerate(pools[p2]):
                        cand1 = pools[p1][:i1] + [item2] + pools[p1][i1 + 1:]
                        cand2 = pools[p2][:i2] + [item1] + pools[p2][i2 + 1:]
                        c12_after = (
                            _pool_cost(cand1, mean_gravy, mean_charge, gravy_scale, charge_scale)
                            + _pool_cost(cand2, mean_gravy, mean_charge, gravy_scale, charge_scale)
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

    for pool_num, pool in enumerate(pools, start=1):
        for prof in pool:
            construct = constructs[prof.index]
            construct.components['pool'] = pool_num
            construct.components['pool_profile'] = {
                'strong_alleles': list(prof.strong_alleles),
                'medium_alleles': list(prof.medium_alleles),
                'net_charge_ph7': prof.net_charge,
            }

    return summarize_peptide_pools(
        constructs,
        num_pools_requested=num_pools,
    )


def summarize_peptide_pools(
    constructs,
    *,
    num_pools_requested: Optional[int] = None,
) -> Optional[dict]:
    """Summarize pool assignments, MHC competition reduction, and balance."""
    pooled = [c for c in constructs if c.components.get('pool') is not None]
    if not pooled:
        return None

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
            if 'net_charge_ph7' in profile:
                charges.append(float(profile['net_charge_ph7']))
            else:
                charges.append(net_charge_ph7(
                    c.sequence,
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
        'pools': pool_rows,
    }
