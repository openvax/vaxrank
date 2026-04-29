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
These "chimeric" peptides dilute the immune response.

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
from dataclasses import dataclass, field

from .vaccine_library import JUNCTION_SWAP_CANDIDATES, get_linker

logger = logging.getLogger(__name__)


# Default presentation-rank thresholds. mhcflurry-presentation reports
# percentile rank where lower = stronger predicted presentation. The
# convention in immunoinformatics is "rank ≤ 2.0" for "presented
# binder" and "rank ≤ 0.5" for "strong binder". The optimizer
# minimizes count of rank ≤ STRONG first, then rank ≤ MILD as
# tie-breaker, then maximizes the worst-case rank.
RANK_STRONG = 0.5   # high-confidence presentation
RANK_MILD = 2.0     # any presentation


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
    """Outcome of the per-junction linker swap optimization."""
    chosen_linker_per_junction: list  # list[Linker]
    burden: int  # count of rank ≤ RANK_MILD across all junctions, post-swap
    strong_burden: int  # count of rank ≤ RANK_STRONG, post-swap
    chimeric_predictions: list = field(default_factory=list)  # all post-swap

    def linker_names(self):
        return [link.name for link in self.chosen_linker_per_junction]


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


def _score_kmers(kmers, alleles, predictor):
    """Run the predictor and return rows of (kmer, allele, rank)."""
    if not kmers:
        return []
    # Predictor API: mhctools.BasePredictor returns BindingPrediction
    # objects with .peptide / .allele / .percentile_rank attributes.
    predictions = predictor.predict_peptides(kmers)
    rows = []
    allele_set = set(alleles) if alleles else None
    for p in predictions:
        if allele_set and p.allele not in allele_set:
            continue
        rank = getattr(p, 'percentile_rank', None)
        if rank is None:
            continue
        rows.append((p.peptide, p.allele, float(rank)))
    return rows


def _burden_key(rows, rank_strong=RANK_STRONG, rank_mild=RANK_MILD):
    """Tuple of (strong_count, mild_count, worst_rank_negated) used to
    pick the linker with lowest junction burden. Lower = better.

    Tie-breaking order:
      1. Fewer rank ≤ rank_strong hits (high-confidence presentation)
      2. Fewer rank ≤ rank_mild hits
      3. Higher worst-case rank (push the worst hit's rank as high as
         possible — i.e. minimize negated worst-rank)
    """
    strong = sum(1 for r in rows if r[2] <= rank_strong)
    mild = sum(1 for r in rows if r[2] <= rank_mild)
    worst = min((r[2] for r in rows), default=100.0)
    return (strong, mild, -worst)


def optimize_junction_linkers(
        antigen_aas, alleles, predictor,
        candidate_names=JUNCTION_SWAP_CANDIDATES,
        k_lengths=(8, 9, 10, 11),
        reference_proteome=None,
        rank_strong=RANK_STRONG, rank_mild=RANK_MILD):
    """Pick the candidate linker at each junction that minimizes
    predicted MHC presentation of chimeric k-mers.

    Parameters
    ----------
    antigen_aas : list[str]
        Antigen amino-acid sequences in the order they'll be
        concatenated. There are ``len(antigen_aas) - 1`` junctions.
    alleles : list[str]
        Patient HLA alleles, e.g. ``['HLA-A*02:01', 'HLA-B*07:02']``.
    predictor : mhctools.BasePredictor
        An MHC predictor with a ``predict_peptides(peptides)`` method
        (e.g. ``mhcflurry.MHCflurryPredictor``). The predictor's
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

    Returns
    -------
    JunctionSwapResult
    """
    if len(antigen_aas) <= 1:
        return JunctionSwapResult(
            chosen_linker_per_junction=[],
            burden=0, strong_burden=0, chimeric_predictions=[])

    candidates = [get_linker(n) for n in candidate_names]

    chosen = []
    all_predictions = []
    total_burden = 0
    total_strong = 0

    for j in range(len(antigen_aas) - 1):
        left_aa = antigen_aas[j]
        right_aa = antigen_aas[j + 1]
        best = None
        for cand in candidates:
            kmers = junction_kmers(
                left_aa, cand.amino_acids, right_aa, k_lengths,
                reference_proteome=reference_proteome)
            rows = _score_kmers(kmers, alleles, predictor)
            key = _burden_key(rows, rank_strong, rank_mild)
            if best is None or key < best[0]:
                best = (key, cand, rows)
        key, cand, rows = best
        strong, mild, _ = key
        chosen.append(cand)
        total_burden += mild
        total_strong += strong
        for kmer, allele, rank in rows:
            all_predictions.append(JunctionPrediction(
                junction_index=j, linker_name=cand.name,
                kmer=kmer, allele=allele, rank=rank))
        logger.info(
            "Junction %d: chose linker %s (strong=%d, mild=%d).",
            j, cand.name, strong, mild)

    return JunctionSwapResult(
        chosen_linker_per_junction=chosen,
        burden=total_burden,
        strong_burden=total_strong,
        chimeric_predictions=all_predictions,
    )
