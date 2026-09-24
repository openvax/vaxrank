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

"""
Final construct ranking and the rules for selecting vaccine-peptide windows.

``rank_constructs`` orders selected constructs across all input sources.
The rule registry builds the lexicographic window-selection sort key on
``VaccinePeptide``. Each rule takes a ``VaccinePeptide`` and
returns a single numeric value to minimize — Python's ascending sort on
tuples then gives us the ordering we want.

``DEFAULT_RANKING_RULES`` is pinned byte-for-byte to the legacy
lexicographic sort produced by ``VaccinePeptide.lexicographic_sort_key``
before rules became data-driven. Do not reorder without updating the
parametrized parity test in ``tests/test_ranking.py``.

The special sentinel ``"manufacturability"`` expands in place to the
peptide's manufacturability sort tuple (respecting whatever
``manufacturability_rules`` / ``manufacturability_thresholds`` the peptide
was built with). That lets users compose peptide-level ranking criteria
with manufacturability criteria in a single flat list, e.g. to
demote manufacturability below tie-breakers, or drop it entirely.
"""

from itertools import groupby
from math import isfinite


def _rna_tiebreak_evidence(peptide):
    """Return a count and its stated semantics, or None when unavailable."""
    fragment = peptide.mutant_protein_fragment
    if fragment is None:
        return None
    method = getattr(fragment, "rna_evidence_method", "")
    subject = getattr(fragment, "rna_evidence_subject", "")
    count = getattr(fragment, "n_rna_alt", None)
    if not method or not subject or count is None or not isfinite(count):
        return None
    return (method, subject), count


def rank_constructs(source_peptides):
    """Rank already-selected constructs identically for every input source.

    Each entry is ``(source, [VaccinePeptide, ...])``. The first peptide is
    the representative chosen by upstream occurrence/window selection; this
    function does not reorder its alternative windows. Its configured
    ``combined_score`` determines final descending rank.

    Within an exact score tie, use descending RNA support only if every
    representative has a count with the same stated method and unit. Missing
    evidence, different units/derivations, and source-agnostic antigens skip
    that tier for the whole tie, rather than comparing incomparable counts
    or replacing unknown support with zero. Descending target-epitope score
    breaks remaining ties; complete ties preserve input order. Empty peptide
    lists sort last, including when valid constructs have negative scores.
    """
    scored, empty = [], []
    for entry in source_peptides:
        if entry[1]:
            scored.append((entry[1][0].combined_score, entry))
        else:
            empty.append(entry)
    scored.sort(key=lambda item: item[0], reverse=True)
    ranked = []
    for _, group in groupby(scored, key=lambda item: item[0]):
        entries = [entry for _, entry in group]
        evidence = [_rna_tiebreak_evidence(entry[1][0]) for entry in entries]
        comparable = all(item is not None for item in evidence) and len({
            item[0] for item in evidence if item is not None
        }) == 1
        ordered = sorted(
            zip(entries, evidence),
            key=lambda item: (
                (item[1][1], item[0][1][0].target_epitope_score)
                if comparable else (item[0][1][0].target_epitope_score,)
            ),
            reverse=True,
        )
        ranked.extend(entry for entry, _ in ordered)
    return ranked + empty


def _target_epitope_score_rule(peptide):
    # Sum of normalized MHC binding scores across kept mutant epitopes.
    # Rounded to 6 digits so floating-point noise can't act as a tiebreaker.
    return -round(peptide.target_epitope_score, 6)


def _n_alt_reads_rule(peptide):
    # Number of reads supporting the variant (RNA evidence).
    return -peptide.mutant_protein_fragment.n_alt_reads


def _n_alt_reads_supporting_rule(peptide):
    # Reads spanning the specific protein-coding sequence selected for
    # this vaccine peptide (a subset of n_alt_reads).
    return -peptide.mutant_protein_fragment.n_alt_reads_supporting_protein_sequence


def _self_epitope_score_rule(peptide):
    # Non-mutant MHC binding score — we want this SMALL (already positive-signed).
    return round(peptide.self_epitope_score, 6)


def _n_mutant_amino_acids_rule(peptide):
    # Prefer peptides containing more mutant residues.
    return -peptide.mutant_protein_fragment.n_mutant_amino_acids


def _mutation_distance_from_edge_rule(peptide):
    # All else equal, center the mutation inside the vaccine peptide.
    return -peptide.mutant_protein_fragment.mutation_distance_from_edge


RANKING_RULE_REGISTRY = {
    "target_epitope_score": _target_epitope_score_rule,
    "n_alt_reads": _n_alt_reads_rule,
    "n_alt_reads_supporting": _n_alt_reads_supporting_rule,
    "self_epitope_score": _self_epitope_score_rule,
    "n_mutant_amino_acids": _n_mutant_amino_acids_rule,
    "mutation_distance_from_edge": _mutation_distance_from_edge_rule,
    # "manufacturability" is a sentinel handled by compute_ranking_tuple
    # — it expands inline using the peptide's own manufacturability
    # rules + thresholds. Listed here so registry introspection and
    # validation can still confirm it's a recognized name.
    "manufacturability": None,
}

MANUFACTURABILITY_SENTINEL = "manufacturability"

# Legacy ordering preserved byte-for-byte: this is the exact concatenation
# `lexicographic_sort_key` produced before rules were made configurable —
# essential (2 rules) + manufacturability (expanded inline) + extra (4 rules).
# Do not reorder without updating the parametrized parity test.
DEFAULT_RANKING_RULES = (
    "target_epitope_score",
    "n_alt_reads",
    MANUFACTURABILITY_SENTINEL,
    "n_alt_reads_supporting",
    "self_epitope_score",
    "n_mutant_amino_acids",
    "mutation_distance_from_edge",
)

MUTATION_SPECIFIC_RANKING_RULES = frozenset({
    "n_alt_reads",
    "n_alt_reads_supporting",
    "n_mutant_amino_acids",
    "mutation_distance_from_edge",
})


def compute_ranking_tuple(peptide, rules=None):
    """Apply the ordered rule list against a ``VaccinePeptide`` and return
    the tuple used as a stable lexicographic sort key.

    ``rules`` defaults to ``DEFAULT_RANKING_RULES`` (the legacy order).
    Each entry must be a key in ``RANKING_RULE_REGISTRY``. The special
    sentinel ``"manufacturability"`` expands inline to the peptide's
    manufacturability sort tuple.
    """
    if rules is None:
        rules = DEFAULT_RANKING_RULES
    values = []
    for rule_name in rules:
        if rule_name == MANUFACTURABILITY_SENTINEL:
            # Let the peptide build its manufacturability sort tuple with
            # whatever rules / thresholds it was configured with — keep
            # the two registries decoupled.
            manufacturability_tuple = peptide.peptide_synthesis_difficulty_score_tuple(
                rules=peptide.manufacturability_rules,
                **peptide.manufacturability_thresholds,
            )
            values.extend(manufacturability_tuple)
            continue
        try:
            fn = RANKING_RULE_REGISTRY[rule_name]
        except KeyError:
            raise ValueError(
                f"Unknown ranking rule '{rule_name}'. "
                f"Available: {sorted(RANKING_RULE_REGISTRY)}"
            ) from None
        if fn is None:
            # Defensive: only "manufacturability" has a None mapping and it's
            # handled above. Any other None is a bug in the registry.
            raise ValueError(
                f"Ranking rule '{rule_name}' has no scoring function — "
                f"likely a registry setup bug."
            )
        values.append(fn(peptide))
    return tuple(values)
