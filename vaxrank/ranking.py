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
Registry of vaccine-peptide ranking rules used to build the lexicographic
sort key on ``VaccinePeptide``. Each rule takes a ``VaccinePeptide`` and
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


def _mutant_epitope_score_rule(peptide):
    # Sum of normalized MHC binding scores across kept mutant epitopes.
    # Rounded to 6 digits so floating-point noise can't act as a tiebreaker.
    return -round(peptide.mutant_epitope_score, 6)


def _n_alt_reads_rule(peptide):
    # Number of reads supporting the variant (RNA evidence).
    return -peptide.mutant_protein_fragment.n_alt_reads


def _n_alt_reads_supporting_rule(peptide):
    # Reads spanning the specific protein-coding sequence selected for
    # this vaccine peptide (a subset of n_alt_reads).
    return -peptide.mutant_protein_fragment.n_alt_reads_supporting_protein_sequence


def _wildtype_epitope_score_rule(peptide):
    # Non-mutant MHC binding score — we want this SMALL (already positive-signed).
    return round(peptide.wildtype_epitope_score, 6)


def _n_mutant_amino_acids_rule(peptide):
    # Prefer peptides containing more mutant residues.
    return -peptide.mutant_protein_fragment.n_mutant_amino_acids


def _mutation_distance_from_edge_rule(peptide):
    # All else equal, center the mutation inside the vaccine peptide.
    return -peptide.mutant_protein_fragment.mutation_distance_from_edge


RANKING_RULE_REGISTRY = {
    "mutant_epitope_score": _mutant_epitope_score_rule,
    "n_alt_reads": _n_alt_reads_rule,
    "n_alt_reads_supporting": _n_alt_reads_supporting_rule,
    "wildtype_epitope_score": _wildtype_epitope_score_rule,
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
    "mutant_epitope_score",
    "n_alt_reads",
    MANUFACTURABILITY_SENTINEL,
    "n_alt_reads_supporting",
    "wildtype_epitope_score",
    "n_mutant_amino_acids",
    "mutation_distance_from_edge",
)


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
