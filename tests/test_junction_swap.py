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

"""Tests for junction-aware linker swap (vaxrank issue #247)."""

from types import SimpleNamespace

from vaxrank.junction_swap import (
    JunctionSwapResult,
    junction_kmers,
    optimize_junction_linkers,
)


class StubPredictor:
    """Fake predictor that scores k-mers against a fixed dict.

    rank_table maps (kmer, allele) -> percentile_rank. Anything not in
    the table gets a high (= weak presentation) rank.
    """
    def __init__(self, rank_table=None, default_rank=99.0):
        self.rank_table = rank_table or {}
        self.default_rank = default_rank

    def predict_peptides(self, peptides):
        out = []
        seen_alleles = set()
        for kmer, allele in self.rank_table.keys():
            seen_alleles.add(allele)
        if not seen_alleles:
            seen_alleles = {"HLA-A*02:01"}
        for kmer in peptides:
            for allele in seen_alleles:
                rank = self.rank_table.get(
                    (kmer, allele), self.default_rank)
                out.append(SimpleNamespace(
                    peptide=kmer, allele=allele, percentile_rank=rank))
        return out


def test_junction_kmers_only_spans_junction():
    out = junction_kmers("AAAAAAAAAA", "GGGGS", "BBBBBBBBBB", k_lengths=(9,))
    # Each k-mer must touch the linker. With 9-mers and 8-aa context on
    # each side + 5-aa linker, there should be 13 spanning k-mers.
    assert all("G" in k or "S" in k for k in out)
    # No k-mer entirely inside the antigens
    assert all(not k.startswith("AAAAAAAAA") for k in out)
    assert all(not k.startswith("BBBBBBBBB") for k in out)


def test_junction_kmers_filters_reference_proteome():
    out_no_filter = junction_kmers(
        "AAAAAAAAAA", "GGGGS", "BBBBBBBBBB", k_lengths=(9,))
    # Drop one of the spanning k-mers via reference proteome
    drop = out_no_filter[0]
    out_with_filter = junction_kmers(
        "AAAAAAAAAA", "GGGGS", "BBBBBBBBBB", k_lengths=(9,),
        reference_proteome={drop})
    assert drop in out_no_filter
    assert drop not in out_with_filter
    assert len(out_with_filter) == len(out_no_filter) - 1


def test_optimize_picks_zero_burden_linker_when_one_exists():
    # Construct: two antigens. The (G4S)2 junction produces a strong
    # binder; AAA does not. Optimizer should pick AAA.
    antigens = ["KLQGHSAPVL", "DVIVNCDESLLAS"]
    # Compute the (G4S)2 junction k-mers and mark one as a strong hit.
    g4s2_kmers = junction_kmers(
        antigens[0], "GGGGSGGGGS", antigens[1], k_lengths=(9,))
    strong_kmer = g4s2_kmers[0]
    rank_table = {(strong_kmer, "HLA-A*02:01"): 0.1}
    predictor = StubPredictor(rank_table=rank_table)

    result = optimize_junction_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("(G4S)2", "AAA"),
        k_lengths=(9,))
    assert isinstance(result, JunctionSwapResult)
    assert len(result.chosen_linker_per_junction) == 1
    chosen = result.chosen_linker_per_junction[0]
    assert chosen.name == "AAA", \
        "Optimizer should swap to AAA when (G4S)2 has a strong hit"
    assert result.strong_burden == 0


def test_optimize_keeps_default_when_all_candidates_equal():
    # Predictor returns nothing strong for any k-mer. All candidates
    # tie at zero burden; the loop returns the first candidate scanned.
    antigens = ["KLQGH", "MNNVD"]
    predictor = StubPredictor(rank_table={})  # all default rank 99
    result = optimize_junction_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("(G4S)2", "AAA"),
        k_lengths=(9,))
    assert result.strong_burden == 0
    assert result.burden == 0


def test_optimize_handles_single_antigen():
    # No junctions → empty result.
    result = optimize_junction_linkers(
        ["KLQGH"], alleles=["HLA-A*02:01"], predictor=StubPredictor(),
        candidate_names=("(G4S)2",))
    assert result.chosen_linker_per_junction == []
    assert result.burden == 0


def test_optimize_records_per_junction_choices():
    # Three antigens → two junctions. Set up so the first junction
    # prefers AAA (strong (G4S)2 hit) and the second prefers (G4S)2
    # (strong AAA hit). Verify each junction gets the right pick.
    antigens = ["KLQGHSAP", "VLDVIVNC", "DESLLASD"]
    g4s2_j1 = junction_kmers(antigens[0], "GGGGSGGGGS", antigens[1], (9,))
    aaa_j2 = junction_kmers(antigens[1], "AAA", antigens[2], (9,))
    rank_table = {
        (g4s2_j1[0], "HLA-A*02:01"): 0.05,  # bad for (G4S)2 at j=0
        (aaa_j2[0], "HLA-A*02:01"): 0.05,   # bad for AAA at j=1
    }
    predictor = StubPredictor(rank_table=rank_table)
    result = optimize_junction_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("(G4S)2", "AAA"),
        k_lengths=(9,))
    names = result.linker_names()
    assert len(names) == 2
    assert names[0] == "AAA"
    assert names[1] == "(G4S)2"


def test_optimize_minimizes_strong_first_then_mild():
    # Linker A: 0 strong, 5 mild
    # Linker B: 1 strong, 0 mild
    # Strong-first ranking should pick A (fewer strong hits) even
    # though it has more mild hits.
    antigens = ["KKKKKK", "RRRRRR"]
    # Pick a candidate set where (G4S) and AAA give different mild/strong
    # patterns — we'll fake this by manipulating the rank table.
    g4s_kmers = junction_kmers(antigens[0], "GGGGS", antigens[1], (9,))
    aaa_kmers = junction_kmers(antigens[0], "AAA", antigens[1], (9,))
    rank_table = {}
    # G4S: 0 strong, 5 mild (rank 1.5 each)
    for k in g4s_kmers[:5]:
        rank_table[(k, "HLA-A*02:01")] = 1.5
    # AAA: 1 strong (rank 0.3), no others mild
    rank_table[(aaa_kmers[0], "HLA-A*02:01")] = 0.3
    predictor = StubPredictor(rank_table=rank_table)
    result = optimize_junction_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("G4S", "AAA"),
        k_lengths=(9,))
    chosen = result.chosen_linker_per_junction[0]
    assert chosen.name == "G4S", \
        "Strong-first ranking should beat mild-count ranking"


def test_optimize_filters_reference_proteome_kmers():
    # If a chimeric k-mer matches the reference proteome, it must be
    # excluded from the burden count.
    antigens = ["KLQGH", "MNNVD"]
    g4s_kmers = junction_kmers(antigens[0], "GGGGS", antigens[1], (9,))
    rank_table = {(g4s_kmers[0], "HLA-A*02:01"): 0.05}  # would be strong
    predictor = StubPredictor(rank_table=rank_table)
    result = optimize_junction_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("G4S",),
        k_lengths=(9,),
        reference_proteome={g4s_kmers[0]})  # mark as already-tolerated
    # Strong burden should be 0 because the only strong hit was filtered.
    assert result.strong_burden == 0
