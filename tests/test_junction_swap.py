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

"""Tests for per-junction linker optimization (vaxrank issue #247)."""

from types import SimpleNamespace

from vaxrank.junction_swap import (
    JunctionSwapResult,
    junction_kmers,
    optimize_linkers,
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

    result = optimize_linkers(
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
    result = optimize_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("(G4S)2", "AAA"),
        k_lengths=(9,))
    assert result.strong_burden == 0
    assert result.burden == 0


def test_optimize_handles_single_antigen():
    # No junctions → empty result.
    result = optimize_linkers(
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
    result = optimize_linkers(
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
    result = optimize_linkers(
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
    result = optimize_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("G4S",),
        k_lengths=(9,),
        reference_proteome={g4s_kmers[0]})  # mark as already-tolerated
    # Strong burden should be 0 because the only strong hit was filtered.
    assert result.strong_burden == 0


# ---- end-to-end integration with assemble_mrna_constructs ----------------

def test_assemble_with_optimizer_uses_chosen_linkers():
    """End-to-end: assemble_mrna_constructs(optimize_linkers=True,
    ...) must actually emit a construct whose protein has the
    per-junction chosen linkers between the right antigens, AND the
    manifest must record the swap.
    """
    from varcode import Variant
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    # Two antigens; default (G4S)2 produces a strong hit, AAA does not.
    a1 = "KLQGHSAPVL"
    a2 = "DVIVNCDESLLAS"
    g4s2_kmers = junction_kmers(a1, "GGGGSGGGGS", a2, k_lengths=(9,))
    rank_table = {(g4s2_kmers[0], "HLA-A*02:01"): 0.05}
    predictor = StubPredictor(rank_table=rank_table)

    fragment_a = SimpleNamespace(
        amino_acids=a1, gene_name='GENEA',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=len(a1))
    fragment_b = SimpleNamespace(
        amino_acids=a2, gene_name='GENEB',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=len(a2))
    pep_a = SimpleNamespace(
        mutant_protein_fragment=fragment_a, mutant_epitopes=[])
    pep_b = SimpleNamespace(
        mutant_protein_fragment=fragment_b, mutant_epitopes=[])
    pairs = [
        (Variant('1', 100, 'A', 'T'), [pep_a]),
        (Variant('2', 200, 'A', 'T'), [pep_b]),
    ]

    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        optimize_linkers=True, junction_swap_candidates=("(G4S)2", "AAA"),
        junction_kmer_lengths=(9,),
        antigens_per_construct=2, max_constructs=1,
        max_antigen_length_aa=20,
        utr_3p='HBB',
    )
    [c] = assemble_mrna_constructs(
        pairs, options=options,
        mhc_predictor=predictor, mhc_alleles=["HLA-A*02:01"])

    # The protein must contain AAA (chosen) between the two antigens.
    protein = c.components['protein']
    # Strip the prepended start M.
    body = protein[1:] if protein.startswith("M") else protein
    assert body == a1 + "AAA" + a2, \
        "Protein should join %s + AAA + %s, got %r" % (a1, a2, body)
    # Manifest annotation must report the swap.
    assert 'junction_swap' in c.components
    swap_meta = c.components['junction_swap']
    assert swap_meta['enabled'] is True
    assert swap_meta['chosen'] == ['AAA']
    assert swap_meta['burden_strong'] == 0
    assert swap_meta['default_burden_strong'] >= 1


def test_assemble_with_optimizer_falls_back_when_predictor_missing(caplog):
    """optimize_linkers=True without a predictor warns and uses
    the shared linker — no exception."""
    import logging
    from varcode import Variant
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    fragment_a = SimpleNamespace(
        amino_acids="KLQGH", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=5)
    fragment_b = SimpleNamespace(
        amino_acids="MNNVD", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=5)
    pep_a = SimpleNamespace(
        mutant_protein_fragment=fragment_a, mutant_epitopes=[])
    pep_b = SimpleNamespace(
        mutant_protein_fragment=fragment_b, mutant_epitopes=[])
    pairs = [
        (Variant('1', 100, 'A', 'T'), [pep_a]),
        (Variant('2', 200, 'A', 'T'), [pep_b]),
    ]

    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        optimize_linkers=True,
        antigens_per_construct=2, max_constructs=1,
        max_antigen_length_aa=10,
        utr_3p='HBB',
    )
    with caplog.at_level(logging.WARNING):
        [c] = assemble_mrna_constructs(pairs, options=options)
    # Should have warned and fallen back.
    assert any("falling back" in r.message.lower()
               for r in caplog.records), \
        "Expected a fallback warning when predictor is missing"
    assert c.components['junction_swap']['enabled'] is False


def test_assemble_default_kept_when_candidates_dont_help():
    """When no candidate beats the default, the construct keeps the
    default linker and the manifest annotates this."""
    from varcode import Variant
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    # Predictor returns no hits for any kmer → all candidates tie at
    # zero burden → default wins (no strict improvement).
    predictor = StubPredictor(rank_table={})

    fragment_a = SimpleNamespace(
        amino_acids="KLQGH", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=5)
    fragment_b = SimpleNamespace(
        amino_acids="MNNVD", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=5)
    pep_a = SimpleNamespace(
        mutant_protein_fragment=fragment_a, mutant_epitopes=[])
    pep_b = SimpleNamespace(
        mutant_protein_fragment=fragment_b, mutant_epitopes=[])
    pairs = [
        (Variant('1', 100, 'A', 'T'), [pep_a]),
        (Variant('2', 200, 'A', 'T'), [pep_b]),
    ]

    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        optimize_linkers=True, junction_swap_candidates=("(G4S)2", "AAA"),
        junction_kmer_lengths=(9,),
        antigens_per_construct=2, max_constructs=1,
        max_antigen_length_aa=10,
        utr_3p='HBB',
        linker="(G4S)2",
    )
    [c] = assemble_mrna_constructs(
        pairs, options=options,
        mhc_predictor=predictor, mhc_alleles=["HLA-A*02:01"])
    swap_meta = c.components['junction_swap']
    assert swap_meta['enabled'] is True
    assert "default" in swap_meta.get('note', '').lower()


def test_optimize_includes_mitd_junction_when_provided():
    # With mitd_aa supplied, the optimizer returns one extra entry
    # for the antigen-N → MITD junction.
    antigens = ["KLQGH", "MNNVD"]
    predictor = StubPredictor()
    result = optimize_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("(G4S)2",),
        k_lengths=(9,),
        mitd_aa="IVGIIAGLVLLG")
    # Two antigens + MITD = 2 inter-antigen junctions worth + 1 MITD = total
    # = (len(antigens) - 1) + 1 = 2
    assert len(result.chosen_linker_per_junction) == 2


def test_optimize_default_linker_burden_tracked_in_sweep():
    # When default_linker_name is supplied, the result must expose
    # default_strong_burden and default_burden so the caller can
    # decide whether to swap without a second sweep.
    antigens = ["KLQGH", "MNNVD"]
    predictor = StubPredictor()
    result = optimize_linkers(
        antigens, alleles=["HLA-A*02:01"], predictor=predictor,
        candidate_names=("(G4S)2", "AAA"),
        k_lengths=(9,),
        default_linker_name="(G4S)2")
    assert hasattr(result, 'default_strong_burden')
    assert hasattr(result, 'default_burden')


# ---- review-fix coverage -----------------------------------------------------

def test_score_kmers_warn_no_rank_fires_only_once(caplog):
    """A predictor that returns predictions without percentile_rank
    should produce exactly one warning per process — not one per
    junction × candidate."""
    import logging

    from vaxrank.junction_swap import (
        _reset_score_kmers_warnings,
        _score_kmers,
    )

    class NoRankPredictor:
        def predict_peptides(self, peptides):
            return [SimpleNamespace(peptide=p, allele="HLA-A*02:01")
                    for p in peptides]

    _reset_score_kmers_warnings()
    with caplog.at_level(logging.WARNING):
        for _ in range(5):
            _score_kmers(["KLQGHSAPV"], ["HLA-A*02:01"], NoRankPredictor())
    no_rank_warnings = [r for r in caplog.records
                        if "usable percentile_rank" in r.message]
    assert len(no_rank_warnings) == 1, (
        "Expected exactly one no-rank warning across 5 calls, got %d"
        % len(no_rank_warnings))


def test_score_kmers_warns_when_allele_filter_drops_all(caplog):
    """If the predictor returns only alleles outside the patient set,
    warn once (the optimizer would otherwise silently pick the first
    candidate)."""
    import logging

    from vaxrank.junction_swap import (
        _reset_score_kmers_warnings,
        _score_kmers,
    )

    class WrongAllelePredictor:
        def predict_peptides(self, peptides):
            # Returns only HLA-B*07:02 predictions; caller asks for A*02:01
            return [SimpleNamespace(
                peptide=p, allele="HLA-B*07:02", percentile_rank=10.0)
                for p in peptides]

    _reset_score_kmers_warnings()
    with caplog.at_level(logging.WARNING):
        rows = _score_kmers(
            ["KLQGHSAPV"], ["HLA-A*02:01"], WrongAllelePredictor())
        # Second call must NOT re-emit the warning
        _score_kmers(
            ["KLQGHSAPV"], ["HLA-A*02:01"], WrongAllelePredictor())
    assert rows == [], "All predictions filtered out → no rows"
    filter_warnings = [r for r in caplog.records
                       if "allele filter dropped all" in r.message]
    assert len(filter_warnings) == 1, (
        "Expected exactly one allele-filter warning, got %d"
        % len(filter_warnings))


def test_assemble_manifest_chosen_uses_canonical_linker_names():
    """Manifest's `junction_swap.chosen` field must contain canonical
    Linker.name strings in BOTH branches: the swap-improved branch
    and the default-tied branch. Downstream consumers shouldn't see
    a mix of user-input linker strings and canonical names."""
    from varcode import Variant

    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
    from vaxrank.vaccine_library import get_linker

    fragment_a = SimpleNamespace(
        amino_acids="KLQGH", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=5)
    fragment_b = SimpleNamespace(
        amino_acids="MNNVD", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=5)
    pep_a = SimpleNamespace(
        mutant_protein_fragment=fragment_a, mutant_epitopes=[])
    pep_b = SimpleNamespace(
        mutant_protein_fragment=fragment_b, mutant_epitopes=[])
    pairs = [
        (Variant('1', 100, 'A', 'T'), [pep_a]),
        (Variant('2', 200, 'A', 'T'), [pep_b]),
    ]

    # Default-tied branch: predictor returns no hits → all candidates
    # tie at zero burden → default wins.
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        optimize_linkers=True,
        junction_swap_candidates=("(G4S)2", "AAA"),
        junction_kmer_lengths=(9,),
        antigens_per_construct=2, max_constructs=1,
        max_antigen_length_aa=10,
        utr_3p='HBB',
        # Pass the linker in lowercase so the manifest can't accidentally
        # echo the user's input string. The canonical name from
        # get_linker should be returned instead.
        linker="(g4s)2",
    )
    [c] = assemble_mrna_constructs(
        pairs, options=options,
        mhc_predictor=StubPredictor(rank_table={}),
        mhc_alleles=["HLA-A*02:01"])
    chosen = c.components['junction_swap']['chosen']
    canonical = get_linker("(g4s)2").name
    assert all(name == canonical for name in chosen), (
        "Manifest 'chosen' should contain canonical Linker.name "
        "(%r), got %r" % (canonical, chosen))

    # Swap-improved branch: AAA must beat (G4S)2.
    a1 = "KLQGHSAPVL"
    a2 = "DVIVNCDESLLAS"
    g4s2_kmers = junction_kmers(a1, "GGGGSGGGGS", a2, k_lengths=(9,))
    rank_table = {(g4s2_kmers[0], "HLA-A*02:01"): 0.05}
    pairs2 = [
        (Variant('1', 100, 'A', 'T'),
         [SimpleNamespace(
             mutant_protein_fragment=SimpleNamespace(
                 amino_acids=a1, gene_name='G',
                 mutant_amino_acid_start_offset=0,
                 mutant_amino_acid_end_offset=len(a1)),
             mutant_epitopes=[])]),
        (Variant('2', 200, 'A', 'T'),
         [SimpleNamespace(
             mutant_protein_fragment=SimpleNamespace(
                 amino_acids=a2, gene_name='G',
                 mutant_amino_acid_start_offset=0,
                 mutant_amino_acid_end_offset=len(a2)),
             mutant_epitopes=[])]),
    ]
    options2 = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        optimize_linkers=True,
        junction_swap_candidates=("(G4S)2", "AAA"),
        junction_kmer_lengths=(9,),
        antigens_per_construct=2, max_constructs=1,
        max_antigen_length_aa=20,
        utr_3p='HBB',
        linker="(g4s)2",
    )
    [c2] = assemble_mrna_constructs(
        pairs2, options=options2,
        mhc_predictor=StubPredictor(rank_table=rank_table),
        mhc_alleles=["HLA-A*02:01"])
    chosen2 = c2.components['junction_swap']['chosen']
    # Names must already be canonical (as resolved by get_linker), no
    # raw user strings.
    for name in chosen2:
        assert name == get_linker(name).name, (
            "Manifest 'chosen' entry %r is not the canonical name" % name)


def test_packing_linker_aa_uses_max_candidate_length():
    """The bin-packer must size junctions against the longest possible
    per-junction substitution, not just the shared linker. Otherwise
    a short shared linker + long candidate could overflow the cap at
    swap time."""
    from vaxrank.mrna import RNAConstructConfig, _packing_linker_aa
    from vaxrank.vaccine_library import get_linker

    # Shared linker is short (G2S = 3 aa). Candidate set includes
    # (G4S)2 (10 aa). Optimizer is on. The packer should bill against
    # the 10-aa candidate, not the 3-aa shared linker.
    options = RNAConstructConfig(
        optimize_linkers=True,
        junction_swap_candidates=("G2S", "(G4S)2"),
    )
    shared = get_linker("G2S")
    packing_aa = _packing_linker_aa(options, shared)
    assert len(packing_aa) == 10, (
        "Packer should use max(shared=3, longest_candidate=10) = 10, "
        "got len=%d (%r)" % (len(packing_aa), packing_aa))

    # When optimizer is off, packer uses the shared linker.
    options_off = RNAConstructConfig(
        optimize_linkers=False,
        junction_swap_candidates=("G2S", "(G4S)2"),
    )
    packing_aa_off = _packing_linker_aa(options_off, shared)
    assert packing_aa_off == "GGS", (
        "With optimizer off, packer should use shared linker "
        "as-is, got %r" % packing_aa_off)


def test_packing_uses_longest_candidate_for_size_cap():
    """End-to-end: when optimize_linkers + a long candidate would push
    the construct past max_length_nt, the packer must split *before*
    the optimizer runs (not after, when it's too late)."""
    from varcode import Variant

    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    # Three small antigens. Set max_length_nt tight enough that picking
    # a 10-aa linker (= 30 nt) at every junction overruns, but a 3-aa
    # linker fits — the packer must side with the 10-aa worst case.
    pairs = []
    for k in range(3):
        frag = SimpleNamespace(
            amino_acids="KLQGHSAPVL", gene_name='G%d' % k,
            mutant_amino_acid_start_offset=0,
            mutant_amino_acid_end_offset=10)
        pep = SimpleNamespace(
            mutant_protein_fragment=frag, mutant_epitopes=[])
        pairs.append((Variant(str(k + 1), 100 + k, 'A', 'T'), [pep]))

    # Budget: HBB UTRs (~50 + 132 nt) + start codon 3 nt + 3 antigens
    # × 30 nt + STOP 3 nt + 2 junctions × N nt = ~218 + 60 + 2N.
    # Pick max_length_nt so 2 × 30 nt linker (worst-case) overruns
    # and 2 × 9 nt (shared = G2S) would not.
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        optimize_linkers=True,
        junction_swap_candidates=("G2S", "(G4S)2"),
        linker="G2S",
        antigens_per_construct=10,  # pure length-cap test
        max_constructs=10,
        max_antigen_length_aa=10,
        utr_3p='HBB',
        max_length_nt=295,  # tight: fits with 9-nt linker, not 30-nt
    )
    constructs = assemble_mrna_constructs(
        pairs, options=options,
        mhc_predictor=StubPredictor(rank_table={}),
        mhc_alleles=["HLA-A*02:01"])
    # The conservative packer should split rather than emit a single
    # 3-antigen construct that the optimizer could later overrun.
    assert len(constructs) >= 2, (
        "Conservative packer should split when a long candidate "
        "linker would overrun max_length_nt; got %d construct(s)"
        % len(constructs))
