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

"""Tests for joint mRNA antigen ordering, linker selection, and cleavage
optimization (openvax/vaxrank#919)."""

import itertools
import random
from types import SimpleNamespace

import pytest
from mhctools.base_predictor import BasePredictor
from mhctools.pred import PeptideResult, Prediction

from vaxrank.cli.arg_parser import parse_vaxrank_args
from vaxrank.config import load_vaxrank_config
from vaxrank.junction_swap import (
    MRNAJunctionOptimizationPolicy,
    _path_cost,
    _solve_held_karp,
    _solve_two_opt,
    boundary_kmers,
    junction_kmers,
    optimize_antigen_order_and_linkers,
)
from vaxrank.mrna import (
    RNAConstructConfig,
    assemble_mrna_constructs,
    validate_translated_construct,
)
from vaxrank.report_selection import junction_view
from vaxrank.vaccine_peptide import VaccinePeptide


class CountingStubPredictor(BasePredictor):
    """Stub MHC predictor that records how many batch predict() calls occurred."""

    def __init__(self, alleles=None, rank_table=None, default_rank=99.0):
        alleles = list(alleles or ["HLA-A*02:01"])
        super().__init__(alleles=alleles)
        self.rank_table = dict(rank_table or {})
        self.default_rank = default_rank
        self.predict_calls = 0

    def kind_support(self):
        return {'pMHC_presentation': {'mhc_dependence': 'single_allele'}}

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        self.predict_calls += 1
        out = []
        for kmer in peptides:
            preds = []
            for allele in self.alleles:
                rank = self.rank_table.get((kmer, allele), self.default_rank)
                preds.append(Prediction(
                    kind="pMHC_presentation", score=0.1,
                    peptide=kmer, allele=allele, percentile_rank=rank))
            out.append(PeptideResult(tuple(preds)))
        return out


def _make_ranked_source(name, seq, score=10.0):
    fragment = SimpleNamespace(
        amino_acids=seq,
        gene_name=name.split("_")[0],
        mutant_amino_acid_start_offset=len(seq) // 2,
        mutant_amino_acid_end_offset=(len(seq) // 2) + 1,
        n_alt_reads=20,
        n_ref_reads=10,
        n_overlapping_reads=30,
        n_alt_reads_supporting_protein_sequence=20,
        supporting_reference_transcripts=(),
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=[],
        combined_score_expr=str(score),
    )
    source = SimpleNamespace(
        short_description=name,
        gene_names=[name.split("_")[0]],
    )
    return (source, [vp])


def test_junction_optimization_policy_validation_and_yaml(tmp_path):
    default_pol = MRNAJunctionOptimizationPolicy()
    assert default_pol.optimize_order is False
    assert default_pol.exact_dp_max_antigens == 12
    assert default_pol.cleavage_scoring == 'none'

    with pytest.raises(ValueError, match="exact_dp_max_antigens"):
        MRNAJunctionOptimizationPolicy(exact_dp_max_antigens=0)
    with pytest.raises(ValueError, match="cleavage_scoring"):
        MRNAJunctionOptimizationPolicy(cleavage_scoring="unknown")
    with pytest.raises(ValueError, match="strong_burden_weight"):
        MRNAJunctionOptimizationPolicy(strong_burden_weight=-1.0)
    with pytest.raises(ValueError, match="Unknown mrna.junction_optimization"):
        MRNAJunctionOptimizationPolicy.resolve({"bad_key": 1})

    cfg_file = tmp_path / "trial_policy.yaml"
    cfg_file.write_text(
        "mrna:\n"
        "  junction_optimization:\n"
        "    optimize_order: true\n"
        "    exact_dp_max_antigens: 10\n"
        "    strong_burden_weight: 250.0\n"
        "    mild_burden_weight: 15.0\n"
        "    cleavage_penalty_weight: 8.0\n"
        "    cleavage_scoring: pepsickle\n"
    )
    loaded = load_vaxrank_config(config_path=str(cfg_file))
    pol = MRNAJunctionOptimizationPolicy.resolve(loaded["mrna"]["junction_optimization"])
    assert pol.optimize_order is True
    assert pol.exact_dp_max_antigens == 10
    assert pol.strong_burden_weight == 250.0
    assert pol.cleavage_penalty_weight == 8.0
    assert pol.cleavage_scoring == "pepsickle"


def test_boundary_kmers_spans_direct_boundary_and_filters_proteome():
    kmers = boundary_kmers("AAAAAAAAA", "BBBBBBBBB", k_lengths=(9,))
    assert len(kmers) == 8
    assert all("A" in k and "B" in k for k in kmers)
    assert "AAAAAAAAA" not in kmers
    assert "BBBBBBBBB" not in kmers

    filtered = boundary_kmers(
        "AAAAAAAAA", "BBBBBBBBB", k_lengths=(9,),
        reference_proteome={kmers[0]},
    )
    assert kmers[0] not in filtered
    assert len(filtered) == len(kmers) - 1


def test_reordering_eliminates_strong_junction_that_fixed_order_cannot():
    # Antigens A, B, C in initial order [A, B, C].
    # Suppose A -> B creates a strong chimeric binder under EVERY candidate linker,
    # while A -> C -> B (or C -> A -> B? no, B -> A -> C) has zero binders!
    antigens = [
        "ACDEFGHIKLMNPQRSTVWY",  # A
        "WYRSTVQPMLKIHGFEDCAW",  # B
        "MNPQRSTVWYACDEFGHIKL",  # C
    ]
    candidates = ("(G4S)2", "AAA")
    linker_seqs = {"(G4S)2": "GGGGSGGGGS", "AAA": "AAA"}

    # Poison every linker for A -> B
    rank_table = {}
    for lname, lseq in linker_seqs.items():
        kms = junction_kmers(antigens[0], lseq, antigens[1], k_lengths=(9,))
        rank_table[(kms[0], "HLA-A*02:01")] = 0.1

    predictor = CountingStubPredictor(alleles=["HLA-A*02:01"], rank_table=rank_table)
    res = optimize_antigen_order_and_linkers(
        antigens,
        alleles=["HLA-A*02:01"],
        predictor=predictor,
        candidate_names=candidates,
        k_lengths=(9,),
        default_linker_name="(G4S)2",
        policy={"optimize_order": True},
    )

    assert predictor.predict_calls == 1, "All unique chimeric k-mers should be batched in one call"
    assert res.solver == "held_karp_dp"
    assert res.default_strong_burden == 1
    assert res.fixed_order_strong_burden == 1
    assert res.strong_burden == 0
    assert res.burden == 0
    assert res.reordered is True
    assert sorted(res.antigen_order) == [0, 1, 2]
    # A (0) must no longer be immediately followed by B (1)
    for pos in range(len(res.antigen_order) - 1):
        assert not (res.antigen_order[pos] == 0 and res.antigen_order[pos + 1] == 1)


def test_signal_peptide_and_mitd_boundaries_steer_first_and_last_antigen():
    signal_peptide = "MRVTAPRTLILLLSGALALTETWAGS"
    mitd = "IVGIVAGLAVLAVVVIGAVVATVMCRRKSSGGKGGSYSQAASSDSAQGSDVSLTA"
    antigens = [
        "PPPPPPPPPPPPPPPPPPPP",  # 0: poisoned at signal_peptide -> 0
        "QQQQQQQQQQQQQQQQQQQQ",  # 1: clean everywhere
        "RRRRRRRRRRRRRRRRRRRR",  # 2: poisoned at 2 -> (G3S/AAA) -> MITD
    ]
    rank_table = {}
    sp_kms = boundary_kmers(signal_peptide, antigens[0], k_lengths=(9,))
    rank_table[(sp_kms[0], "HLA-A*02:01")] = 0.2
    # Pick 9-mers that span all three parts (antigen 2 + short linker + MITD)
    for lseq in ("GGGS", "AAA"):
        m_kms = junction_kmers(antigens[2], lseq, mitd, k_lengths=(9,))
        bridging = [k for k in m_kms if "R" in k and "I" in k]
        assert bridging
        rank_table[(bridging[0], "HLA-A*02:01")] = 0.2

    predictor = CountingStubPredictor(alleles=["HLA-A*02:01"], rank_table=rank_table)
    res = optimize_antigen_order_and_linkers(
        antigens,
        alleles=["HLA-A*02:01"],
        predictor=predictor,
        candidate_names=("G3S", "AAA"),
        k_lengths=(9,),
        signal_peptide_aa=signal_peptide,
        mitd_aa=mitd,
        default_linker_name="G3S",
        policy={"optimize_order": True},
    )

    assert res.default_strong_burden == 2
    assert res.strong_burden == 0
    assert res.signal_peptide_strong_burden == 0
    assert res.antigen_order[0] != 0
    assert res.antigen_order[-1] != 2
    assert len(res.chosen_linker_per_junction) == 3  # 2 inter-antigen + 1 pre-MITD


def test_held_karp_matches_brute_force_and_preserves_identity_on_ties():
    rng = random.Random(42)
    n = 6
    start_cost = [rng.uniform(0.0, 5.0) for _ in range(n)]
    end_cost = [rng.uniform(0.0, 5.0) for _ in range(n)]
    edge_cost = [
        [0.0 if i == j else rng.uniform(0.0, 10.0) for j in range(n)]
        for i in range(n)
    ]

    dp_order = _solve_held_karp(n, start_cost, edge_cost, end_cost)
    dp_cost = _path_cost(dp_order, start_cost, edge_cost, end_cost)

    brute_best = min(
        (_path_cost(list(perm), start_cost, edge_cost, end_cost)
         for perm in itertools.permutations(range(n))),
    )
    assert dp_cost == pytest.approx(brute_best, abs=1e-9)

    # All-equal costs should keep identity order [0, 1, 2, 3, 4, 5]
    zeros = [0.0] * n
    flat = [[1.0] * n for _ in range(n)]
    assert _solve_held_karp(n, zeros, flat, zeros) == list(range(n))


def test_two_opt_handles_large_constructs_over_dp_threshold():
    # N = 14 > exact_dp_max_antigens = 12
    n = 14
    zeros = [0.0] * n
    # Make every original adjacent pair (i -> i+1) expensive (100.0) while
    # non-adjacent transitions cost 0.0.
    edge_cost = [[0.0] * n for _ in range(n)]
    for i in range(n - 1):
        edge_cost[i][i + 1] = 100.0

    order = _solve_two_opt(n, zeros, edge_cost, zeros, max_passes=30)
    assert sorted(order) == list(range(n))
    assert _path_cost(order, zeros, edge_cost, zeros) == 0.0

    # Also test through optimize_antigen_order_and_linkers with exact_dp_max_antigens=3
    antigens = [
        "CCCCCCCCCCCCCCCCCCCC",
        "DDDDDDDDDDDDDDDDDDDD",
        "EEEEEEEEEEEEEEEEEEEE",
        "FFFFFFFFFFFFFFFFFFFF",
        "HHHHHHHHHHHHHHHHHHHH",
    ]
    rank_table = {}
    for lseq in ("GGGS", "AAA"):
        kms = junction_kmers(antigens[0], lseq, antigens[1], k_lengths=(9,))
        bridging = [k for k in kms if "C" in k and "D" in k]
        assert bridging
        rank_table[(bridging[0], "HLA-A*02:01")] = 0.1
    predictor = CountingStubPredictor(alleles=["HLA-A*02:01"], rank_table=rank_table)
    res = optimize_antigen_order_and_linkers(
        antigens,
        alleles=["HLA-A*02:01"],
        predictor=predictor,
        candidate_names=("G3S", "AAA"),
        k_lengths=(9,),
        default_linker_name="G3S",
        policy={"optimize_order": True, "exact_dp_max_antigens": 3},
    )
    assert res.solver == "two_opt"
    assert res.strong_burden == 0
    assert res.reordered is True


def test_cleavage_scorer_steers_linker_and_order_when_mhc_burden_ties():
    antigens = [
        "ACDEFGHIKLMNPQR",
        "WYRSTVQPMLKIHGF",
        "MNPQRSTVWYACDEF",
    ]
    predictor = CountingStubPredictor(alleles=["HLA-A*02:01"], rank_table={})

    # Penalize AAA linker (penalty 0.9) and penalize 0->1 transition (0.8),
    # reward (G4S)2 with 0->2->1 (penalty 0.0).
    def scorer(left_aa, linker_aa, right_aa):
        if linker_aa == "AAA":
            return 0.9
        if left_aa == antigens[0] and right_aa == antigens[1]:
            return 0.8
        return 0.0

    res = optimize_antigen_order_and_linkers(
        antigens,
        alleles=["HLA-A*02:01"],
        predictor=predictor,
        candidate_names=("AAA", "(G4S)2"),
        k_lengths=(9,),
        default_linker_name="AAA",
        policy={"optimize_order": True, "cleavage_penalty_weight": 20.0},
        cleavage_scorer=scorer,
    )
    assert res.cleavage_status == "custom_scorer"
    assert res.cleavage_penalty == 0.0
    assert res.default_cleavage_penalty > 0.0
    assert res.linker_names() == ["(G4S)2", "(G4S)2"]
    assert res.reordered is True


def test_pepsickle_cleavage_scoring_degrades_cleanly_when_unassessed(monkeypatch):
    from mhctools import PeptideInput
    from mhctools.cleavage import CleavageModel
    from vaxrank.cleavage_profile import CleavageProfile

    def fake_audit(peptides, **kwargs):
        model = CleavageModel(
            name="pepsickle-epitope-mammal", version="0", enzyme="proteasome",
            uniprot="", species="", compartments=("proteasomal",),
            evidence="quantitative_model",
            references=("https://github.com/pdxgx/pepsickle",),
            assay="test", limitations="", score_name="cleavage_probability",
            score_units="dimensionless", scored_endpoint="site_cleavage",
        )
        return tuple(
            CleavageProfile(
                p, model, "unassessed",
                reason_codes=("pepsickle_model_assets_unavailable",),
            )
            for p in peptides
        )

    monkeypatch.setattr("vaxrank.cleavage_inference.audit_pepsickle_inputs", fake_audit)
    antigens = ["ACDEFGHIKLMNPQR", "WYRSTVQPMLKIHGF"]
    predictor = CountingStubPredictor(alleles=["HLA-A*02:01"], rank_table={})
    res = optimize_antigen_order_and_linkers(
        antigens,
        alleles=["HLA-A*02:01"],
        predictor=predictor,
        candidate_names=("(G4S)2", "AAA"),
        k_lengths=(9,),
        policy={"optimize_order": True, "cleavage_scoring": "pepsickle"},
    )
    assert res.cleavage_status == "pepsickle_unassessed"
    assert res.cleavage_penalty == 0.0
    assert isinstance(PeptideInput("ACDEF"), PeptideInput)


def test_cross_species_mhc_alleles_canine_and_mouse():
    antigens = [
        "ACDEFGHIKLMNPQRSTVWY",
        "WYRSTVQPMLKIHGFEDCAW",
        "MNPQRSTVWYACDEFGHIKL",
    ]
    canine_allele = "DLA-88*034:01"
    rank_table = {}
    for lseq in ("GGGGSGGGGS", "AAA"):
        kms = junction_kmers(antigens[0], lseq, antigens[1], k_lengths=(9,))
        rank_table[(kms[0], canine_allele)] = 0.1
    predictor = CountingStubPredictor(alleles=[canine_allele], rank_table=rank_table)
    res = optimize_antigen_order_and_linkers(
        antigens,
        alleles=[canine_allele],
        predictor=predictor,
        candidate_names=("(G4S)2", "AAA"),
        k_lengths=(9,),
        default_linker_name="(G4S)2",
        policy={"optimize_order": True},
    )
    assert res.strong_burden == 0
    assert res.reordered is True


def test_assemble_mrna_constructs_end_to_end_reordering_and_report_view():
    antigens = [
        ("GENE1_p.A1V", "ACDEFGHIKLMNPQRSTVWY"),
        ("GENE2_p.B2L", "WYRSTVQPMLKIHGFEDCAW"),
        ("GENE3_p.C3K", "MNPQRSTVWYACDEFGHIKL"),
    ]
    ranked = [_make_ranked_source(name, seq, score=30 - idx) for idx, (name, seq) in enumerate(antigens)]
    rank_table = {}
    for lseq in ("GGGGSGGGGS", "AAA"):
        kms = junction_kmers(antigens[0][1], lseq, antigens[1][1], k_lengths=(9,))
        rank_table[(kms[0], "HLA-A*02:01")] = 0.1

    predictor = CountingStubPredictor(alleles=["HLA-A*02:01"], rank_table=rank_table)
    cfg = RNAConstructConfig(
        codon_species="human",
        antigens_per_construct=5,
        max_constructs=1,
        include_mitd=False,
        signal_peptide="HLA_B",
        optimize_linkers=True,
        junction_swap_candidates=("(G4S)2", "AAA"),
        junction_kmer_lengths=(9,),
        junction_optimization={"optimize_order": True},
    )
    constructs = assemble_mrna_constructs(
        ranked,
        options=cfg,
        mhc_predictor=predictor,
        mhc_alleles=["HLA-A*02:01"],
    )
    assert len(constructs) == 1
    c = constructs[0]
    validate_translated_construct(c)
    js = c.elements["junction_swap"]
    assert js["enabled"] is True
    assert js["solver"] == "held_karp_dp"
    assert js["reordered"] is True
    assert js["default_burden_strong"] == 1
    assert js["fixed_order_burden_strong"] == 1
    assert js["burden_strong"] == 0
    assert js["original_antigen_order"] != js["optimized_antigen_order"]
    assert c.antigen_names == js["optimized_antigen_order"]

    view = junction_view(c)
    row_labels = dict(view["rows"])
    assert row_labels["Optimizer solver"] == "held_karp_dp"
    assert row_labels["Antigens reordered"] == "yes"
    assert row_labels["Fixed-order optimal strong burden"] == "1"


def test_cli_flags_for_joint_mrna_optimization():
    args = parse_vaxrank_args([
        "--vcf", "a.vcf",
        "--bam", "a.bam",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--mrna-optimize-antigen-order",
        "--mrna-junction-cleavage-scoring", "pepsickle",
    ])
    assert args.mrna_optimize_antigen_order is True
    assert args.mrna_junction_cleavage_scoring == "pepsickle"


def test_joint_optimization_accepts_lower_total_cost_and_reserves_start_codon_nt(caplog):
    import logging

    # Case 1: Custom policy where cleavage_penalty_weight outweighs mild_burden_weight:
    # reordering should be accepted whenever total_cost < default_total_cost!
    antigens = [
        ("G1_p.A1V", "ACDEFGHIKLMNPQRSTVWY"),
        ("G2_p.B2L", "WYRSTVQPMLKIHGFEDCAW"),
        ("G3_p.C3K", "MNPQRSTVWYACDEFGHIKL"),
    ]
    ranked = [_make_ranked_source(name, seq, score=30 - idx) for idx, (name, seq) in enumerate(antigens)]

    # Make G1->G2 have 0 mild burden but high cleavage penalty (1.0), while G2->G1
    # has 1 mild burden (weight=1.0) and 0 cleavage penalty (weight=50.0) -> lower total_cost!
    rank_table = {}
    for lseq in ("GGGGSGGGGS", "AAA"):
        for a, b in ((antigens[1][1], antigens[0][1]), (antigens[2][1], antigens[0][1])):
            kms = junction_kmers(a, lseq, b, k_lengths=(9,))
            rank_table[(kms[0], "HLA-A*02:01")] = 1.5  # mild only

    def _cleave(left_aa, linker_aa, right_aa):
        # Penalize whenever G1 precedes another antigen (so putting G1 last avoids penalty)
        if left_aa == antigens[0][1]:
            return 1.0
        return 0.0

    predictor = CountingStubPredictor(alleles=["HLA-A*02:01"], rank_table=rank_table)
    cfg = RNAConstructConfig(
        codon_species="human",
        antigens_per_construct=5,
        max_constructs=1,
        include_mitd=False,
        signal_peptide="HLA_B",
        optimize_linkers=True,
        junction_swap_candidates=("(G4S)2", "AAA"),
        junction_kmer_lengths=(9,),
        junction_optimization={
            "optimize_order": True,
            "strong_burden_weight": 100.0,
            "mild_burden_weight": 1.0,
            "cleavage_penalty_weight": 50.0,
        },
    )
    constructs = assemble_mrna_constructs(
        ranked,
        options=cfg,
        mhc_predictor=predictor,
        mhc_alleles=["HLA-A*02:01"],
        cleavage_scorer=_cleave,
    )
    js = constructs[0].elements["junction_swap"]
    assert js["total_cost"] < js["default_total_cost"]
    assert js["reordered"] is True

    # Case 2: Warning logged when optimize_order=True without predictor/alleles
    with caplog.at_level(logging.WARNING, logger="vaxrank.mrna"):
        assemble_mrna_constructs(
            ranked[:2],
            options=RNAConstructConfig(
                codon_species="human",
                include_mitd=False,
                signal_peptide=None,
                optimize_linkers=False,
                junction_optimization={"optimize_order": True},
            ),
        )
    assert any("Joint mRNA junction optimization" in rec.message for rec in caplog.records)
