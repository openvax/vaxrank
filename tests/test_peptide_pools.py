"""Tests for multi-pool SLP formulation partitioning (#917)."""

import csv
import json
from collections import Counter

import pytest
from mhctools import Prediction

from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.cli.arg_parser import make_vaxrank_arg_parser
from vaxrank.config.loader import load_vaxrank_config
from vaxrank.patient_info import PatientInfo
from vaxrank.peptide import (
    PeptideConstructConfig,
    assemble_peptide_constructs,
    write_peptide_outputs,
)
from vaxrank.peptide_pools import (
    net_charge_ph7,
    partition_peptide_constructs,
    summarize_peptide_pools,
)
from vaxrank.report import TemplateDataCreator
from vaxrank.vaccine_antigen import (
    AminoAcidInterval,
    TargetableMask,
    TumorSpecificityAttestation,
    VaccineAntigen,
)
from vaxrank.vaccine_peptide import VaccinePeptide

from .test_product_reports import render


def _make_candidate(
    idx: int,
    sequence: str,
    allele_ranks: dict[str, float],
    *,
    gene: str | None = None,
):
    gene_name = gene or f"GENE{idx + 1}"
    antigen = VaccineAntigen(
        kind="mutation",
        amino_acids=sequence,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status="admitted",
            evidence_kind="synthetic_fixture",
            evidence_source="test",
            patient_specific=True,
        ),
        gene_name=gene_name,
        source_identifier=f"var:{idx + 1}",
    )
    preds = [
        Prediction(
            allele=allele,
            peptide=sequence[:9],
            kind="pMHC_affinity",
            score=1.0,
            value=20.0 if rank <= 0.5 else 200.0,
            percentile_rank=rank,
            predictor_name="netmhcpan",
            predictor_version="4.2c",
        )
        for allele, rank in allele_ranks.items()
    ]
    per_allele_scores = {
        allele: (0.9 if rank <= 0.5 else 0.3)
        for allele, rank in allele_ranks.items()
    }
    ep = CandidateEpitope(
        sequence=sequence[:9],
        source_sequence=sequence,
        offset=0,
        overlaps_targetable=True,
        predictions=preds,
        per_allele_scores=per_allele_scores,
    )
    vp = VaccinePeptide(antigen=antigen, epitopes=[ep])
    return (antigen, [vp])


def test_net_charge_ph7_respects_side_chains_and_termini():
    # Neutral alanine peptide with free N (+1) and free C (-1) -> 0.0
    assert net_charge_ph7("AAAAAAAAAAAAAAA") == 0.0
    # N-terminal acetylation removes +1 -> -1.0
    assert net_charge_ph7("AAAAAAAAAAAAAAA", n_terminal_acetylation=True) == -1.0
    # C-terminal amidation removes -1 -> +1.0
    assert net_charge_ph7("AAAAAAAAAAAAAAA", c_terminal_amidation=True) == 1.0
    # Both capped -> 0.0
    assert net_charge_ph7(
        "AAAAAAAAAAAAAAA", n_terminal_acetylation=True, c_terminal_amidation=True
    ) == 0.0
    # 3 basic (K, R, K) and 1 acidic (D) -> +2.0
    assert net_charge_ph7("AAKRKADAAAAAAAA") == 2.0


def test_multi_pool_separates_same_hla_competitors_to_zero_collisions():
    # 4 peptides binding HLA-A*02:01 + 4 peptides binding HLA-B*07:02.
    # In 1 pool there are C(4,2) + C(4,2) = 12 same-allele pairs.
    # Across 4 pools, optimal partitioning achieves 0 same-allele pairs!
    bases = [
        "ALQSTVNGKPRDFYMA",
        "GLQSTVNGKPRDFYMA",
        "SLQSTVNGKPRDFYMA",
        "TLQSTVNGKPRDFYMA",
        "NLQSTVNGKPRDFYMA",
        "QLQSTVNGKPRDFYMA",
        "PLQSTVNGKPRDFYMA",
        "VLQSTVNGKPRDFYMA",
    ]
    ranked = [
        _make_candidate(i, bases[i], {"HLA-A*02:01": 0.2 if i < 4 else 5.0,
                                      "HLA-B*07:02": 0.2 if i >= 4 else 5.0})
        for i in range(8)
    ]
    audit = {}
    opts = PeptideConstructConfig(max_constructs=8, num_pools=4)
    constructs = assemble_peptide_constructs(
        ranked, opts, target_alleles=["HLA-A*02:01", "HLA-B*07:02"], selection_audit=audit
    )
    assert len(constructs) == 8
    assert constructs[0].pool == 1  # canonical numbering: top construct in Pool 1

    summary = audit["pool_partition"]
    assert summary["num_pools_assigned"] == 4
    assert summary["single_pool_same_allele_pairs"] == 12
    assert summary["partitioned_same_allele_pairs"] == 0
    assert [p["count"] for p in summary["pools"]] == [2, 2, 2, 2]
    for prow in summary["pools"]:
        assert prow["same_allele_pairs"] == 0
        assert prow["allele_counts"] == {"HLA-A*02:01": 1, "HLA-B*07:02": 1}


def test_multi_pool_works_with_non_human_mhc_alleles():
    seqs = [
        "SIINFEKLAAAAAAAA",
        "SIINFEKLGGGGGGGG",
        "ASNENMETMAAAAAAA",
        "ASNENMETMGGGGGGG",
    ]
    ranked = [
        _make_candidate(0, seqs[0], {"H2-Kb": 0.3}),
        _make_candidate(1, seqs[1], {"H2-Kb": 0.4}),
        _make_candidate(2, seqs[2], {"H2-Db": 0.2}),
        _make_candidate(3, seqs[3], {"H2-Db": 0.3}),
    ]
    constructs = assemble_peptide_constructs(
        ranked,
        PeptideConstructConfig(num_pools=2),
        target_alleles=["H2-Kb", "H2-Db"],
    )
    summary = summarize_peptide_pools(constructs)
    assert summary["single_pool_same_allele_pairs"] == 2
    assert summary["partitioned_same_allele_pairs"] == 0


def test_multi_pool_separates_cysteines_and_balances_gravy_and_charge():
    # All 8 peptides have distinct or no allele competition so physicochemical
    # balancing drives the partition: 4 Cys-containing + 4 Cys-free, and
    # alternating hydrophobic/acidic vs hydrophilic/basic.
    seqs = [
        "CIIIIIIIIIIIIIII",  # Cys + very hydrophobic
        "CLLLLLLLLLLLLLLL",  # Cys + very hydrophobic
        "CDDDDDDDDDDDDDDD",  # Cys + acidic + hydrophilic
        "CKKKKKKKKKKKKKKK",  # Cys + basic + hydrophilic
        "AIIIIIIIIIIIIIII",  # no Cys + very hydrophobic
        "ALLLLLLLLLLLLLLL",  # no Cys + very hydrophobic
        "ADDDDDDDDDDDDDDD",  # no Cys + acidic
        "AKKKKKKKKKKKKKKK",  # no Cys + basic
    ]
    ranked = [_make_candidate(i, seqs[i], {}) for i in range(8)]
    constructs = assemble_peptide_constructs(
        ranked, PeptideConstructConfig(num_pools=4)
    )
    summary = summarize_peptide_pools(constructs)
    assert [p["count"] for p in summary["pools"]] == [2, 2, 2, 2]
    # Exactly 1 cysteine-containing peptide per pool -> zero intra-pool Cys pairs!
    assert [p["cys_peptides"] for p in summary["pools"]] == [1, 1, 1, 1]
    # No pool should get two extreme-hydrophobic peptides together.
    gravys = [p["mean_max_gravy"] for p in summary["pools"]]
    assert max(gravys) - min(gravys) < 3.5


def test_cardinality_bounds_edge_cases_and_validation():
    aas = "ACDEFGHIKLMNPQRSTVWY"
    seqs = [f"A{aas[i]}QSTVNGKPRDFYMA" for i in range(10)]
    ranked = [_make_candidate(i, seqs[i], {"HLA-A*01:01": 0.4}) for i in range(10)]

    constructs = assemble_peptide_constructs(
        ranked, PeptideConstructConfig(num_pools=4)
    )
    counts = sorted(Counter(c.pool for c in constructs).values())
    assert counts == [2, 2, 3, 3]

    # Determinism check
    constructs_2 = assemble_peptide_constructs(
        ranked, PeptideConstructConfig(num_pools=4)
    )
    assert [(c.name, c.pool) for c in constructs] == [
        (c.name, c.pool) for c in constructs_2
    ]

    # num_pools > len(constructs) caps cleanly to len(constructs)
    small = assemble_peptide_constructs(
        ranked[:2], PeptideConstructConfig(num_pools=5)
    )
    assert [c.pool for c in small] == [1, 2]

    # Invalid num_pools rejected
    for bad in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="num_pools must be a positive integer"):
            PeptideConstructConfig(num_pools=bad)
        with pytest.raises(ValueError, match="num_pools must be a positive integer"):
            partition_peptide_constructs(small, bad)


def test_outputs_reports_and_cli_yaml_integration(tmp_path):
    seqs = [
        "ALQSTVNGKPRDFYMA",
        "GLQSTVNGKPRDFYMA",
        "SLQSTVNGKPRDFYMA",
        "TLQSTVNGKPRDFYMA",
    ]
    ranked = [
        _make_candidate(i, seqs[i], {"HLA-A*02:01": 0.2 if i < 2 else 5.0,
                                     "HLA-B*07:02": 0.2 if i >= 2 else 5.0})
        for i in range(4)
    ]

    # Default num_pools=1: no pool column in order_form.csv, no pool key in manifest
    single_opts = PeptideConstructConfig(num_pools=1)
    single_constructs = assemble_peptide_constructs(ranked, single_opts)
    assert all(c.pool is None for c in single_constructs)
    write_peptide_outputs(
        single_constructs,
        fasta_path=str(tmp_path / "single.fasta"),
        manifest_path=str(tmp_path / "single.json"),
        order_form_path=str(tmp_path / "single.csv"),
        options=single_opts,
    )
    with open(tmp_path / "single.csv") as f:
        single_header = next(csv.reader(f))
    assert "pool" not in single_header
    single_manifest = json.loads((tmp_path / "single.json").read_text())
    assert all("pool" not in entry for entry in single_manifest)

    # Multi-pool num_pools=2: FASTA, manifest, order_form, HTML & ASCII reports
    multi_opts = PeptideConstructConfig(num_pools=2)
    multi_constructs = assemble_peptide_constructs(
        ranked, multi_opts, target_alleles=["HLA-A*02:01", "HLA-B*07:02"]
    )
    write_peptide_outputs(
        multi_constructs,
        fasta_path=str(tmp_path / "multi.fasta"),
        manifest_path=str(tmp_path / "multi.json"),
        order_form_path=str(tmp_path / "multi.csv"),
        options=multi_opts,
    )
    fasta_text = (tmp_path / "multi.fasta").read_text()
    assert "pool=1 " in fasta_text and "pool=2 " in fasta_text

    multi_manifest = json.loads((tmp_path / "multi.json").read_text())
    assert {entry["pool"] for entry in multi_manifest} == {1, 2}

    with open(tmp_path / "multi.csv") as f:
        rows = list(csv.DictReader(f))
    assert list(rows[0].keys())[:2] == ["name", "pool"]
    assert {int(r["pool"]) for r in rows} == {1, 2}

    # Reports
    data = TemplateDataCreator(
        ranked, PatientInfo("pool-test"), final_review="", reviewers="",
        args_for_report={"manufacturability": True, "wt_epitopes": False},
        input_json_file=None, products={"peptide": multi_constructs},
    ).compute_template_data()
    html = render(data, "html")
    ascii_txt = render(data, "txt")
    assert "Formulation pools (2)" in html
    assert "Pool 1</span>" in html and "Pool 2</span>" in html
    assert "from 2 in a single pool to 0" in html
    assert "Formulation pools: 2 (same-allele binder pairs: 2 single-pool -> 0 partitioned)." in ascii_txt
    assert "[Pool 1]" in ascii_txt and "[Pool 2]" in ascii_txt

    # CLI & YAML config round-trip
    parser = make_vaxrank_arg_parser()
    parsed = parser.parse_args(["--mhc-predictor", "random", "--peptide-num-pools", "4"])
    assert parsed.peptide_num_pools == 4

    cfg_path = tmp_path / "pools.yaml"
    cfg_path.write_text("peptide:\n  num_pools: 4\n")
    merged = load_vaxrank_config(config_path=str(cfg_path))
    assert merged["peptide"]["num_pools"] == 4


def test_cli_end_to_end_and_offline_replay_preserves_pools(tmp_path):
    from vaxrank.cli.entry_point import run_cli
    from vaxrank.epitope_dataset import ANTIGEN_COLUMN
    from vaxrank.native_serialization import to_native_json
    from .test_epitope_dataset import evidence

    sequences = ["AASIINFEKLKKKKKKK", "AAGILGFVFTLKKKKKK"]
    antigens = [
        VaccineAntigen(
            kind="viral",
            amino_acids=seq,
            gene_name=f"SOURCE_{i}",
            source_identifier=f"source-{i}",
            targetable_mask=TargetableMask((AminoAcidInterval(0, len(seq)),)),
            tumor_specificity=TumorSpecificityAttestation(
                status="admitted",
                evidence_kind="synthetic_fixture",
                evidence_source="test_peptide_pools",
                patient_specific=True,
                rationale_code="test_only",
            ),
        )
        for i, seq in enumerate(sequences)
    ]
    source = tmp_path / "source.tsv"
    evidence(**{
        ANTIGEN_COLUMN: [to_native_json(a) for a in antigens],
        "source_sequence": sequences,
        "peptide_offset": [2, 2],
        "host_species": ["human", "human"],
        "sample_name": ["patient", "patient"],
    }).to_tsv(source)

    out1 = tmp_path / "run1"
    native = tmp_path / "native.tsv"
    run_cli([
        "--input-topiary", str(source),
        "--vaccine-type", "peptide",
        "--peptide-num-pools", "2",
        "--output-dir", str(out1),
        "--output-epitopes", str(native),
        "--pdf-backend", "weasyprint",
        "--no-processing-aware-annotation",
        "--config-text", "epitopes.score_expr=1.0",
    ])
    with open(out1 / "order_form.csv") as f:
        rows1 = list(csv.DictReader(f))
    assert {int(r["pool"]) for r in rows1} == {1, 2}

    # Replay preserves exact pool assignments byte-for-byte
    out2 = tmp_path / "replay"
    run_cli([
        "--input-epitopes", str(native),
        "--output-dir", str(out2),
        "--pdf-backend", "weasyprint",
        "--no-processing-aware-annotation",
    ])
    for fname in ("order_form.csv", "manifest.json", "vaccine.fasta"):
        assert (out1 / fname).read_bytes() == (out2 / fname).read_bytes()

    # Overriding --peptide-num-pools 1 on replay drops the pool partition cleanly
    out3 = tmp_path / "replay_single"
    run_cli([
        "--input-epitopes", str(native),
        "--peptide-num-pools", "1",
        "--output-dir", str(out3),
        "--pdf-backend", "weasyprint",
        "--no-processing-aware-annotation",
    ])
    with open(out3 / "order_form.csv") as f:
        header3 = next(csv.reader(f))
    assert "pool" not in header3

