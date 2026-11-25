"""Tests for synthesis-safe mRNA CDS optimization, QC, and LinearDesign/MFE hooks (#918)."""

import json
import stat
import sys
from types import SimpleNamespace

from Bio.Seq import Seq
import pytest

from vaxrank.cli.arg_parser import make_vaxrank_arg_parser
from vaxrank.config.loader import load_vaxrank_config
from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs, codon_optimize
from vaxrank.mrna_cds import (
    MRNACDSOptimizationPolicy,
    T7_CLASS_II_TERMINATOR,
    compute_mrna_sequence_qc,
    optimize_cds_dna,
    resolve_restriction_motif,
)
from vaxrank.patient_info import PatientInfo
from vaxrank.report import TemplateDataCreator
from vaxrank.vaccine_antigen import (
    AminoAcidInterval,
    TargetableMask,
    TumorSpecificityAttestation,
    VaccineAntigen,
)
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.vaccine_peptide import VaccinePeptide

from .test_product_reports import render


def _make_ranked(n=3):
    seqs = [
        "FLSFMLRAVSFFVKDAVLYSGATLD",
        "LKHGKRFHATISFDTDTGLKQALET",
        "PIRGIRMKFVVLTGLVEVGEVSNRD",
    ]
    out = []
    for i in range(n):
        seq = seqs[i % len(seqs)]
        ag = VaccineAntigen(
            kind="mutation",
            amino_acids=seq,
            targetable_mask=TargetableMask((AminoAcidInterval(10, 11),)),
            tumor_specificity=TumorSpecificityAttestation(
                status="admitted",
                evidence_kind="synthetic_fixture",
                evidence_source="test_mrna_cds",
                patient_specific=True,
            ),
            gene_name=f"GENE{i + 1}",
            source_identifier=f"var:{i + 1}",
        )
        ep = CandidateEpitope(
            sequence=seq[6:15],
            source_sequence=seq,
            offset=6,
            overlaps_targetable=True,
            per_allele_scores={"HLA-A*02:01": 0.9 - 0.1 * i},
        )
        out.append((ag, [VaccinePeptide(antigen=ag, epitopes=[ep])]))
    return out


def test_synthesis_safe_cds_enforces_gc_homopolymer_t7_restriction_and_u_depletion():
    # Protein containing motifs that can encode T7 terminator (I-C-S -> ATC TGT)
    # and BsaI (G-L -> GGT CTC) plus leucine/serine runs.
    aa = "MFLSFMLRAVSFFVKDAVLYSGATLDICSLKHGKRFHATISFDTDTGLKQALETGLPIRGIRMKFVVLTGLVEVGEVSNRD"
    policy = MRNACDSOptimizationPolicy.synthesis_safe(
        avoid_restriction_sites=("BsaI", "SapI"),
        uridine_depletion=True,
        uridine_depletion_weight=2.0,
    )
    dna = codon_optimize(aa, species="h_sapiens", cds_optimization=policy)
    assert str(Seq(dna).translate()) == aa

    qc = compute_mrna_sequence_qc(dna + "TGA", species="h_sapiens", policy=policy)
    assert qc["min_window_gc"] >= 0.35
    assert qc["max_window_gc"] <= 0.68
    assert qc["max_homopolymer_length"] <= 5
    assert qc["t7_terminator_count"] == 0
    assert T7_CLASS_II_TERMINATOR not in dna
    assert resolve_restriction_motif("BsaI") not in dna
    assert resolve_restriction_motif("SapI") not in dna
    assert qc["synthesis_flags"] == []
    assert qc["cai"] is not None and 0.0 < qc["cai"] <= 1.0


def test_impossible_gc_run_falls_back_gracefully_and_flags_qc():
    # Poly-glycine (GGN) has minimum 66.7% GC; forcing gc_max=0.45 is impossible
    # without breaking translation, so optimizer must preserve translation and flag QC.
    aa = "M" + ("G" * 25)
    policy = MRNACDSOptimizationPolicy(gc_min=0.30, gc_max=0.45, gc_window_nt=30)
    dna = codon_optimize(aa, species="h_sapiens", cds_optimization=policy)
    assert str(Seq(dna).translate()) == aa
    qc = compute_mrna_sequence_qc(dna + "TGA", species="h_sapiens", policy=policy)
    assert any("max_30nt_gc" in flag for flag in qc["synthesis_flags"])


def test_lineardesign_hook_with_frozen_segment_and_error_checks(tmp_path):
    # Create a tiny fake LinearDesign CLI script that outputs a valid RNA sequence.
    fake_ld = tmp_path / "fake_lineardesign"
    fake_ld.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "from dnachisel.biotools import reverse_translate\n"
        "aa = sys.stdin.read().strip()\n"
        "rna = reverse_translate(aa).replace('T', 'U')\n"
        "print('mRNA sequence:  ' + rna)\n"
        "print('mRNA structure: ' + '.' * len(rna))\n"
        "print('mRNA folding free energy: -42.50 kcal/mol; mRNA CAI: 0.910')\n"
    )
    fake_ld.chmod(fake_ld.stat().st_mode | stat.S_IXUSR)

    aa = "MFLSFMLRAVSFFVKDAVLY"
    # Pin middle 3 residues (offset 6..9 = 'LRA') to a specific blessed DNA
    blessed = "TTACGCGCG"
    assert str(Seq(blessed).translate()) == aa[6:9]
    policy = MRNACDSOptimizationPolicy(
        optimizer="lineardesign",
        lineardesign_executable=str(fake_ld),
        lineardesign_lambda=0.5,
    )
    dna = codon_optimize(
        aa, species="h_sapiens", frozen_segments=[(6, 9, blessed)], cds_optimization=policy
    )
    assert str(Seq(dna).translate()) == aa
    assert dna[18:27] == blessed

    # Missing binary raises a clear RuntimeError
    missing_policy = MRNACDSOptimizationPolicy(
        optimizer="lineardesign",
        lineardesign_executable=str(tmp_path / "does_not_exist"),
    )
    with pytest.raises(RuntimeError, match="not found on PATH"):
        codon_optimize(aa, species="h_sapiens", cds_optimization=missing_policy)


def test_optional_viennarna_mfe_hook_and_graceful_absence(monkeypatch):
    cds = "ATGTTCCTGAGCTTCATGCTGAGAGCCGTGAGCTTCTTCGTGAAGGACGCCGTGCTGTGA"

    # Simulate ViennaRNA installed
    fake_rna = SimpleNamespace(fold=lambda seq: ("." * len(seq), -18.4))
    monkeypatch.setitem(sys.modules, "RNA", fake_rna)
    qc = compute_mrna_sequence_qc(
        cds,
        species="h_sapiens",
        policy=MRNACDSOptimizationPolicy(compute_mfe=True, initiation_window_nt=30),
    )
    assert qc["mfe_status"] == "viennarna"
    assert qc["cds_mfe_kcal_mol"] == -18.4
    assert qc["initiation_window_mfe_kcal_mol"] == -18.4

    # Simulate ViennaRNA absent
    monkeypatch.setitem(sys.modules, "RNA", None)
    qc_missing = compute_mrna_sequence_qc(
        cds,
        species="h_sapiens",
        policy=MRNACDSOptimizationPolicy(compute_mfe=True),
    )
    assert qc_missing["mfe_status"] == "viennarna_not_installed"
    assert qc_missing["cds_mfe_kcal_mol"] is None


def test_mrna_assembly_reports_yaml_policy_and_cli_integration(tmp_path):
    trial_yaml = tmp_path / "trial_mrna.yaml"
    trial_yaml.write_text(
        "mrna:\n"
        "  codon_species: h_sapiens\n"
        "  cds_optimization:\n"
        "    gc_min: 0.35\n"
        "    gc_max: 0.68\n"
        "    gc_window_nt: 50\n"
        "    max_homopolymer_length: 5\n"
        "    avoid_t7_terminator: true\n"
        "    avoid_restriction_sites:\n"
        "      - BsaI\n"
        "      - SapI\n"
        "    uridine_depletion: true\n"
        "    uridine_depletion_weight: 1.5\n"
    )
    merged = load_vaxrank_config(config_path=["builtin:openvax-v1", str(trial_yaml)])
    cds_cfg = merged["mrna"]["cds_optimization"]
    assert cds_cfg["gc_min"] == 0.35
    assert cds_cfg["avoid_restriction_sites"] == ["BsaI", "SapI"]

    ranked = _make_ranked(3)
    opts = RNAConstructConfig(
        codon_species="h_sapiens",
        optimize_linkers=False,
        cds_optimization=cds_cfg,
    )
    constructs = assemble_mrna_constructs(ranked, options=opts)
    assert len(constructs) == 1
    qc = constructs[0].elements["sequence_qc"]
    assert qc["cds_length_nt"] == len(constructs[0].cds_nt)
    assert qc["t7_terminator_count"] == 0
    assert constructs[0].elements["cds_optimization"]["uridine_depletion"] is True

    # Reports include CDS sequence QC
    data = TemplateDataCreator(
        ranked, PatientInfo("mrna-qc-test"), final_review="", reviewers="",
        args_for_report={"manufacturability": False, "wt_epitopes": False},
        input_json_file=None, products={"mrna": constructs},
    ).compute_template_data()
    html = render(data, "html")
    txt = render(data, "txt")
    assert "CDS sequence QC:" in html
    assert "Max homopolymer:" in html
    assert "CDS sequence QC: GC " in txt

    # CLI flag parsing
    parser = make_vaxrank_arg_parser()
    parsed = parser.parse_args([
        "--mhc-predictor", "random",
        "--mrna-cds-synthesis-safe",
        "--mrna-codon-optimizer", "lineardesign",
    ])
    assert parsed.mrna_cds_synthesis_safe is True
    assert parsed.mrna_codon_optimizer == "lineardesign"

    # Validation checks
    with pytest.raises(ValueError, match="Unknown mrna.cds_optimization keys"):
        MRNACDSOptimizationPolicy.resolve({"bad_key": True})
    with pytest.raises(ValueError, match="gc_min"):
        MRNACDSOptimizationPolicy(gc_min=0.7, gc_max=0.4)
    with pytest.raises(ValueError, match="Unknown restriction enzyme"):
        MRNACDSOptimizationPolicy(avoid_restriction_sites=("NotAnEnzyme!",))

    # JSON roundtrip of sequence_qc & cds_optimization
    json.dumps(constructs[0].elements)


def test_uridine_depletion_only_penalizes_sense_strand_u_not_a():
    # Poly-lysine (AAA vs AAG) has zero sense-strand T/U codons. Before adding
    # strand=1, AvoidPattern("T") matched antisense T (sense A) and forced AAG.
    # With strand=1 and uridine depletion on, a poly-phenylalanine (TTT vs TTC)
    # MUST prefer TTC (fewer sense-strand T), while lysine AAA is NOT penalized
    # as a sense-strand U.
    aa_f = "F" * 12
    dna_f = optimize_cds_dna(
        aa_f,
        species="h_sapiens",
        method="match_codon_usage",
        policy=MRNACDSOptimizationPolicy(
            uridine_depletion=True,
            uridine_depletion_weight=5.0,
            avoid_uu_dinucleotides=True,
        ),
    )
    assert dna_f.count("TTC") == 12
    assert dna_f.count("TTT") == 0

    # Verify stop-codon boundary prevents terminal homopolymer / PacI across CDS+TAA
    # Terminal IN* with ATTAAT + TAA = ATTAATTAA (PacI = TTAATTAA!)
    aa_paci = "M" + "G" * 8 + "IN"
    dna_paci = optimize_cds_dna(
        aa_paci,
        species="h_sapiens",
        policy=MRNACDSOptimizationPolicy(
            avoid_restriction_sites=("PacI",),
            max_homopolymer_length=4,
        ),
        stop_codon="TAA",
    )
    assert "TTAATTAA" not in (dna_paci + "TAA")
    qc = compute_mrna_sequence_qc(
        dna_paci + "TAA",
        species="h_sapiens",
        policy=MRNACDSOptimizationPolicy(
            avoid_restriction_sites=("PacI",),
            max_homopolymer_length=4,
        ),
    )
    assert qc["synthesis_flags"] == []
    assert qc["max_homopolymer_length"] <= 4


def test_mrna_cds_synthesis_safe_cli_overrides_default_yaml_false():
    from vaxrank.mrna_cds import MRNACDSOptimizationPolicy

    # Simulate default.yaml populating avoid_t7_terminator=False, uridine_depletion=False
    merged = load_vaxrank_config(config_path="builtin:default")
    yaml_cds = dict(merged["mrna"]["cds_optimization"])
    assert yaml_cds["avoid_t7_terminator"] is False
    assert yaml_cds["uridine_depletion"] is False

    baseline = MRNACDSOptimizationPolicy().to_dict()
    preset = MRNACDSOptimizationPolicy.synthesis_safe().to_dict()
    for k, v in preset.items():
        cur = yaml_cds.get(k)
        if k not in yaml_cds or cur is None or cur == baseline.get(k):
            yaml_cds[k] = v
    resolved = MRNACDSOptimizationPolicy.resolve(yaml_cds)
    assert resolved.avoid_t7_terminator is True
    assert resolved.uridine_depletion is True
    assert resolved.max_homopolymer_length == 5


def test_lineardesign_warns_on_non_human_species_and_active_constraints(monkeypatch, caplog):
    import logging

    monkeypatch.setattr("vaxrank.mrna_cds.shutil.which", lambda _: "/usr/local/bin/lineardesign")
    monkeypatch.setattr(
        "vaxrank.mrna_cds.subprocess.run",
        lambda *a, **kw: SimpleNamespace(returncode=0, stdout="mRNA sequence: ATGTTC\n", stderr=""),
    )
    with caplog.at_level(logging.WARNING, logger="vaxrank.mrna_cds"):
        dna = optimize_cds_dna(
            "MF",
            species="m_musculus",
            policy=MRNACDSOptimizationPolicy(
                optimizer="lineardesign",
                avoid_t7_terminator=True,
            ),
        )
    assert dna == "ATGTTC"
    assert any("ignores requested codon_species" in rec.message for rec in caplog.records)
    assert any("does not enforce DnaChisel synthesis" in rec.message for rec in caplog.records)
