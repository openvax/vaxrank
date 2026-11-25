"""Tests for cross-species bulk RNA MHC ASE (upstream tools + direct k-mer +
isovar-assisted) and scRNA tumor-cluster nTPM + MHC ASE (#936)."""

import gzip
import json
from types import SimpleNamespace

import pytest

from vaxrank.cli.entry_point import run_cli
from vaxrank.config import load_vaxrank_config
from vaxrank.mhc_allele_status import (
    load_mhc_allele_status,
    merge_mhc_allele_status_bundles,
    resolve_effective_status,
)
from vaxrank.mhc_expression import (
    MHCExpressionPolicy,
    estimate_bulk_mhc_ase_from_isovar,
    estimate_bulk_mhc_ase_from_reads,
    load_mhc_status_from_tool,
)
from vaxrank.scrna_tumor import (
    SCRNATumorClusterPolicy,
    analyze_scrna_tumor_clusters,
    load_cluster_assignments,
    load_scrna_counts,
    load_scrna_mhc_allele_umis,
    load_somatic_barcode_calls,
)


def test_upstream_companion_tools_human_canine_and_merging(tmp_path):
    # 1. T1K headerless genotype file (human HLA)
    t1k_human = tmp_path / "patient_genotype.tsv"
    t1k_human.write_text(
        "HLA-A\t2\tHLA-A*02:01\t2.0\t40\tHLA-A*24:02\t98.0\t60\n"
        "HLA-B\t2\tHLA-B*07:02\t50.0\t55\tHLA-B*44:02\t50.0\t55\n"
    )
    b_t1k = load_mhc_status_from_tool(t1k_human, "t1k", tumor_purity=0.75)
    by_a = b_t1k.by_allele()
    assert by_a["HLA-A*02:01"].bulk_rna_allele_ratio == pytest.approx(0.02)
    assert by_a["HLA-A*02:01"].dna_loh is None  # strictly RNA axis
    assert resolve_effective_status(by_a["HLA-A*02:01"]).rna_status == "silenced"
    assert resolve_effective_status(by_a["HLA-A*02:01"]).dna_status == "unknown"
    assert resolve_effective_status(by_a["HLA-A*24:02"]).effective_status == "retained"

    # 2. T1K headered file (canine DLA)
    t1k_dog = tmp_path / "dog_t1k.tsv"
    t1k_dog.write_text(
        "allele\ttpm\n"
        "DLA-88*034:01\t3.0\n"
        "DLA-88*508:01\t97.0\n"
    )
    b_dog = load_mhc_status_from_tool(t1k_dog, "t1k", tumor_purity=0.80)
    assert resolve_effective_status(b_dog.by_allele()["DLA-88*034:01"]).effective_status == "silenced"

    # 3. LOHHLA DNA LOH file merged with T1K RNA -> dna_and_rna mechanism!
    lohhla = tmp_path / "patient.HLAlossPrediction_CI.tsv"
    lohhla.write_text(
        "HLA_A_type1\tHLA_A_type2\tHLA_type1copyNum_withBAFBin\tHLA_type2copyNum_withBAFBin\tPVal_unique\n"
        "hla_a_02_01_01_01\thla_a_24_02_01_01\t0.05\t1.95\t0.0001\n"
    )
    b_loh = load_mhc_status_from_tool(lohhla, "lohhla")
    merged = merge_mhc_allele_status_bundles(b_loh, b_t1k)
    assert merged is not None
    dec = resolve_effective_status(merged.by_allele()["HLA-A*02:01"])
    assert dec.dna_status == "lost"
    assert dec.rna_status == "silenced"
    assert dec.loss_mechanism == "dna_and_rna"
    assert dec.masked is True

    # 4. arcasHLA, HLApers, DASH, LILAC, POLYSOLVER
    arcas = tmp_path / "arcas.quant.tsv"
    arcas.write_text("allele\ttpm\nHLA-A*02:01\t5.0\nHLA-A*01:01\t95.0\n")
    assert load_mhc_status_from_tool(arcas, "arcashla").by_allele()["HLA-A*02:01"].bulk_rna_allele_ratio == pytest.approx(0.05)

    hlapers = tmp_path / "hlapers_quant.tsv"
    hlapers.write_text("target_id\test_counts\ttpm\nHLA-B*07:02\t10\t4.0\nHLA-B*08:01\t190\t96.0\n")
    assert load_mhc_status_from_tool(hlapers, "hlapers").by_allele()["HLA-B*07:02"].bulk_rna_allele_ratio == pytest.approx(0.04)

    dash = tmp_path / "dash.csv"
    dash.write_text("allele,predicted_deletion,deletion_probability,copy_number\nHLA-A*02:01,True,0.98,0.1\n")
    assert load_mhc_status_from_tool(dash, "dash").by_allele()["HLA-A*02:01"].dna_loh is True

    lilac = tmp_path / "sample.lilac.tsv"
    lilac.write_text(
        "Allele\tTumorCopyNumber\tSomaticMissense\tSomaticNonsenseOrFrameshift\tSomaticSplice\tSomaticInframeIndel\n"
        "HLA-A*02:01\t1.0\t0\t1\t0\t0\n"
        "HLA-B*07:02\t0.1\t0\t0\t0\t0\n"
    )
    b_lilac = load_mhc_status_from_tool(lilac, "lilac")
    assert resolve_effective_status(b_lilac.by_allele()["HLA-A*02:01"]).dna_status == "truncating_lof"
    assert resolve_effective_status(b_lilac.by_allele()["HLA-B*07:02"]).dna_status == "lost"

    poly = tmp_path / "polysolver.tsv"
    poly.write_text("allele\tconsequence\thgvs_p\nhla_c_07_02_01_01\tnonsense\tp.Q88*\n")
    b_poly = load_mhc_status_from_tool(poly, "polysolver")
    assert resolve_effective_status(b_poly.by_allele()["HLA-C*07:02"]).dna_status == "truncating_lof"


def test_direct_kmer_bulk_mhc_ase_estimator_canine_and_purity_guard(tmp_path):
    # Personalized canine DLA-88 alleles with distinct polymorphic exon-2 regions
    shared_prefix = "ATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCG"
    poly_a = "GCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGC"
    poly_b = "CGTAGGTCAGCTTGCAGTCCCGGAATTCCGATAGCGTTAACGGCTAGCTTACCGGTAGCCAT"
    shared_suffix = "TTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGC"

    fasta = tmp_path / "dog_mhc.fa"
    fasta.write_text(
        f">DLA-88*034:01\n{shared_prefix}{poly_a}{shared_suffix}\n"
        f">DLA-88*508:01\n{shared_prefix}{poly_b}{shared_suffix}\n"
    )

    # Synthetic bulk RNA-seq FASTQ.gz: 2 reads from DLA-88*034:01, 48 reads from DLA-88*508:01
    fq = tmp_path / "tumor_rna.fastq.gz"
    with gzip.open(fq, "wt", encoding="utf-8") as fh:
        for idx in range(2):
            fh.write(f"@r_a_{idx}\n{poly_a}\n+\n{'I' * len(poly_a)}\n")
        for idx in range(48):
            fh.write(f"@r_b_{idx}\n{poly_b}\n+\n{'I' * len(poly_b)}\n")

    # High purity (0.80) -> DLA-88*034:01 cleanly resolved as silenced!
    b_high = estimate_bulk_mhc_ase_from_reads(fasta, fq, tumor_purity=0.80)
    dec_a = resolve_effective_status(b_high.by_allele()["DLA-88*034:01"])
    dec_b = resolve_effective_status(b_high.by_allele()["DLA-88*508:01"])
    assert dec_a.effective_status == "silenced"
    assert dec_a.masked is True
    assert dec_b.effective_status == "retained"
    assert dec_b.masked is False

    # Low purity (0.20) -> purity guard prevents false silencing call (unresolvable!)
    b_low = estimate_bulk_mhc_ase_from_reads(fasta, fq, tumor_purity=0.20)
    dec_low = resolve_effective_status(b_low.by_allele()["DLA-88*034:01"])
    assert dec_low.effective_status == "unresolvable"
    assert dec_low.masked is False
    assert "bulk_rna_purity_confounded" in dec_low.caveats


def test_isovar_assisted_mhc_local_assembly_ase():
    var_map = {
        "12:100500:A>G": {
            "ref_allele": "DLA-88*034:01",
            "alt_allele": "DLA-88*508:01",
        }
    }
    iso_res = [
        SimpleNamespace(
            variant=SimpleNamespace(short_description="12:100500:A>G"),
            num_ref_supporting_reads=2,
            num_alt_supporting_reads=58,
        )
    ]
    bundle = estimate_bulk_mhc_ase_from_isovar(var_map, iso_res, tumor_purity=0.75)
    by_a = bundle.by_allele()
    assert by_a["DLA-88*034:01"].bulk_rna_allele_ratio == pytest.approx(2 / 60, abs=1e-4)
    assert by_a["DLA-88*508:01"].bulk_rna_allele_ratio == pytest.approx(58 / 60, abs=1e-4)
    assert resolve_effective_status(by_a["DLA-88*034:01"]).effective_status == "silenced"
    assert resolve_effective_status(by_a["DLA-88*508:01"]).effective_status == "retained"


def test_scrna_tumor_cluster_extraction_ntpm_apm_and_5p_vs_3p_mhc_ase(tmp_path):
    # Build synthetic 10x MTX directory with 3 clusters (6 cells each):
    # C0 = tumor (EPCAM high, PTPRC=0, B2M silenced, somatic variants present)
    # C1 = immune control (PTPRC/CD3D/CD8A high, B2M high)
    # C2 = stromal CAF (COL1A1/DCN high, PTPRC=0)
    mtx_dir = tmp_path / "filtered_feature_bc_matrix"
    mtx_dir.mkdir()
    genes = ["EPCAM", "PTPRC", "CD3D", "COL1A1", "DCN", "B2M", "TAP1", "TP53", "GAPDH"]
    barcodes = [f"TUM_{i}" for i in range(6)] + [f"IMM_{i}" for i in range(6)] + [f"STR_{i}" for i in range(6)]
    (mtx_dir / "features.tsv").write_text("\n".join(f"G{i}\t{g}\tGene Expression" for i, g in enumerate(genes)) + "\n")
    (mtx_dir / "barcodes.tsv").write_text("\n".join(barcodes) + "\n")

    entries = []
    # MatrixMarket is 1-based (row=gene, col=barcode, val)
    for col_idx, bc in enumerate(barcodes, start=1):
        if bc.startswith("TUM_"):
            row_vals = {1: 400, 6: 1, 7: 30, 8: 250, 9: 319}  # EPCAM, tiny B2M, TAP1, TP53, GAPDH
        elif bc.startswith("IMM_"):
            row_vals = {2: 300, 3: 200, 6: 150, 7: 40, 8: 20, 9: 290}  # PTPRC, CD3D, B2M high
        else:
            row_vals = {4: 450, 5: 250, 6: 80, 9: 220}  # COL1A1, DCN stromal
        for r, v in row_vals.items():
            entries.append(f"{r} {col_idx} {v}")
    (mtx_dir / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n"
        f"{len(genes)} {len(barcodes)} {len(entries)}\n"
        + "\n".join(entries) + "\n"
    )

    clusters_tsv = tmp_path / "clusters.tsv"
    clusters_tsv.write_text(
        "barcode\tcluster\n"
        + "\n".join(
            f"{bc}\t{'C0' if bc.startswith('TUM_') else 'C1' if bc.startswith('IMM_') else 'C2'}"
            for bc in barcodes
        ) + "\n"
    )
    somatic_tsv = tmp_path / "vartrix.tsv"
    somatic_tsv.write_text(
        "barcode\talt_count\tref_count\n"
        + "\n".join(f"TUM_{i}\t4\t2" for i in range(5)) + "\n"
    )

    # scHLAcount-style directory: 10x 5' resolvable UMIs (HLA-A*02:01 silenced in C0, intact in C1)
    schla_dir = tmp_path / "schlacount_5p"
    schla_dir.mkdir()
    (schla_dir / "alleles.tsv").write_text("HLA-A*02:01\nHLA-B*07:02\n")
    (schla_dir / "barcodes.tsv").write_text("\n".join(barcodes) + "\n")
    mhc_entries = []
    for col_idx, bc in enumerate(barcodes, start=1):
        if bc.startswith("TUM_"):
            mhc_entries.append(f"2 {col_idx} 25")  # HLA-B*07:02 expressed; HLA-A*02:01 = 0
        elif bc.startswith("IMM_"):
            mhc_entries.append(f"1 {col_idx} 20")  # HLA-A*02:01 intact in immune control (120 UMIs)
            mhc_entries.append(f"2 {col_idx} 20")  # HLA-B*07:02 intact in immune control
    (schla_dir / "count_matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n"
        f"2 {len(barcodes)} {len(mhc_entries)}\n"
        + "\n".join(mhc_entries) + "\n"
    )

    counts = load_scrna_counts(mtx_dir)
    clusters = load_cluster_assignments(clusters_tsv)
    somatic = load_somatic_barcode_calls(somatic_tsv)
    mhc_umis = load_scrna_mhc_allele_umis(schla_dir)

    res_5p = analyze_scrna_tumor_clusters(
        counts, clusters,
        somatic_barcodes=somatic,
        mhc_allele_umis=mhc_umis,
        policy={"chemistry": "10x_5p", "tumor_markers": ["EPCAM"]},
    )
    assert res_5p.tumor_cluster_ids == ("C0",)
    assert res_5p.immune_control_cluster_ids == ("C1",)
    # Tumor-cluster TP53 nTPM is high (~250,000 nTPM) vs immune control (~20,000)
    assert res_5p.tumor_gene_ntpm["TP53"] > res_5p.immune_control_gene_ntpm["TP53"] * 10
    # APM gene B2M is silenced on RNA axis ONLY (dna_status stays unknown!)
    bundle_5p = res_5p.to_mhc_allele_status_bundle()
    assert bundle_5p.apm is not None
    b2m_gene = {g.gene: g for g in bundle_5p.apm.genes}["B2M"]
    assert b2m_gene.dna_status == "unknown"
    assert b2m_gene.rna_status == "silenced"
    assert "B2M" in bundle_5p.apm.disrupted_genes()
    # MHC allele ASE: HLA-A*02:01 silenced & masked; HLA-B*07:02 retained!
    dec_a = resolve_effective_status(bundle_5p.by_allele()["HLA-A*02:01"])
    dec_b = resolve_effective_status(bundle_5p.by_allele()["HLA-B*07:02"])
    assert dec_a.effective_status == "silenced" and dec_a.masked is True
    assert dec_b.effective_status == "retained" and dec_b.masked is False

    # Now test 10x 3' low-UMI guard (only 3 resolvable UMIs -> unresolvable, NOT silenced!)
    low_umi_tsv = tmp_path / "schla_3p.tsv"
    low_umi_tsv.write_text("barcode\tallele\tumi_count\nIMM_0\tHLA-A*02:01\t3\n")
    res_3p = analyze_scrna_tumor_clusters(
        counts, clusters,
        somatic_barcodes=somatic,
        mhc_allele_umis=load_scrna_mhc_allele_umis(low_umi_tsv),
        policy={"chemistry": "10x_3p"},
    )
    dec_3p = resolve_effective_status(res_3p.to_mhc_allele_status_bundle().by_allele()["HLA-A*02:01"])
    assert dec_3p.effective_status == "unresolvable"
    assert dec_3p.masked is False
    assert "scrna_3p_exon2_3_coverage_gap" in dec_3p.caveats

    # Write tumor_cluster_ntpm.tsv and verify columns
    ntpm_out = tmp_path / "tumor_ntpm.tsv"
    res_5p.write_tumor_cluster_ntpm_tsv(ntpm_out)
    text = ntpm_out.read_text()
    assert "tumor_cluster_ntpm\timmune_control_ntpm" in text
    assert "TP53\t" in text


def test_scrna_cross_species_mouse_and_cli_mhc_status_subcommand(tmp_path):
    # Mouse wide TSV matrix (mouse capitalization: Ptprc, Cd3d, B2m, H2-Kb/H2-Db)
    counts_tsv = tmp_path / "mouse_counts.tsv"
    bcs = [f"M_T{i}" for i in range(5)] + [f"M_I{i}" for i in range(5)]
    counts_tsv.write_text(
        "gene\t" + "\t".join(bcs) + "\n"
        "Ptprc\t" + "\t".join(["0"] * 5 + ["400"] * 5) + "\n"
        "Cd3d\t" + "\t".join(["0"] * 5 + ["300"] * 5) + "\n"
        "B2m\t" + "\t".join(["200"] * 5 + ["200"] * 5) + "\n"
        "Kras\t" + "\t".join(["800"] * 5 + ["100"] * 5) + "\n"
    )
    cl_tsv = tmp_path / "mouse_clusters.tsv"
    cl_tsv.write_text(
        "barcode\tcluster\n"
        + "\n".join(f"{bc}\t{'tumor_like' if bc.startswith('M_T') else 'immune'}" for bc in bcs) + "\n"
    )
    umi_tsv = tmp_path / "mouse_mhc_umis.tsv"
    umi_tsv.write_text(
        "barcode\tallele\tumi_count\n"
        + "\n".join(f"M_I{i}\tH2-Kb\t10\nM_I{i}\tH2-Db\t10" for i in range(5)) + "\n"
        + "\n".join(f"M_T{i}\tH2-Db\t12" for i in range(5)) + "\n"
    )
    out_status = tmp_path / "mhc_status.json"
    out_ntpm = tmp_path / "mouse_tumor_ntpm.tsv"
    out_summary = tmp_path / "mouse_scrna_summary.json"

    run_cli([
        "mhc-status",
        "--scrna-counts", str(counts_tsv),
        "--scrna-clusters", str(cl_tsv),
        "--scrna-mhc-umis", str(umi_tsv),
        "--scrna-chemistry", "10x_5p",
        "--output", str(out_status),
        "--output-tumor-ntpm", str(out_ntpm),
        "--output-scrna-summary", str(out_summary),
    ])
    assert out_status.is_file()
    assert out_ntpm.is_file()
    assert out_summary.is_file()

    loaded = load_mhc_allele_status(out_status)
    assert resolve_effective_status(loaded.by_allele()["H2-K*b"]).effective_status == "silenced"
    assert resolve_effective_status(loaded.by_allele()["H2-D*b"]).effective_status == "retained"
    summary = json.loads(out_summary.read_text())
    assert summary["tumor_cluster_ids"] == ["tumor_like"]
    assert summary["immune_control_cluster_ids"] == ["immune"]

    # YAML schema & policy validation
    cfg = load_vaxrank_config(config_path="builtin:default")
    assert MHCExpressionPolicy.resolve(cfg["mhc_expression"]).kmer_size == 21
    assert SCRNATumorClusterPolicy.resolve(cfg["scrna_tumor"]).min_scrna_resolvable_umis == 20
    with pytest.raises(ValueError, match="Unknown mhc_expression keys"):
        MHCExpressionPolicy.resolve({"bad_knob": 1})
    with pytest.raises(ValueError, match="Unknown scrna_tumor keys"):
        SCRNATumorClusterPolicy.resolve({"bad_knob": 1})
