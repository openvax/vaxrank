"""Comprehensive tests for cross-species, multi-modality AntigenClonality (#4)."""

from types import SimpleNamespace

import pytest

from vaxrank.clonality import (
    CallerClonalityEstimate,
    ClonalityBundle,
    ClonalityPolicy,
    SampleClonalityObservation,
    apply_clonality_to_ranked,
    combine_genomic_callers,
    combine_sample_clonality,
    combine_spatial_and_longitudinal,
    detect_and_load_clonality_file,
    estimate_direct_dna_ccf,
    estimate_scrna_allele_expressional_clonality,
    estimate_scrna_gene_expressional_clonality,
    load_canonical_clonality_table,
    load_cellsnp_lite_dir,
    load_facets_maf,
    load_flames_transcript_counts,
    load_isoquant_read_info,
    load_purecn_variants,
    load_purple_purity_tsv,
    load_purple_vcf,
    load_pyclone_vi,
    load_sequenza_segments,
    load_vartrix_matrices,
    merge_clonality_records,
)


# ---------------------------------------------------------------------------
# 1. Direct Cross-Species DNA CCF + Mitochondrial Boundary
# ---------------------------------------------------------------------------

def test_direct_dna_ccf_clonal_subclonal_and_low_depth_uncertain():
    """Beta-Binomial CCF distinguishes clonal, subclonal, and shallow-depth uncertain variants."""
    # Purity 0.60, diploid CN=2, m=1 -> clonal expected VAF = 0.30 (30/100)
    clonal = estimate_direct_dna_ccf(alt_reads=30, total_depth=100, tumor_purity=0.60, contig="chr7")
    assert clonal.source_tool == "direct_diploid"
    assert clonal.assumed_diploid is True
    assert "assumed_diploid_cn" in clonal.notes
    assert clonal.ccf >= 0.85
    assert clonal.clonal_probability >= 0.65

    # Subclonal at VAF 0.06 (6/100) with purity 0.60 -> CCF ~ 0.20
    subclonal = estimate_direct_dna_ccf(alt_reads=6, total_depth=100, tumor_purity=0.60, contig="chr7")
    assert subclonal.ccf <= 0.30
    assert subclonal.ci_upper < 0.60
    assert subclonal.clonal_probability < 0.05

    # Shallow depth (2/8) -> wide CI flagged uncertain
    shallow = estimate_direct_dna_ccf(alt_reads=2, total_depth=8, tumor_purity=0.60, contig="chr7")
    assert "low_depth_uncertain" in shallow.notes
    assert (shallow.ci_upper - shallow.ci_lower) > 0.45


def test_direct_dna_ccf_copy_number_wgd_loh_and_cross_species_sex_ploidy():
    """Accounts for CN gain/LOH and cross-species hemizygous X chromosome ploidy (dog/mouse/human)."""
    # Copy-neutral LOH (CNt=2, major=2, minor=0) at purity 0.50: clonal m=2 expected VAF = 0.50
    cnloh = estimate_direct_dna_ccf(
        alt_reads=50, total_depth=100, tumor_purity=0.50,
        contig="chr17", cn_total=2.0, cn_major=2.0, cn_minor=0.0,
    )
    assert cnloh.source_tool == "direct_cn"
    assert "local_loh_or_cnloh" in cnloh.notes
    assert cnloh.multiplicity == 2.0
    assert cnloh.ccf >= 0.85

    # Male murine/canine hemizygous chrX (normal ploidy = 1, tumor CNt = 1, m = 1)
    # At purity 0.60, clonal expected VAF = 0.60 / (0.60*1 + 0.40*1) = 0.60
    male_x = estimate_direct_dna_ccf(
        alt_reads=60, total_depth=100, tumor_purity=0.60,
        contig="chrX", sex="male",
    )
    assert male_x.ccf >= 0.85
    assert male_x.clonal_probability >= 0.70


def test_mitochondrial_boundary_refuses_diploid_ccf():
    """#4 osteosarc MT-ND5 boundary: never compute nuclear diploid CCF on MT/chrM contigs."""
    mt = estimate_direct_dna_ccf(
        alt_reads=42, total_depth=100, tumor_purity=0.50,
        contig="MT", normal_alt_reads=1, normal_depth=200,
    )
    assert mt.compartment == "mitochondrial"
    assert mt.ccf is None
    assert mt.clonal_probability is None
    assert mt.mt_heteroplasmy_tumor == pytest.approx(0.42)
    assert mt.mt_heteroplasmy_normal == pytest.approx(0.005)
    assert "mitochondrial_non_diploid_ccf_not_applicable" in mt.notes

    genomic = combine_genomic_callers([mt])
    assert genomic.dna_ccf is None
    assert genomic.compartment == "mitochondrial"
    combined = combine_sample_clonality("MT:12345:A:G", "mutation", genomic=genomic)
    assert combined.effective_clonality is None
    assert combined.clonality_status == "unassessed"
    assert combined.combination_rule == "mt_heteroplasmy_non_diploid"


# ---------------------------------------------------------------------------
# 2. Pure-Python Upstream DNA Format Parsers + Multi-Caller Reconciliation
# ---------------------------------------------------------------------------

def test_upstream_dna_parsers_and_auto_detection(tmp_path):
    """Realistic synthetic fixtures for PyClone-VI, PURPLE, FACETS, PureCN, and Sequenza."""
    # PyClone-VI TSV
    pyclone = tmp_path / "pyclone_vi.tsv"
    pyclone.write_text(
        "mutation_id\tsample_id\tcluster_id\tcellular_prevalence\tcellular_prevalence_std\tcluster_assignment_prob\n"
        "7:140753336:A:T\tT1\t0\t0.94\t0.04\t0.99\n"
        "12:25245350:C:A\tT1\t1\t0.22\t0.05\t0.95\n",
        encoding="utf-8",
    )
    pc = load_pyclone_vi(pyclone)
    assert pc["7:140753336:A:T"][0].ccf == pytest.approx(0.94)
    assert pc["7:140753336:A:T"][0].cluster_id == "0"
    assert detect_and_load_clonality_file(pyclone)["12:25245350:C:A"].clonality_status == "subclonal"

    # PURPLE somatic VCF + purity TSV (including an MT variant)
    purity_tsv = tmp_path / "patient.purple.purity.tsv"
    purity_tsv.write_text(
        "purity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\twholeGenomeDuplication\tqcStatus\n"
        "0.68\t1.01\t0.42\t0.70\t2.15\tFEMALE\tNORMAL\t0.10\tfalse\tPASS\n",
        encoding="utf-8",
    )
    assert load_purple_purity_tsv(purity_tsv)["purity"] == pytest.approx(0.68)

    purple_vcf = tmp_path / "patient.purple.somatic.vcf"
    purple_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t140753336\t.\tA\tT\t.\tPASS\tPURPLE_CN=2.01;PURPLE_MACN=1.02;PURPLE_VCN=1.03;SUBCL=0.02;PURPLE_AF=0.35\n"
        "chr12\t25245350\t.\tC\tA\t.\tPASS\tPURPLE_CN=2.00;PURPLE_MACN=1.00;PURPLE_VCN=0.24;SUBCL=0.91;PURPLE_AF=0.08\n"
        "chrM\t13105\t.\tA\tG\t.\tPASS\tPURPLE_AF=0.61\n",
        encoding="utf-8",
    )
    pv = load_purple_vcf(purple_vcf, purity_tsv_path=purity_tsv)
    assert pv["7:140753336:A:T"][0].ccf == pytest.approx(1.0)
    assert pv["7:140753336:A:T"][0].clonal_probability == pytest.approx(0.98)
    assert pv["12:25245350:C:A"][0].ccf == pytest.approx(0.24)
    assert pv["MT:13105:A:G" if "MT:13105:A:G" in pv else "M:13105:A:G"][0].compartment == "mitochondrial"
    assert pv["M:13105:A:G"][0].ccf is None

    # FACETS / facets-suite annotated MAF
    facets_maf = tmp_path / "patient.facets.maf"
    facets_maf.write_text(
        "#version 2.4\n"
        "Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\t"
        "t_alt_count\tt_ref_count\ttcn.em\tlcn.em\tccf_expected_copies\tccf_expected_copies_lower\t"
        "ccf_expected_copies_upper\tclonality\tpurity\n"
        "7\t140753336\tA\tT\tT1\t34\t66\t2\t1\t0.92\t0.81\t1.00\tCLONAL\t0.68\n",
        encoding="utf-8",
    )
    fm = load_facets_maf(facets_maf)
    assert fm["7:140753336:A:T"][0].ccf == pytest.approx(0.92)
    assert fm["7:140753336:A:T"][0].local_cn_total == 2.0
    assert detect_and_load_clonality_file(facets_maf)["7:140753336:A:T"].clonality_status == "clonal"

    # PureCN *_variants.csv
    purecn_csv = tmp_path / "patient_variants.csv"
    purecn_csv.write_text(
        "Sampleid,chr,start,end,REF,ALT,AR,DEPTH,ML.C,ML.M.SEGMENT,ML.M,CELLFRACTION,CELLFRACTION.95.LOWER,CELLFRACTION.95.UPPER,POSTERIOR.SOMATIC\n"
        "T1,chr7,140753336,140753336,A,T,0.34,100,2,1,1,0.95,0.82,1.0,0.99\n",
        encoding="utf-8",
    )
    pn = load_purecn_variants(purecn_csv)
    assert pn["7:140753336:A:T"][0].ccf == pytest.approx(0.95)
    assert detect_and_load_clonality_file(purecn_csv)["7:140753336:A:T"].clonality_status == "clonal"

    # Sequenza segments + confints feeding into direct CCF
    seq_segs = tmp_path / "canine_segments.txt"
    seq_segs.write_text(
        "chromosome\tstart.pos\tend.pos\tCNt\tA\tB\n"
        "chr17\t100000\t5000000\t3\t2\t1\n",
        encoding="utf-8",
    )
    seq_cp = tmp_path / "canine_confints_CP.txt"
    seq_cp.write_text(
        "cellularity\tploidy.estimate\n"
        "0.55\t2.1\n"
        "0.62\t2.3\n"
        "0.70\t2.5\n",
        encoding="utf-8",
    )
    segs, purity, ploidy = load_sequenza_segments(seq_segs, confints_path=seq_cp)
    assert purity == pytest.approx(0.62)
    assert ploidy == pytest.approx(2.3)
    est = estimate_direct_dna_ccf(
        alt_reads=25, total_depth=100, tumor_purity=purity,
        contig="chr17", position=500000, cn_segments=segs,
    )
    assert est.local_cn_total == 3.0
    assert est.local_cn_minor == 1.0
    assert est.ccf >= 0.85


def test_multi_caller_reconciliation_precedence_and_discordance_envelope():
    """Multiple DNA callers on the same reads follow precedence and widen CI when discordant."""
    # Concordant callers -> PyClone-VI wins by precedence, no discordance
    est_pyclone = CallerClonalityEstimate(source_tool="pyclone_vi", ccf=0.93, ci_lower=0.85, ci_upper=1.0)
    est_purple = CallerClonalityEstimate(source_tool="purple", ccf=0.90, ci_lower=0.80, ci_upper=1.0)
    concordant = combine_genomic_callers([est_purple, est_pyclone])
    assert concordant.primary_source == "pyclone_vi"
    assert concordant.caller_discordance is False
    assert concordant.caller_spread == pytest.approx(0.03)

    # Discordant callers (0.92 vs 0.35 -> spread 0.57 > 0.25) -> widens CI envelope & downgrades to uncertain!
    est_facets_low = CallerClonalityEstimate(source_tool="facets", ccf=0.35, ci_lower=0.22, ci_upper=0.48)
    discordant = combine_genomic_callers([est_facets_low, est_pyclone])
    assert discordant.caller_discordance is True
    assert discordant.ci_lower == pytest.approx(0.22)
    assert discordant.ci_upper == pytest.approx(1.0)
    rec = combine_sample_clonality("7:140753336:A:T", "mutation", genomic=discordant)
    assert rec.clonality_status == "uncertain"
    assert "caller_discordance" in rec.divergences

    # merge_clonality_records reconciles multiple loaded maps instead of overwriting
    m1 = {"7:140753336:A:T": combine_sample_clonality("7:140753336:A:T", "mutation", genomic=combine_genomic_callers([est_pyclone]))}
    m2 = {"7:140753336:A:T": combine_sample_clonality("7:140753336:A:T", "mutation", genomic=combine_genomic_callers([est_facets_low]))}
    merged = merge_clonality_records([m1, m2])
    assert merged["7:140753336:A:T"].genomic.caller_discordance is True
    assert merged["7:140753336:A:T"].clonality_status == "uncertain"


# ---------------------------------------------------------------------------
# 3. Short-Read & Long-Read scRNA Parsers + Dropout/Abundance Confounding
# ---------------------------------------------------------------------------

def test_scrna_gene_bimodal_subclonal_vs_low_abundance_subcluster_fallback():
    """High-UMI bimodal CTA is called subclonal; low-UMI clonal CTA falls back to tumor subcluster pseudobulk."""
    # Case 1: Bimodal subclonal CTA — 12 UMIs/cell in 20 tumor cells, 0 in 80 tumor cells
    bimodal_counts = {f"BC{i}": (12 if i < 20 else 0) for i in range(100)}
    subclusters = {f"BC{i}": f"T_sub_{i // 25}" for i in range(100)}
    bimodal = estimate_scrna_gene_expressional_clonality(
        bimodal_counts, cell_to_subcluster=subclusters, assay_platform="short_read_3prime",
    )
    assert bimodal.is_bimodal_heterogeneous is True
    assert bimodal.abundance_confounded is False
    assert bimodal.detection_power > 0.99
    assert bimodal.selected_scale == "per_cell"
    assert bimodal.selected_clonality == pytest.approx(0.20, abs=0.01)
    rec_bi = combine_sample_clonality("MAGEA3", "CTA", expressional=bimodal)
    assert rec_bi.clonality_status == "subclonal"

    # Case 2: Low-abundance clonal CTA — 1 UMI in 15% of tumor cells evenly across all 4 subclusters
    # Per-cell mean among positives is ~1.0 -> very low detection power -> must NOT falsely call subclonal!
    low_abund_counts = {f"BC{i}": (1 if i % 7 == 0 else 0) for i in range(100)}
    low_abund = estimate_scrna_gene_expressional_clonality(
        low_abund_counts, cell_to_subcluster=subclusters, assay_platform="short_read_3prime",
    )
    assert low_abund.abundance_confounded is True
    assert low_abund.is_bimodal_heterogeneous is False
    assert low_abund.selected_scale == "subcluster_pseudobulk"
    assert low_abund.expressing_subcluster_fraction == pytest.approx(1.0)
    assert low_abund.selected_clonality == pytest.approx(1.0)
    rec_low = combine_sample_clonality("CTAG1B", "CTA", expressional=low_abund)
    assert rec_low.clonality_status != "subclonal"


def test_short_read_and_long_read_scrna_parsers_and_bulk_separation(tmp_path):
    """Verify VarTrix, cellsnp-lite, IsoQuant (sc vs bulk), and FLAMES parsers end-to-end."""
    # VarTrix 1-based MatrixMarket
    (tmp_path / "barcodes.tsv").write_text("BC1\nBC2\nBC3\nBC4\nBC5\nBC6\n", encoding="utf-8")
    (tmp_path / "vars.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t140753336\t.\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    (tmp_path / "alt.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n1 6 4\n1 1 3\n1 2 2\n1 3 4\n1 4 2\n", encoding="utf-8")
    (tmp_path / "ref.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n1 6 5\n1 1 2\n1 2 3\n1 3 2\n1 4 3\n1 5 5\n", encoding="utf-8")
    vt = load_vartrix_matrices(
        tmp_path / "alt.mtx", tmp_path / "ref.mtx",
        tmp_path / "barcodes.tsv", tmp_path / "vars.vcf",
        tumor_barcodes={"BC1", "BC2", "BC3", "BC4", "BC5"},
    )
    expr_vt = estimate_scrna_allele_expressional_clonality(
        vt["7:140753336:A:T"], n_total_tumor_cells=10, assay_platform="short_read_5prime",
    )
    assert expr_vt.n_covered_tumor_cells == 5
    assert expr_vt.n_expressing_tumor_cells == 4
    assert expr_vt.locus_coverage_fraction == pytest.approx(0.5)
    assert expr_vt.detection_power > 0.90
    assert expr_vt.selected_clonality >= 0.78

    # cellsnp-lite directory
    cs_dir = tmp_path / "cellsnp"
    cs_dir.mkdir()
    (cs_dir / "cellSNP.samples.tsv").write_text("BC1\nBC2\nBC3\nBC4\nBC5\n", encoding="utf-8")
    (cs_dir / "cellSNP.base.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr12\t25245350\t.\tC\tA\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    (cs_dir / "cellSNP.tag.AD.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n1 5 5\n1 1 4\n1 2 3\n1 3 5\n1 4 2\n1 5 4\n", encoding="utf-8")
    (cs_dir / "cellSNP.tag.DP.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n1 5 5\n1 1 8\n1 2 6\n1 3 9\n1 4 5\n1 5 8\n", encoding="utf-8")
    cs = load_cellsnp_lite_dir(cs_dir)
    assert cs["12:25245350:C:A"]["BC1"] == (4, 4)

    # IsoQuant single-cell long-read read_info.tsv (UMI deduplicated!)
    iq_sc = tmp_path / "isoquant_sc.read_info.tsv"
    rows = ["read_id\tchr\tstrand\tgene_id\tisoform_id\tbarcode\tumi"]
    for i in range(10):
        # All 10 tumor cells express cryptic isoform ISO_CRYPTIC with 3 unique UMIs + 2 canonical UMIs
        for u in range(3):
            rows.append(f"r_{i}_{u}_a\tchr1\t+\tGENE1\tISO_CRYPTIC\tBC{i}\tUMI_C_{u}")
            rows.append(f"r_{i}_{u}_dup\tchr1\t+\tGENE1\tISO_CRYPTIC\tBC{i}\tUMI_C_{u}")  # PCR duplicate
        for u in range(2):
            rows.append(f"r_{i}_{u}_b\tchr1\t+\tGENE1\tISO_CANON\tBC{i}\tUMI_W_{u}")
    iq_sc.write_text("\n".join(rows) + "\n", encoding="utf-8")
    iq_counts, is_sc = load_isoquant_read_info(iq_sc)
    assert is_sc is True
    assert iq_counts["ISO_CRYPTIC"]["BC0"] == (3, 2)  # deduplicated UMIs!
    expr_lr = estimate_scrna_allele_expressional_clonality(
        iq_counts["ISO_CRYPTIC"], n_total_tumor_cells=10,
        assay_platform="long_read_ont", measurement_level="isoform",
    )
    assert expr_lr.locus_coverage_fraction == pytest.approx(1.0)
    assert expr_lr.selected_clonality >= 0.95

    # IsoQuant BULK long-read (barcode == '.') must NOT produce single-cell clonality
    iq_bulk = tmp_path / "isoquant_bulk.read_info.tsv"
    iq_bulk.write_text(
        "read_id\tchr\tstrand\tgene_id\tisoform_id\tbarcode\tumi\n"
        "r1\tchr1\t+\tGENE1\tISO_CRYPTIC\t.\t.\n"
        "r2\tchr1\t+\tGENE1\tISO_CANON\t.\t.\n",
        encoding="utf-8",
    )
    bulk_counts, is_sc_bulk = load_isoquant_read_info(iq_bulk)
    assert is_sc_bulk is False
    expr_bulk = estimate_scrna_allele_expressional_clonality(
        bulk_counts["ISO_CRYPTIC"], assay_platform="bulk_long_read",
    )
    assert expr_bulk.selected_clonality is None
    assert expr_bulk.bulk_expression_value == pytest.approx(0.5)

    # FLAMES transcript_count.csv
    flames = tmp_path / "transcript_count.csv"
    flames.write_text(
        "transcript_id,gene_id,BC0,BC1,BC2,BC3,BC4,BC5\n"
        "ENST_SPLICE_1,GENE1,5,4,6,5,3,4\n",
        encoding="utf-8",
    )
    fl = load_flames_transcript_counts(flames)
    expr_fl = estimate_scrna_gene_expressional_clonality(
        fl["ENST_SPLICE_1"], assay_platform="long_read_pacbio", measurement_level="isoform",
    )
    assert expr_fl.selected_clonality >= 0.95


# ---------------------------------------------------------------------------
# 4. Cross-Axis Combination Without Double-Counting + Silencing Divergence
# ---------------------------------------------------------------------------

def test_cross_axis_combination_no_double_counting_and_silencing_detection():
    """Mutation-specific scRNA measures joint P(G,E) — never multiplied by C_DNA!"""
    # 50% subclone with DNA CCF = 0.50 AND long-read scRNA mutant fraction = 0.50
    # True presenting fraction is ~0.50 (NOT 0.25!)
    dna_50 = combine_genomic_callers([
        CallerClonalityEstimate(source_tool="pyclone_vi", ccf=0.50, ci_lower=0.40, ci_upper=0.60),
    ])
    lr_cells_50 = {f"BC{i}": ((4, 4) if i < 20 else (0, 8)) for i in range(40)}
    expr_50 = estimate_scrna_allele_expressional_clonality(
        lr_cells_50, n_total_tumor_cells=40, assay_platform="long_read_ont",
        measurement_level="allele_or_junction",
    )
    fused = combine_sample_clonality("chr1:100:A:G", "mutation", genomic=dna_50, expressional=expr_50)
    assert fused.combination_rule == "precision_weighted_dna_and_allele_scrna"
    assert fused.effective_clonality == pytest.approx(0.50, abs=0.05)
    assert (fused.ci_upper - fused.ci_lower) < (dna_50.ci_upper - dna_50.ci_lower)

    # Clonal DNA (0.96) + well-powered scRNA silencing (mutant only in 15% of covered cells)
    dna_clonal = combine_genomic_callers([
        CallerClonalityEstimate(source_tool="pyclone_vi", ccf=0.96, ci_lower=0.88, ci_upper=1.0),
    ])
    lr_silenced = {f"BC{i}": ((5, 4) if i < 6 else (0, 9)) for i in range(40)}
    expr_silenced = estimate_scrna_allele_expressional_clonality(
        lr_silenced, n_total_tumor_cells=40, assay_platform="long_read_pacbio",
        measurement_level="allele_or_junction",
    )
    silenced = combine_sample_clonality("chr2:200:C:T", "mutation", genomic=dna_clonal, expressional=expr_silenced)
    assert silenced.combination_rule == "expression_silencing_bottleneck"
    assert "transcriptional_silencing_or_nmd_suspected" in silenced.divergences
    assert silenced.effective_clonality <= 0.25
    assert silenced.clonality_status == "subclonal"

    # Gene-level scRNA bottleneck (min vs product policy)
    gene_expr = estimate_scrna_gene_expressional_clonality(
        {f"BC{i}": (10 if i < 12 else 0) for i in range(40)},
        assay_platform="short_read_3prime", measurement_level="gene",
    )
    rec_min = combine_sample_clonality("chr3:300:G:A", "mutation", genomic=dna_clonal, expressional=gene_expr)
    assert rec_min.combination_rule == "min_dna_and_gene_level_scrna"
    assert "gene_level_expression_bottleneck" in rec_min.divergences
    assert rec_min.effective_clonality == pytest.approx(0.30, abs=0.02)

    rec_prod = combine_sample_clonality(
        "chr3:300:G:A", "mutation", genomic=dna_clonal, expressional=gene_expr,
        policy=ClonalityPolicy(gene_level_combination="product"),
    )
    assert rec_prod.combination_rule == "dna_times_gene_level_scrna"
    assert rec_prod.effective_clonality == pytest.approx(0.96 * 0.30, abs=0.02)


# ---------------------------------------------------------------------------
# 5. Spatial Multi-Region + Directed Longitudinal T0 -> T1 -> T2 -> T3 (Sid-style)
# ---------------------------------------------------------------------------

def test_longitudinal_directed_trajectories_and_multi_sequencing_per_timepoint():
    """Sid-style multi-sequencing-per-timepoint + directed T0->T1->T2->T3 trajectories."""
    def _obs(sid, tp_i, tp_lbl, region, ccf, tool="purple", depth=80, present=None):
        genomic = (
            combine_genomic_callers([CallerClonalityEstimate(source_tool=tool, ccf=ccf, depth=depth, sample_id=sid)])
            if ccf is not None else None
        )
        return SampleClonalityObservation(
            sample_id=sid, timepoint_index=tp_i, timepoint_label=tp_lbl,
            region_label=region, genomic=genomic, depth=depth, present=present,
        )

    # 1) Persistent across T0, T1 (with TWO sequencing runs on T1!), T2, T3
    persistent = combine_spatial_and_longitudinal("DYNC1H1_p.R100Q", "mutation", [
        _obs("sid_wes_t0", 0, "T0", "resection", 0.94, tool="facets"),
        _obs("sid_wes_t1_a", 1, "T1", "met1", 0.92, tool="facets"),
        _obs("sid_wes_t1_b", 1, "T1", "met1", 0.95, tool="pyclone_vi"),  # second sequencing of same T1 sample!
        _obs("sid_wes_t3", 3, "T3", "met2", 0.96, tool="purple"),
    ])
    assert persistent.longitudinal.trajectory == "persistent"
    assert persistent.effective_clonality >= 0.90
    assert persistent.clonality_status == "clonal"
    # T1 had 2 sequencing runs of the SAME region "met1" -> must count as 1 region, not 2!
    t1_spatial = [s for s in persistent.longitudinal.per_timepoint_spatial if s.timepoint_label == "T1"][0]
    assert t1_spatial.regions_sampled == 1

    # 2) Emergent clonal sweep in T3 only (absent at T0 & T1, 0.95 clonal at T3) -> keeps 0.95 clonal!
    emergent = combine_spatial_and_longitudinal("EMERGENT_T3", "mutation", [
        _obs("t0", 0, "T0", "resection", None, depth=90, present=False),
        _obs("t1", 1, "T1", "met1", None, depth=80, present=False),
        _obs("t3", 3, "T3", "met2", 0.95, depth=100, present=True),
    ])
    assert emergent.longitudinal.trajectory == "emergent_latest"
    assert emergent.effective_clonality == pytest.approx(0.95)
    assert emergent.clonality_status == "clonal"

    # 3) Lost in T3 (clonal in T0 & T1, absent at 80x depth in T3) -> effective=0.0, status="lost"!
    lost = combine_spatial_and_longitudinal("LOST_IN_T3", "mutation", [
        _obs("t0", 0, "T0", "resection", 0.93, depth=90, present=True),
        _obs("t1", 1, "T1", "met1", 0.85, depth=80, present=True),
        _obs("t3", 3, "T3", "met2", None, depth=80, present=False),
    ])
    assert lost.longitudinal.trajectory == "lost_in_latest"
    assert lost.effective_clonality == 0.0
    assert lost.clonality_status == "lost"
    assert "lost_in_latest_timepoint" in lost.divergences

    # 4) Low depth at T3 (depth=4, present=None) -> latest_unassessed (never falsely called lost!)
    unassessed_latest = combine_spatial_and_longitudinal("LOW_COV_T3", "mutation", [
        _obs("t0", 0, "T0", "resection", 0.92, depth=90, present=True),
        _obs("t3", 3, "T3", "met2", None, depth=4, present=None),
    ])
    assert unassessed_latest.longitudinal.trajectory == "latest_unassessed"
    assert unassessed_latest.clonality_status == "unassessed"

    # 5) Expanding & Contracting trajectories
    expanding = combine_spatial_and_longitudinal("EXPANDING", "mutation", [
        _obs("t0", 0, "T0", "primary", 0.20),
        _obs("t1", 1, "T1", "met1", 0.55),
        _obs("t3", 3, "T3", "met2", 0.92),
    ])
    assert expanding.longitudinal.trajectory == "expanding"

    contracting = combine_spatial_and_longitudinal("CONTRACTING", "mutation", [
        _obs("t0", 0, "T0", "primary", 0.92),
        _obs("t1", 1, "T1", "met1", 0.55),
        _obs("t3", 3, "T3", "met2", 0.18),
    ])
    assert contracting.longitudinal.trajectory == "contracting"

    # 6) Multi-region primary resection at T0 (3 physical regions, only 1/3 carries variant)
    regional = combine_spatial_and_longitudinal("REGIONAL_SNV", "mutation", [
        _obs("regA", 0, "T0", "region_A", 0.90, depth=80),
        _obs("regB", 0, "T0", "region_B", None, depth=80, present=False),
        _obs("regC", 0, "T0", "region_C", None, depth=80, present=False),
    ])
    assert regional.spatial.regions_sampled == 3
    assert regional.spatial.regions_present == 1
    assert regional.spatial.is_spatially_ubiquitous is False
    assert regional.effective_clonality == pytest.approx(0.30)


# ---------------------------------------------------------------------------
# 6. Mixed-Class Asymmetric Policy + scRNA Tumor Bridge + Config & Dossier
# ---------------------------------------------------------------------------

def test_mixed_antigen_asymmetric_policy_and_dossier_integration(tmp_path):
    """Mixed SNV / SV / CTA / viral run: penalizes verified subclonality, never missing assays."""
    policy_gate = ClonalityPolicy(
        mode="gate",
        min_clonality=0.30,
        unassessed_weight=0.85,
        by_antigen_kind={
            "viral": {"integrated_driver_prior": 0.95},
        },
    )
    snv_clonal = combine_sample_clonality(
        "SNV_CLONAL", "mutation",
        genomic=combine_genomic_callers([CallerClonalityEstimate(source_tool="pyclone_vi", ccf=0.95, ci_lower=0.88, ci_upper=1.0)]),
        policy=policy_gate,
    )
    snv_subclonal = combine_sample_clonality(
        "SNV_SUBCLONAL", "mutation",
        genomic=combine_genomic_callers([CallerClonalityEstimate(source_tool="pyclone_vi", ccf=0.15, ci_lower=0.08, ci_upper=0.25)]),
        policy=policy_gate,
    )
    sv_unassessed = combine_sample_clonality("SV_FUSION", "SV", policy=policy_gate)
    cta_bulk = combine_sample_clonality("MAGEA3", "CTA", policy=policy_gate)
    hpv_e7 = combine_sample_clonality("HPV16_E7", "viral", policy=policy_gate)

    assert sv_unassessed.clonality_status == "unassessed"
    assert cta_bulk.clonality_status == "unassessed"
    assert hpv_e7.combination_rule == "viral_integrated_driver_prior_assumed"
    assert "prior_assumed_not_measured" in hpv_e7.caveats

    # Gate mode drops ONLY proven subclonal SNV; keeps unassessed SV & bulk CTA!
    bundle_gate = ClonalityBundle(
        policy=policy_gate,
        records={r.antigen_id: r for r in (snv_clonal, snv_subclonal, sv_unassessed, cta_bulk, hpv_e7)},
    )
    ranked = [
        (SimpleNamespace(short_description="SNV_SUBCLONAL", antigen_kind="mutation"), [SimpleNamespace(combined_score=10.0)]),
        (SimpleNamespace(short_description="SV_FUSION", antigen_kind="SV"), [SimpleNamespace(combined_score=8.0)]),
        (SimpleNamespace(short_description="MAGEA3", antigen_kind="CTA"), [SimpleNamespace(combined_score=7.5)]),
        (SimpleNamespace(short_description="SNV_CLONAL", antigen_kind="mutation"), [SimpleNamespace(combined_score=7.0)]),
        (SimpleNamespace(short_description="HPV16_E7", antigen_kind="viral"), [SimpleNamespace(combined_score=6.5)]),
    ]
    gated, audit = apply_clonality_to_ranked(ranked, bundle_gate)
    kept_ids = [s.short_description for s, _ in gated]
    assert "SNV_SUBCLONAL" not in kept_ids
    assert set(kept_ids) == {"SV_FUSION", "MAGEA3", "SNV_CLONAL", "HPV16_E7"}
    assert len(audit["dropped"]) == 1

    # Weight mode with unassessed_weight=0.85: equal base score -> clonal > unassessed > subclonal
    policy_weight = ClonalityPolicy(mode="weight", unassessed_weight=0.85)
    assert snv_clonal.score_multiplier(policy_weight) > sv_unassessed.score_multiplier(policy_weight) > snv_subclonal.score_multiplier(policy_weight)

    # Canonical JSON/TSV multi-timepoint table roundtrip
    canon_tsv = tmp_path / "canonical_clonality.tsv"
    canon_tsv.write_text(
        "antigen_id\tantigen_kind\tsample_id\ttimepoint_index\ttimepoint_label\tregion_label\tdna_ccf\tdepth\tpresent\n"
        "SNV1\tmutation\ts0\t0\tT0\tprimary\t0.92\t80\ttrue\n"
        "SNV1\tmutation\ts1\t1\tT1\trelapse\t0.95\t90\ttrue\n",
        encoding="utf-8",
    )
    loaded_canon = load_canonical_clonality_table(canon_tsv)
    assert loaded_canon["SNV1"].longitudinal.trajectory == "persistent"
    assert loaded_canon["SNV1"].clonality_status == "clonal"


def test_scrna_tumor_cluster_bridge_and_config_schema():
    """Bridge from scrna_tumor.py tumor clusters to ExpressionalClonality + YAML config schema."""
    from vaxrank.config import load_vaxrank_config
    from vaxrank.scrna_tumor import (
        SCRNACountMatrix,
        analyze_scrna_tumor_clusters,
        estimate_scrna_tumor_expressional_clonality,
    )

    # 20 tumor cells across 2 tumor clusters + 10 CD45+ immune cells
    clusters = {
        **{f"T1_{i}": "tumor_A" for i in range(10)},
        **{f"T2_{i}": "tumor_B" for i in range(10)},
        **{f"I_{i}": "immune" for i in range(10)},
    }
    counts_by_barcode = {
        **{f"T1_{i}": {"SOX9": 100.0, "HPV16_E7": 10.0, "MAGEA3": 15.0} for i in range(10)},
        **{f"T2_{i}": {"SOX9": 100.0, "HPV16_E7": 12.0} for i in range(10)},
        **{f"I_{i}": {"PTPRC": 80.0, "CD3D": 60.0} for i in range(10)},
    }
    matrix = SCRNACountMatrix(
        genes=("PTPRC", "CD3D", "SOX9", "HPV16_E7", "MAGEA3"),
        barcodes=tuple(clusters.keys()),
        counts_by_barcode=counts_by_barcode,
    )
    tumor_res = analyze_scrna_tumor_clusters(
        matrix, clusters, policy={"tumor_markers": ["SOX9"]},
    )
    clon = estimate_scrna_tumor_expressional_clonality(
        matrix, clusters, tumor_res, ["HPV16_E7", "MAGEA3"],
        assay_platform="long_read_ont",
    )
    assert clon["HPV16_E7"].selected_clonality == pytest.approx(1.0)
    assert clon["MAGEA3"].selected_clonality == pytest.approx(0.5, abs=0.05)

    # YAML schema accepts clonality: overrides cleanly
    cfg = load_vaxrank_config(set_overrides=["clonality.mode=weight", "clonality.unassessed_weight=0.85"])
    assert cfg["clonality"]["mode"] == "weight"
    assert cfg["clonality"]["unassessed_weight"] == pytest.approx(0.85)


def test_tumor_board_dossier_embeds_clonality_and_flags_subclonal_or_lost():
    """Tumor board dossier includes Layer 3 clonality and emits soft review flags."""
    from vaxrank.tumor_board import build_tumor_board_dossier

    fake_vp = SimpleNamespace(
        amino_acids="SIINFEKLAAAAAAAAAAAA",
        combined_score=4.2,
        target_epitope_score=3.5,
        self_epitope_score=0.0,
        expression_score=1.2,
        target_epitopes=(),
        self_epitopes=(),
        antigen=None,
        mutant_protein_fragment=SimpleNamespace(
            gene_name="KRAS",
            n_alt_reads=12,
            n_ref_reads=18,
            n_other_reads=0,
            dna_vaf=0.10,
            rna_vaf=0.40,
            is_dna_only_proxy=False,
        ),
        manufacturability_scores=None,
    )
    src = SimpleNamespace(short_description="chr12:25245350:C:A", gene_name="KRAS")
    subclonal = combine_sample_clonality(
        "chr12:25245350:C:A", "mutation",
        genomic=combine_genomic_callers([
            CallerClonalityEstimate(source_tool="pyclone_vi", ccf=0.18, ci_lower=0.10, ci_upper=0.28),
        ]),
    )
    bundle = ClonalityBundle(records={"chr12:25245350:C:A": subclonal})
    dossier = build_tumor_board_dossier([(src, [fake_vp])], clonality_bundle=bundle)
    cand = dossier["candidates"][0]
    assert cand["tumor_evolution_and_mhc"]["clonality"]["clonality_status"] == "subclonal"
    flag_names = {f["flag"] for f in cand["review_flags"]}
    assert "subclonal_antigen" in flag_names
    assert cand["default_recommendation"] == "Review"
