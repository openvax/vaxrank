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

"""End-to-end tests covering the LENS / pVACseq → vaccine-design path.

Pre-PR-#253: external-input mode hard-short-circuited and the
construct writers were never reached. These tests pin the new
architecture: a single shared ``ranked_variants_with_vaccine_peptides``
intermediate that *both* the VCF/BAM pipeline and the LENS / pVACseq
loaders produce, then a shared dispatch (``_emit_outputs``) into
modality-specific construct writers.
"""

import json
import os

from vaxrank.epitope_io import load_lens
from vaxrank.external_input import (
    _mut_offsets_in_context,
    _parse_variant_coords,
    ranked_from_lens_predictions,
)

DATA_DIR = os.path.join(
    os.path.dirname(__file__), "data", "epitope_fixtures")


def test_parse_variant_coords_typical_lens_form():
    v = _parse_variant_coords("chr17:7675088:C:T")
    assert v is not None
    assert v.contig == "chr17"
    assert v.start == 7675088
    assert v.ref == "C"
    assert v.alt == "T"


def test_parse_variant_coords_malformed_returns_none():
    assert _parse_variant_coords("garbage") is None
    assert _parse_variant_coords("chr1:notapos:A:T") is None
    assert _parse_variant_coords(None) is None


def test_parse_variant_coords_2_part_chr_pos_form():
    """Real LENS v1.9 reports emit ``chr1:26780312`` (2-part: chr+pos
    only, no ref/alt). Pre-fix this returned None and 3000+ warnings
    fired on a real LENS run. Now: parses with placeholder ref/alt that
    satisfy varcode's nucleotide validator."""
    v = _parse_variant_coords("chr1:26780312")
    assert v is not None
    assert v.contig == "chr1"
    assert v.start == 26780312
    # ref/alt are placeholder nucleotides; the (chr, pos) tuple is
    # what downstream construct assembly keys off, not the alleles.
    assert v.ref in ("A", "C", "G", "T")
    assert v.alt in ("A", "C", "G", "T")
    assert v.ref != v.alt


def test_parse_variant_coords_3_part_with_ref_only():
    v = _parse_variant_coords("chr1:26780312:C")
    assert v is not None
    assert v.ref == "C"
    # alt placeholder, not equal to ref
    assert v.alt != "C"


def test_parse_variant_coords_nan_string():
    """Pandas often reads empty TSV cells as the literal string 'nan'
    (not NaN). Treat both as 'genuinely empty' → None, no warning
    (caller skips silently for non-SNV antigen rows)."""
    assert _parse_variant_coords("nan") is None
    assert _parse_variant_coords("NaN") is None
    assert _parse_variant_coords("") is None
    assert _parse_variant_coords("   ") is None


def test_parse_variant_coords_n_nucleotide_rejected():
    """varcode rejects 'N' as a nucleotide; the parser substitutes
    safe placeholders instead of returning None."""
    v = _parse_variant_coords("chr1:1234:N:N")
    assert v is not None
    assert v.ref != "N"
    assert v.alt != "N"


def test_real_lens_v19_subset_produces_ranked_entries():
    """Regression test for the real-LENS-v1.9 path. The subset fixture
    came from a Hugo IPRES Pt10 dump and tripped both #259 (variant_coords
    parser) and #260 (duplicate-row handling). Pin that the path now
    produces a non-empty ranked list end-to-end."""
    path = os.path.join(
        DATA_DIR, "real_lens_subsets", "lens_v1.9_real_subset.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)
    # The fixture has at least one parseable variant_coords row.
    # Pre-fix: 0. Now: > 0.
    assert len(ranked) > 0, (
        "Expected real LENS v1.9 fixture to produce a non-empty ranked "
        "list after the variant_coords parser fix")


def test_mut_offsets_in_context_finds_peptide():
    # AASVVGSSSSSGTR contains SVVGSSSSS at offset 2, length 9
    start, end = _mut_offsets_in_context("SVVGSSSSS", "AASVVGSSSSSGTR")
    assert start == 2
    assert end == 11


def test_mut_offsets_in_context_falls_back_to_full_window():
    # When the peptide isn't found in the context, the whole window is
    # treated as the mutation span (so window-centering doesn't crop
    # it out).
    start, end = _mut_offsets_in_context("XYZXYZ", "AASVVGSSSSSGTR")
    assert start == 0
    assert end == len("AASVVGSSSSSGTR")


def test_lens_picks_strongest_binder_when_multiple_rows_per_variant():
    """Multi-row LENS fixture: one variant has three peptides at very
    different IC50s. The representative pick must be the strongest
    binder (lowest IC50), not the file-order first row.

    Pre-fix: ``_pick_representative`` looked for non-existent columns
    (``ic50``, ``percentile_rank``) and tied every row at (2, 0.0),
    so it returned the alphabetically-first allele's row regardless
    of binding strength. This test pins the fix.
    """
    path = os.path.join(DATA_DIR, "lens_multi_row_per_variant.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)
    assert len(ranked) == 1, (
        "Fixture has one variant with 3 peptides; got %d entries" % len(ranked))
    _, peptides = ranked[0]
    fragment = peptides[0].mutant_protein_fragment
    # Strongest binder is STRNGLVLLL (mhcflurry IC50 = 18.0). Its
    # pep_context is "AASTRNGLVLLLGTR".
    assert "STRNGLVLLL" in fragment.amino_acids, (
        "Expected pep_context to come from the strongest binder row "
        "(STRNGLVLLL, IC50=18); got %r" % fragment.amino_acids)


def test_lens_to_ranked_vaccine_peptides_round_trip():
    """Loading the LENS fixture and converting to ranked vaccine peptides
    should produce one entry per unique variant, each carrying the
    pep_context as the antigen fragment."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)
    assert len(ranked) == 3, (
        "LENS fixture has 3 unique variants; got %d" % len(ranked))
    # Each ranked entry: (Variant, [VaccinePeptide])
    for variant, peptides in ranked:
        assert len(peptides) == 1
        vp = peptides[0]
        # Fragment AA = pep_context (the SLP-style window from LENS)
        assert vp.mutant_protein_fragment.amino_acids
        # All three fixture rows have pep_context longer than peptide
        assert (len(vp.mutant_protein_fragment.amino_acids)
                > len(vp.epitope_predictions[0].peptide_sequence))
        # Mutation span is non-empty inside the fragment
        assert (vp.mutant_protein_fragment.mutant_amino_acid_end_offset
                > vp.mutant_protein_fragment.mutant_amino_acid_start_offset)


def test_lens_drives_mrna_construct_assembly_end_to_end():
    """LENS file → ranked vaccine peptides → assemble_mrna_constructs
    must produce a non-empty construct. Pre-#253 the construct
    writers were unreachable from external input — this test pins
    that the path is now wired."""
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)
    assert ranked, "Expected non-empty ranked list from LENS fixture"

    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        utr_3p='HBB', poly_a_length=10,
        antigens_per_construct=3, max_constructs=1,
        max_antigen_length_aa=14,
        optimize_linkers=False,
    )
    constructs = assemble_mrna_constructs(ranked, options=options)
    assert len(constructs) == 1
    c = constructs[0]
    # Three antigens packed into one construct
    assert len(c.antigen_names) == 3
    # cds_aa contains all three antigen AA strings (or windowed
    # subwindows centered on the mutation)
    for vp_pair, name in zip(ranked, c.antigen_names):
        # The antigen name carries the gene_name; LENS fixture genes
        # are TP53, BRAF, KRAS.
        assert any(g in name for g in ("TP53", "BRAF", "KRAS")), name
    assert c.cds_nt.endswith("TAA")
    assert c.full_nt.endswith("A" * 10)


def test_lens_drives_peptide_construct_assembly_end_to_end():
    """LENS file → peptide-vaccine constructs (SLP mode)."""
    from vaxrank.peptide import (
        PeptideConstructConfig, assemble_peptide_constructs)

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)

    options = PeptideConstructConfig(mode='slp')
    constructs = assemble_peptide_constructs(ranked, options=options)
    assert len(constructs) == 3, (
        "SLP mode: one construct per ranked vaccine peptide; "
        "got %d" % len(constructs))


# ---- pVACseq path coverage ----------------------------------------------

def test_pvacseq_to_ranked_vaccine_peptides_round_trip():
    """Smoke test for the pVACseq path: fixture has 3 unique variants
    (TP53/BRAF/KRAS) with dashed-form IDs (chr-start-end-ref-alt).
    The path was previously degenerate — _PVACSEQ_IC50_COLS was
    untested and the ID parser only handled dotted form, producing
    an empty ranked list end-to-end.
    """
    from vaxrank.epitope_io import load_pvacseq
    from vaxrank.external_input import ranked_from_pvacseq_predictions

    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    report_df, predictions = load_pvacseq(path)
    ranked = ranked_from_pvacseq_predictions(predictions, path)
    assert len(ranked) == 3, (
        "pVACseq fixture has 3 unique variants; got %d" % len(ranked))
    genes = sorted(
        peps[0].mutant_protein_fragment.gene_name for _, peps in ranked)
    assert genes == ['BRAF', 'KRAS', 'TP53']
    # IDs are dashed-form (chr1-100000-100001-A-T), parsed correctly
    contigs = sorted(v.contig for v, _ in ranked)
    assert contigs == ['chr1', 'chr2', 'chr3']


def test_pvacseq_id_parser_handles_dashed_and_dotted():
    """The parser accepts both modern (chr1-100000-100001-A-T) and
    legacy (1.123.A.T) ID forms, plus the 4-part dashed variant
    (chr1-100000-A-T) without an explicit end position."""
    from vaxrank.external_input import _parse_pvacseq_id

    assert _parse_pvacseq_id("chr1-100000-100001-A-T") == ("chr1", 100000, "A", "T")
    assert _parse_pvacseq_id("chr1-100000-A-T") == ("chr1", 100000, "A", "T")
    assert _parse_pvacseq_id("1.123.A.T") == ("1", 123, "A", "T")
    # Garbage / wrong shape returns None instead of raising
    assert _parse_pvacseq_id("garbage") is None
    assert _parse_pvacseq_id("chr1-notapos-A-T") is None
    assert _parse_pvacseq_id(None) is None
    assert _parse_pvacseq_id("") is None


def test_pvacseq_drives_mrna_construct_assembly_end_to_end():
    """End-to-end: pVACseq report → ranked → mRNA constructs (mirror
    of the LENS end-to-end test). Pre-fix this produced zero
    constructs because the ranked list was empty."""
    from vaxrank.epitope_io import load_pvacseq
    from vaxrank.external_input import ranked_from_pvacseq_predictions
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _, predictions = load_pvacseq(path)
    ranked = ranked_from_pvacseq_predictions(predictions, path)
    assert ranked

    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        utr_3p='HBB', poly_a_length=10,
        antigens_per_construct=3, max_constructs=1,
        max_antigen_length_aa=12, min_antigen_length_aa=5,
        optimize_linkers=False)
    constructs = assemble_mrna_constructs(ranked, options=options)
    assert len(constructs) == 1
    c = constructs[0]
    assert len(c.antigen_names) == 3
    assert c.cds_nt.endswith("TAA")


# ---- LENS scoring-column edge cases --------------------------------------

def test_lens_warns_when_no_affinity_columns_detected(tmp_path, caplog):
    """When no pMHC_affinity predictor is detected, the
    representative-peptide pick is degenerate (every row ties at
    (2, 0.0); file-order first row "wins"). Pin that we warn
    instead of silently picking arbitrarily.

    Today this can be triggered by a stability-only LENS file
    (netmhcstabpan present without mhcflurry/netmhcpan). It can also
    fire if LENS ever adds tools whose ``kind`` is something other
    than ``pMHC_affinity`` (presentation-only, cleavage,
    immunogenicity, etc.) and the registry classifies them
    accordingly. The warning is kind-agnostic — it reports the
    detected kinds in its message.
    """
    import logging
    from vaxrank.epitope_io import load_lens
    from vaxrank.external_input import ranked_from_lens_predictions

    # Stability-only fixture (the only "no affinity" case the
    # current registry can produce). If new predictor kinds are added
    # to _LENS_PREDICTOR_REGISTRY, add fixtures here covering them.
    no_aff = tmp_path / "no_affinity.tsv"
    no_aff.write_text(
        "allele\tpeptide\tnetmhcstabpan_1.0.halflife_hours\t"
        "netmhcstabpan_1.0.perc_rank_stab\tantigen_source\tmut_aa_pos\t"
        "variant_coords\tgene_name\ttpm\tpep_context\n"
        "HLA-A02:01\tSVVGSSSSS\t12.5\t0.5\tSNV\t245\t"
        "chr17:7675088:C:T\tTP53\t42.5\tAASVVGSSSSSGTR\n")

    _, predictions = load_lens(str(no_aff))
    with caplog.at_level(logging.WARNING):
        ranked = ranked_from_lens_predictions(predictions, str(no_aff))
    # Still produces output (the file-order representative)
    assert len(ranked) == 1
    # Warning was emitted, mentions the detected non-affinity kinds
    msg = next(
        (r.message for r in caplog.records
         if "no pMHC_affinity scoring columns" in r.message), None)
    assert msg is not None, \
        "Expected a warn-once for LENS input lacking affinity predictors"
    assert "pMHC_stability" in msg, \
        "Warning should report the detected non-affinity kinds; got %r" % msg


# ---- _emit_outputs observability -----------------------------------------

def test_emit_outputs_logs_active_and_fired_dispatch(caplog, tmp_path):
    """The dispatch summary log line is part of the contract — pins
    which vaccine types were considered active and which actually fired.
    """
    import logging
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)

    out_dir = str(tmp_path / "mrna_out")
    args = SimpleNamespace(
        output_csv='', output_xlsx_report='',
        output_neoepitope_report='',
        output_peptide='',
        output_mrna=out_dir,
        output_mrna_manifest='', output_mrna_csv='',
        output_mrna_csv_full_rows=True,
        mrna_signal_peptide=None, mrna_linker='(G4S)2',
        mrna_include_mitd=False, mrna_mitd='HLA_A',
        mrna_5p_utr='HBB', mrna_3p_utr='HBB',
        mrna_codon_species='h_sapiens',
        mrna_codon_method='use_best_codon',
        mrna_min_antigen_length_aa=5, mrna_max_antigen_length_aa=14,
        mrna_antigens_per_construct=3, mrna_max_constructs=1,
        mrna_candidates_per_slot=1, mrna_max_length_nt=4000,
        mrna_optimize_linkers=False, mrna_junction_candidates='',
        mrna_junction_rank_strong=0.5, mrna_junction_rank_mild=2.0,
        mrna_poly_a_length=20, mrna_poly_a_segmented=False,
        mrna_poly_a_first_segment=30,
        mrna_poly_a_segment_linker='GCATATGACT',
        vaccine_type=['mrna'],
    )
    with caplog.at_level(logging.INFO):
        _emit_outputs(args, ranked, source='external')
    msgs = [r.message for r in caplog.records]
    dispatch = [m for m in msgs if 'Vaccine-type dispatch' in m]
    assert len(dispatch) == 1
    line = dispatch[0]
    assert "[external]" in line
    assert "active=['mrna']" in line
    assert "wrote=['mrna']" in line


def test_resolve_vaccine_types_default_is_peptide():
    """Vaxrank is a vaccine-ranking library: there's always ranking
    to do. The default vaccine type is peptide; ``--vaccine-type`` is
    multi-valued so future types (DNA, etc.) plug in cleanly."""
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_vaccine_types

    args = SimpleNamespace(vaccine_type=['peptide'])
    assert _resolve_vaccine_types(args) == {'peptide'}

    # Default applied when attr is missing or empty
    args = SimpleNamespace()
    assert _resolve_vaccine_types(args) == {'peptide'}
    args = SimpleNamespace(vaccine_type=[])
    assert _resolve_vaccine_types(args) == {'peptide'}


def test_resolve_vaccine_types_explicit_lists():
    """Multi-valued list returned as a set; order-free; duplicates
    collapse."""
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_vaccine_types

    cases = [
        (['peptide'], {'peptide'}),
        (['mrna'], {'mrna'}),
        (['peptide', 'mrna'], {'peptide', 'mrna'}),
        (['mrna', 'peptide'], {'peptide', 'mrna'}),
        (['mrna', 'mrna'], {'mrna'}),
    ]
    for vtypes, expected in cases:
        args = SimpleNamespace(vaccine_type=vtypes)
        assert _resolve_vaccine_types(args) == expected


def test_resolve_vaccine_types_rejects_unknown():
    """'dna' or any unknown type must raise; users shouldn't get
    silent no-ops while writing CLI scripts that pre-declare future
    types."""
    import pytest as _pytest
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_vaccine_types

    args = SimpleNamespace(vaccine_type=['dna'])
    with _pytest.raises(ValueError, match="Unknown vaccine type"):
        _resolve_vaccine_types(args)


def test_resolve_vaccine_types_bare_string_tolerated():
    """argparse always yields a list, but a downstream caller might
    pass a bare string. Coerce instead of crashing."""
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_vaccine_types

    args = SimpleNamespace(vaccine_type='mrna')
    assert _resolve_vaccine_types(args) == {'mrna'}


def test_emit_outputs_warns_when_output_path_lacks_matching_type(
        caplog, tmp_path):
    """If user passes --output-mrna but doesn't add 'mrna' to
    --vaccine-type, we should warn loudly so they don't think mRNA
    is being written."""
    import logging
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _emit_outputs

    args = SimpleNamespace(
        vaccine_type=['peptide'],
        output_peptide='',
        output_mrna=str(tmp_path / "mrna_out"),
        output_csv='', output_xlsx_report='',
        output_neoepitope_report='',
    )
    with caplog.at_level(logging.WARNING):
        _emit_outputs(args, ranked=[], source='external')
    assert any(
        "--output-mrna" in r.message
        and "not in --vaccine-type" in r.message
        for r in caplog.records), \
        "Expected mismatch warning when --output-mrna is set but "\
        "vaccine-type=['peptide']"


def test_lens_with_vaccine_type_mrna_writes_three_fastas(tmp_path):
    """End-to-end: --input-lens + --output-mrna writes three FASTAs.
    Closes the pre-#253 short-circuit gap that made this fail."""
    from types import SimpleNamespace

    from vaxrank.cli.entry_point import _emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)

    out_dir = str(tmp_path / "mrna_out")
    args = SimpleNamespace(
        # report flags off
        output_csv='', output_xlsx_report='',
        output_neoepitope_report='',
        # peptide modality off
        output_peptide='',
        # mRNA modality on
        output_mrna=out_dir,
        output_mrna_manifest=str(tmp_path / "manifest.json"),
        output_mrna_csv='',
        output_mrna_csv_full_rows=True,
        # mRNA defaults
        mrna_signal_peptide=None,
        mrna_linker='(G4S)2',
        mrna_include_mitd=False,
        mrna_mitd='HLA_A',
        mrna_5p_utr='HBB', mrna_3p_utr='HBB',
        mrna_codon_species='h_sapiens',
        mrna_codon_method='use_best_codon',
        mrna_min_antigen_length_aa=5,
        mrna_max_antigen_length_aa=14,
        mrna_antigens_per_construct=3,
        mrna_max_constructs=1,
        mrna_candidates_per_slot=1,
        mrna_max_length_nt=4000,
        mrna_optimize_linkers=False,
        mrna_junction_candidates='',
        mrna_junction_rank_strong=0.5,
        mrna_junction_rank_mild=2.0,
        mrna_poly_a_length=20,
        mrna_poly_a_segmented=False,
        mrna_poly_a_first_segment=30,
        mrna_poly_a_segment_linker='GCATATGACT',
        # vaccine-type switch (multi-valued; pass a list)
        vaccine_type=['mrna'],
    )
    _emit_outputs(args, ranked, source='external')
    assert os.path.isfile(os.path.join(out_dir, "cds.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "no_polyA.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "full.fasta"))
    with open(args.output_mrna_manifest) as f:
        manifest = json.load(f)
    assert len(manifest) == 1
    assert manifest[0]['modality'] == 'mrna'
    # full ends with polyA
    full_body = ''.join(
        line for line in open(
            os.path.join(out_dir, "full.fasta")).read().split("\n")
        if not line.startswith(">"))
    assert full_body.endswith("A" * 20)


def test_output_path_without_matching_type_does_not_write(tmp_path):
    """Default vaccine-type=['peptide']. Passing --output-mrna without
    'mrna' in vaccine-type should NOT write mRNA constructs (the
    paired warning is covered separately)."""
    from types import SimpleNamespace

    from vaxrank.cli.entry_point import _emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked = ranked_from_lens_predictions(predictions, path)

    out_dir = str(tmp_path / "mrna_out")
    args = SimpleNamespace(
        output_csv='', output_xlsx_report='',
        output_neoepitope_report='',
        output_peptide='',  # peptide writer also no-ops (no path)
        output_mrna=out_dir,
        output_mrna_manifest='', output_mrna_csv='',
        output_mrna_csv_full_rows=True,
        mrna_signal_peptide=None, mrna_linker='(G4S)2',
        mrna_include_mitd=False, mrna_mitd='HLA_A',
        mrna_5p_utr='HBB', mrna_3p_utr='HBB',
        mrna_codon_species='h_sapiens',
        mrna_codon_method='use_best_codon',
        mrna_min_antigen_length_aa=5, mrna_max_antigen_length_aa=14,
        mrna_antigens_per_construct=3, mrna_max_constructs=1,
        mrna_candidates_per_slot=1, mrna_max_length_nt=4000,
        mrna_optimize_linkers=False, mrna_junction_candidates='',
        mrna_junction_rank_strong=0.5, mrna_junction_rank_mild=2.0,
        mrna_poly_a_length=20, mrna_poly_a_segmented=False,
        mrna_poly_a_first_segment=30,
        mrna_poly_a_segment_linker='GCATATGACT',
        vaccine_type=['peptide'],
    )
    _emit_outputs(args, ranked, source='external')
    assert not os.path.exists(out_dir), \
        "vaccine-type didn't include 'mrna'; mRNA writer should "\
        "not have run, but output dir was created"
