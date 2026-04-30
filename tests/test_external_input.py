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


def test_resolve_modality_auto_infers_from_output_flags():
    """The 'auto' sentinel reads which output paths are populated."""
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_modality

    cases = [
        (['auto'], 'out.fasta', '', {'peptide'}),
        (['auto'], '', 'out_dir', {'mrna'}),
        (['auto'], 'p.fa', 'm_dir', {'peptide', 'mrna'}),
        (['auto'], '', '', set()),
    ]
    for modality, op, om, expected in cases:
        args = SimpleNamespace(
            vaccine_modality=modality, output_peptide=op, output_mrna=om)
        assert _resolve_modality(args) == expected, (
            "modality=%r output_peptide=%r output_mrna=%r "
            "expected %r got %r" % (
                modality, op, om, expected, _resolve_modality(args)))


def test_resolve_modality_explicit_lists():
    """Concrete modality lists override 'auto' inference. Listed
    modalities are returned as a set; output flags are not consulted."""
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_modality

    cases = [
        (['peptide'], 'p.fa', 'm_dir', {'peptide'}),
        (['mrna'], 'p.fa', 'm_dir', {'mrna'}),
        (['peptide', 'mrna'], '', '', {'peptide', 'mrna'}),
        (['mrna', 'peptide'], '', '', {'peptide', 'mrna'}),  # order-free
        # duplicates collapse
        (['mrna', 'mrna'], '', '', {'mrna'}),
    ]
    for modality, op, om, expected in cases:
        args = SimpleNamespace(
            vaccine_modality=modality, output_peptide=op, output_mrna=om)
        assert _resolve_modality(args) == expected


def test_resolve_modality_none_sentinel():
    """'none' suppresses construct writers even if output paths set."""
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_modality

    args = SimpleNamespace(
        vaccine_modality=['none'],
        output_peptide='p.fa', output_mrna='m_dir')
    assert _resolve_modality(args) == set()


def test_resolve_modality_rejects_sentinel_mixed_with_concrete():
    """auto/none must be passed alone; mixing with mrna/peptide is a
    user error."""
    import pytest as _pytest
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_modality

    args = SimpleNamespace(
        vaccine_modality=['auto', 'mrna'],
        output_peptide='', output_mrna='m_dir')
    with _pytest.raises(ValueError, match="sentinel"):
        _resolve_modality(args)

    args = SimpleNamespace(
        vaccine_modality=['none', 'peptide'],
        output_peptide='', output_mrna='')
    with _pytest.raises(ValueError, match="sentinel"):
        _resolve_modality(args)


def test_resolve_modality_rejects_unknown_modality():
    """Adding 'dna' before its writer ships should fail clearly so
    users don't get a silent no-op."""
    import pytest as _pytest
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_modality

    args = SimpleNamespace(
        vaccine_modality=['dna'],
        output_peptide='', output_mrna='')
    with _pytest.raises(ValueError, match="Unknown vaccine modality"):
        _resolve_modality(args)


def test_resolve_modality_bare_string_tolerated():
    """argparse always yields a list, but a downstream caller might
    pass a bare string. Coerce instead of crashing."""
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_modality

    args = SimpleNamespace(
        vaccine_modality='mrna',
        output_peptide='', output_mrna='m_dir')
    assert _resolve_modality(args) == {'mrna'}


def test_lens_with_vaccine_modality_mrna_writes_three_fastas(tmp_path):
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
        # modality switch (multi-valued; pass a list)
        vaccine_modality=['mrna'],
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


def test_modality_none_skips_constructs_even_if_output_paths_set(tmp_path):
    """--vaccine-modality=none suppresses construct writers even if
    --output-peptide / --output-mrna are passed."""
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
        vaccine_modality=['none'],
    )
    _emit_outputs(args, ranked, source='external')
    assert not os.path.exists(out_dir), \
        "modality=none should suppress mRNA writer; output dir was created"
