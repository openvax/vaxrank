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

"""Tests for mRNA vaccine construct assembly."""

import json
import os
import tempfile
from types import SimpleNamespace

import pytest
from varcode import Variant

from vaxrank.mrna import (
    RNAConstructConfig,
    assemble_mrna_constructs,
    codon_optimize,
    summarize_mrna_ranking_decisions,
    write_mrna_outputs,
)
from vaxrank.mrna_library import (
    MITD_HLA_A,
    SIGNAL_PEPTIDE_TPA,
    UTR_3P_HBB,
    UTR_5P_HBB,
)
from vaxrank.vaccine_library import get_linker

# (G4S)3 was a static entry in 2.12.x; in 2.13+ it is parsed from the
# compositional grammar. Resolve via get_linker so the test references
# survive the rename.
LINKER_G4S3 = get_linker("(G4S)3")


CODON_TABLE = {
    'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC',
    'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCC', 'Q': 'CAG', 'R': 'CGG',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAC',
}


def _translate(dna):
    """Minimal DNA -> protein for assertions."""
    table = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'TGC': 'C', 'TGT': 'C',
        'GAC': 'D', 'GAT': 'D',
        'GAA': 'E', 'GAG': 'E',
        'TTC': 'F', 'TTT': 'F',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'CAC': 'H', 'CAT': 'H',
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
        'AAA': 'K', 'AAG': 'K',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L',
        'ATG': 'M',
        'AAC': 'N', 'AAT': 'N',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAA': 'Q', 'CAG': 'Q',
        'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'TGG': 'W',
        'TAC': 'Y', 'TAT': 'Y',
        'TAA': '*', 'TAG': '*', 'TGA': '*',
    }
    return ''.join(table.get(dna[i:i + 3], 'X') for i in range(0, len(dna) - 2, 3))


def _peptide_stub(amino_acids, gene_name='GENE'):
    fragment = SimpleNamespace(amino_acids=amino_acids, gene_name=gene_name)
    return SimpleNamespace(mutant_protein_fragment=fragment)


def _variant_pair(amino_acids, contig='1', start=1000, gene_name='GENE'):
    variant = Variant(contig, start, 'A', 'T')
    return (variant, [_peptide_stub(amino_acids, gene_name=gene_name)])


# ---- mRNA ranking-decisions summary (#270) ------------------------------


def _vp_with_score(amino_acids, gene_name, combined_score):
    """Like ``_peptide_stub`` but with a settable ``combined_score``
    so we can verify the ordered list comes through."""
    fragment = SimpleNamespace(amino_acids=amino_acids, gene_name=gene_name)
    return SimpleNamespace(
        mutant_protein_fragment=fragment, combined_score=combined_score)


def test_summarize_mrna_ranking_decisions_splits_at_cap():
    """``selected`` = top ``antigens_per_construct × max_constructs``
    of the ranked variants; ``dropped`` = the rest. ``cap`` and
    ``total_ranked`` are surfaced for the report header."""
    pairs = [
        (Variant('1', 100 + i, 'A', 'T'),
         [_vp_with_score("KLQGHSAPVLDVIVN", "GENE%d" % i, 100.0 - i)])
        for i in range(7)
    ]
    options = RNAConstructConfig(
        antigens_per_construct=2, max_constructs=2)  # cap = 4
    summary = summarize_mrna_ranking_decisions(pairs, options)
    assert summary['cap'] == 4
    assert summary['total_ranked'] == 7
    assert summary['antigens_per_construct'] == 2
    assert summary['max_constructs'] == 2
    assert len(summary['selected']) == 4
    assert len(summary['dropped']) == 3
    # Selected entries are the first 4 in ranked order with rank 1..4.
    assert [a['rank'] for a in summary['selected']] == [1, 2, 3, 4]
    assert [a['gene_name'] for a in summary['selected']] == \
        ['GENE0', 'GENE1', 'GENE2', 'GENE3']
    # Dropped entries continue the ranking from rank 5.
    assert [a['rank'] for a in summary['dropped']] == [5, 6, 7]


def test_summarize_mrna_ranking_decisions_no_drops_when_under_cap():
    """When the ranked list is shorter than the cap, every variant
    is selected and ``dropped`` is empty."""
    pairs = [
        (Variant('1', 100, 'A', 'T'),
         [_vp_with_score("KLQGHSAPVLDVIVN", "GENE_A", 50.0)]),
        (Variant('1', 200, 'A', 'T'),
         [_vp_with_score("MNNVDEILGRWESPV", "GENE_B", 30.0)]),
    ]
    options = RNAConstructConfig(
        antigens_per_construct=5, max_constructs=2)  # cap = 10
    summary = summarize_mrna_ranking_decisions(pairs, options)
    assert summary['cap'] == 10
    assert summary['total_ranked'] == 2
    assert len(summary['selected']) == 2
    assert summary['dropped'] == []


def test_summarize_mrna_ranking_decisions_empty_input():
    """No ranked variants → both lists empty, cap still surfaces so
    the report can still render the header."""
    options = RNAConstructConfig(
        antigens_per_construct=5, max_constructs=2)
    summary = summarize_mrna_ranking_decisions([], options)
    assert summary['cap'] == 10
    assert summary['total_ranked'] == 0
    assert summary['selected'] == []
    assert summary['dropped'] == []


def test_ascii_report_renders_mrna_ranking_section():
    """End-to-end: the ASCII report includes a ``mRNA vaccine
    antigen selection`` section when the
    ``mrna_ranking`` template-data field is populated."""
    from vaxrank.report import _make_report
    import io

    template_data = {
        'patient_info': {'Patient ID': 'TEST'},
        'package_versions': {},
        'reviewers': [], 'final_review': '',
        'input_json_file': None,
        'include_manufacturability': False,
        'include_wt_epitopes': True,
        'args': [],
        'variants': [],
        'vaccine_constructions': {
            'mrna': {
                'antigens_per_construct': 5,
                'max_constructs': 2,
                'cap': 10,
                'total_ranked': 12,
                'selected': [
                    {'rank': 1, 'gene_name': 'GENEA',
                     'description': 'GENEA_1_100_A_T',
                     'combined_score': 12.5,
                     'allele_tiers': {}, 'manufacturability': {}},
                    {'rank': 2, 'gene_name': 'GENEB',
                     'description': 'GENEB_2_200_A_T',
                     'combined_score': 11.1,
                     'allele_tiers': {}, 'manufacturability': {}},
                ],
                'dropped': [
                    {'rank': 11, 'gene_name': 'GENEK',
                     'description': 'GENEK_3_300_A_T',
                     'combined_score': 4.3,
                     'allele_tiers': {}, 'manufacturability': {}},
                ],
                'coverage': [],
            },
        },
    }
    f = io.StringIO()
    _make_report(template_data, f, 'templates/template.txt')
    rendered = f.getvalue()
    # Post-2.25 section header (was "mRNA vaccine antigen selection").
    assert 'Vaccine construction — mrna' in rendered
    assert '2 of 12 ranked' in rendered
    assert 'GENEA_1_100_A_T' in rendered
    assert 'GENEB_2_200_A_T' in rendered
    assert 'Not selected' in rendered
    assert 'GENEK_3_300_A_T' in rendered
    # Cap params named in the header.
    assert '5' in rendered  # antigens_per_construct
    assert '10' in rendered  # cap


def test_ascii_report_skips_construction_section_when_no_modalities():
    """When ``vaccine_constructions`` is empty (no active modality
    or no ranked antigens), no construction sections render."""
    from vaxrank.report import _make_report
    import io
    template_data = {
        'patient_info': {'Patient ID': 'TEST'},
        'package_versions': {},
        'reviewers': [], 'final_review': '',
        'input_json_file': None,
        'include_manufacturability': False,
        'include_wt_epitopes': True,
        'args': [], 'variants': [],
        'vaccine_constructions': {},
        'mrna_ranking': None,
    }
    f = io.StringIO()
    _make_report(template_data, f, 'templates/template.txt')
    rendered = f.getvalue()
    assert 'Vaccine construction' not in rendered


def test_resolve_named_accepts_dashes_and_case_variants():
    """Library keys are written ``HLA_B`` / ``tPA`` / ``IgK`` etc.,
    but users in the field write ``HLA-B`` / ``hla-b`` / ``HLAB``.
    ``_resolve_named`` strips dashes and underscores and lowercases
    so all of those resolve to the same entry."""
    from vaxrank.mrna_library import SIGNAL_PEPTIDES, MITDS, UTRS_5P
    from vaxrank.mrna import _resolve_named

    canonical = _resolve_named(SIGNAL_PEPTIDES, 'HLA_B', 'signal_peptide')
    for variant in ['HLA-B', 'hla_b', 'hla-b', 'HLAB', 'hlab']:
        assert _resolve_named(
            SIGNAL_PEPTIDES, variant, 'signal_peptide') == canonical, variant
    # MITD
    a_canonical = _resolve_named(MITDS, 'HLA_A', 'mitd')
    assert _resolve_named(MITDS, 'HLA-A', 'mitd') == a_canonical
    assert _resolve_named(MITDS, 'hla-a', 'mitd') == a_canonical
    # UTR
    hbb = _resolve_named(UTRS_5P, 'HBB', '5p_utr')
    assert _resolve_named(UTRS_5P, 'hbb', '5p_utr') == hbb
    # Unknown still raises
    import pytest
    with pytest.raises(ValueError, match="Unknown signal_peptide"):
        _resolve_named(SIGNAL_PEPTIDES, 'NOT_A_REAL_KEY', 'signal_peptide')


def test_codon_optimize_preserves_translation():
    aa = "MGSLLPVTLLT"
    dna = codon_optimize(aa, species='h_sapiens', method='use_best_codon')
    assert len(dna) == len(aa) * 3
    assert _translate(dna) == aa


def test_codon_optimize_empty():
    assert codon_optimize("") == ""


def test_assembly_basic_construct():
    pairs = [
        _variant_pair("KLQGHSAPVLDVIVN", contig='1', start=100, gene_name='GENEA'),
        _variant_pair("MNNVDEILGRWESPV", contig='2', start=200, gene_name='GENEB'),
    ]
    options = RNAConstructConfig(signal_peptide='tPA', linker='GS3',
                          include_mitd=False, utr_3p='HBB')
    constructs = assemble_mrna_constructs(pairs, options=options)
    assert len(constructs) == 1
    c = constructs[0]
    assert c.antigen_names == ['GENEA_1_100_A_T', 'GENEB_2_200_A_T']
    assert c.no_polya_nt.startswith(UTR_5P_HBB)
    assert c.no_polya_nt.endswith(UTR_3P_HBB)
    # full_nt has polyA tail beyond UTR_3P_HBB
    assert c.full_nt.startswith(UTR_5P_HBB)
    assert c.full_nt.endswith("A" * 120)
    # CDS sits between UTRs and ends with a stop codon and starts with ATG
    cds = c.no_polya_nt[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    assert cds[:3] == 'ATG'
    assert cds[-3:] == 'TAA'
    protein = _translate(cds[:-3])
    expected = SIGNAL_PEPTIDE_TPA + LINKER_G4S3.amino_acids.join(['KLQGHSAPVLDVIVN', 'MNNVDEILGRWESPV'])
    assert protein == expected


def test_assembly_prepends_atg_when_no_signal_peptide():
    # Antigens that don't naturally start with M should still produce a
    # translatable CDS — assembler prepends an M.
    pairs = [_variant_pair("KLQGHSAPVLDVIVN", gene_name='GENE')]
    options = RNAConstructConfig(signal_peptide=None, linker='GS3', include_mitd=False, utr_3p='HBB')
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.no_polya_nt[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    assert cds[:3] == 'ATG'
    protein = _translate(cds[:-3])
    assert protein == "MKLQGHSAPVLDVIVN"


def test_assembly_does_not_double_prepend_when_antigen_starts_with_m():
    pairs = [_variant_pair("MNNVDEILGRWESPV", gene_name='GENE')]
    options = RNAConstructConfig(signal_peptide=None, linker='GS3', include_mitd=False, utr_3p='HBB')
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.no_polya_nt[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    protein = _translate(cds[:-3])
    assert protein == "MNNVDEILGRWESPV"


def test_assembly_with_mitd_appends_trafficking_domain():
    pairs = [_variant_pair("KLQGHSAPVLDVIVN")]
    options = RNAConstructConfig(signal_peptide=None, linker='GS3',
                          include_mitd=True, mitd='HLA_A', utr_3p='HBB')
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.no_polya_nt[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    protein = _translate(cds[:-3])
    # No signal peptide selected and antigen doesn't start with M, so the
    # assembler prepends one to keep the CDS translatable.
    expected = "M" + "KLQGHSAPVLDVIVN" + LINKER_G4S3.amino_acids + MITD_HLA_A
    assert protein == expected


def test_assembly_splits_on_max_antigens():
    pairs = [_variant_pair("AAAA", start=i) for i in range(5)]
    options = RNAConstructConfig(signal_peptide=None, linker='GS',
                          include_mitd=False, antigens_per_construct=2,
                          max_constructs=10, max_antigen_length_aa=4,
                          utr_3p='HBB')
    constructs = assemble_mrna_constructs(pairs, options=options)
    assert [len(c.antigen_names) for c in constructs] == [2, 2, 1]


def test_assembly_splits_on_max_length():
    long_aa = "A" * 200
    pairs = [_variant_pair(long_aa, start=i) for i in range(4)]
    options = RNAConstructConfig(
        signal_peptide=None, linker='GS', include_mitd=False,
        antigens_per_construct=10, max_constructs=10,
        max_antigen_length_aa=200, utr_3p='HBB',
        max_length_nt=len(UTR_5P_HBB) + len(UTR_3P_HBB) + 3 + 600 * 2)
    constructs = assemble_mrna_constructs(pairs, options=options)
    assert len(constructs) > 1


def test_assembly_max_constructs_caps_output():
    # New default max_constructs=1: extra antigens spill into a second
    # construct only when explicitly allowed.
    pairs = [_variant_pair("A" * 30, start=i) for i in range(3)]
    options = RNAConstructConfig(
        signal_peptide=None, linker='GS', include_mitd=False,
        antigens_per_construct=1, max_constructs=1,
        max_antigen_length_aa=30, utr_3p='HBB')
    [c] = assemble_mrna_constructs(pairs, options=options)
    assert len(c.antigen_names) == 1


def test_assembly_no_antigens_returns_empty():
    assert assemble_mrna_constructs([], options=RNAConstructConfig()) == []


def test_unknown_signal_peptide_raises():
    pairs = [_variant_pair("KLQ")]
    with pytest.raises(ValueError):
        assemble_mrna_constructs(pairs, options=RNAConstructConfig(signal_peptide='nope'))


def test_assembly_p2a_linker_preserves_blessed_dna():
    # When a 2A linker is selected, the blessed DNA from the literature
    # must appear verbatim in the optimized CDS.
    from vaxrank.vaccine_library import LINKER_P2A
    pairs = [
        _variant_pair("KLQGHSAPVLDVIVN", contig='1', start=100, gene_name='GENEA'),
        _variant_pair("MNNVDEILGRWESPV", contig='2', start=200, gene_name='GENEB'),
    ]
    options = RNAConstructConfig(signal_peptide='tPA', linker='P2A',
                          include_mitd=False, utr_3p='HBB')
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.no_polya_nt[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    assert LINKER_P2A.dna in cds
    # And translation matches the AA-level expectation (signal + a1 + P2A + a2)
    protein = _translate(cds[:-3])
    expected = (
        SIGNAL_PEPTIDE_TPA
        + LINKER_P2A.amino_acids.join(['KLQGHSAPVLDVIVN', 'MNNVDEILGRWESPV'])
    )
    assert protein == expected


def test_write_mrna_outputs_fasta_and_manifest():
    pairs = [_variant_pair("KLQGHSAP")]
    constructs = assemble_mrna_constructs(
        pairs,
        options=RNAConstructConfig(signal_peptide=None, include_mitd=False,
                                  utr_3p='HBB', poly_a_length=0))
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        manifest_path = os.path.join(tmp, "out.json")
        write_mrna_outputs(constructs, out_dir, manifest_path)
        # Three FASTAs land in the directory.
        assert os.path.isfile(os.path.join(out_dir, "cds.fasta"))
        assert os.path.isfile(os.path.join(out_dir, "no_polyA.fasta"))
        assert os.path.isfile(os.path.join(out_dir, "full.fasta"))
        with open(os.path.join(out_dir, "no_polyA.fasta")) as f:
            text = f.read()
        assert text.startswith(">mrna_001")
        assert UTR_5P_HBB in text.replace('\n', '')
        with open(manifest_path) as f:
            manifest = json.load(f)
        entry = manifest[0]
        assert entry['modality'] == 'mrna'
        assert entry['length_unit'] == 'nt'
        assert entry['name'] == 'mrna_001'
        assert entry['antigen_names'] == ['GENE_1_1000_A_T']
        # New structured fields
        assert 'lengths' in entry
        assert 'cds' in entry
        assert 'no_polya_nt' in entry
        assert 'full_nt' in entry
        assert 'antigens' in entry
        assert 'elements' in entry
        # Back-compat 2.12 schema preserved
        assert 'components' in entry
        assert 'manufacturability' in entry


def test_mrna_fasta_uses_semicolon_antigen_separator():
    """The FASTA description's ``antigens=`` field must use ``;`` as
    the separator — comma collides with FASTA-parsing tools that
    treat ``,`` as a token boundary in description fields. Pin the
    contract so a future maintainer can't accidentally regress to
    comma-joined.
    """
    pairs = [
        _variant_pair("KLQGHSAPVLDVIVN", contig='1', start=100, gene_name='GENEA'),
        _variant_pair("MNNVDEILGRWESPV", contig='2', start=200, gene_name='GENEB'),
    ]
    constructs = assemble_mrna_constructs(
        pairs,
        options=RNAConstructConfig(
            signal_peptide=None, include_mitd=False,
            utr_3p='HBB', poly_a_length=0))
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        write_mrna_outputs(constructs, out_dir)
        with open(os.path.join(out_dir, "cds.fasta")) as f:
            header = f.readline().rstrip()
    assert "antigens=" in header
    # Both antigens land in one description, separated by ';'.
    assert "GENEA_1_100_A_T;GENEB_2_200_A_T" in header
    # And the legacy comma-joined form must not survive.
    assert "GENEA_1_100_A_T,GENEB_2_200_A_T" not in header


def test_mrna_max_constructs_logs_topk_selection_at_info(caplog):
    """When ``--mrna-max-constructs`` caps off the input list, the
    log line must (a) be INFO not WARNING (top-k selection is the
    whole point of vaxrank, not an error), and (b) name the kept
    count + total-ranked count so the operator can tell what fraction
    landed in the vaccine."""
    import logging
    # Build 6 antigens against a cap of 1 construct × 2 antigens =
    # 2 selected, 4 not selected.
    pairs = [
        _variant_pair("KLQGHSAPVLDVIVN", contig=str(i + 1),
                      start=1000 + i, gene_name="GENE%d" % i)
        for i in range(6)
    ]
    with caplog.at_level(logging.INFO):
        assemble_mrna_constructs(
            pairs,
            options=RNAConstructConfig(
                signal_peptide=None, include_mitd=False,
                utr_3p='HBB', poly_a_length=0,
                antigens_per_construct=2, max_constructs=1))
    # The cap event surfaces at INFO, never WARNING.
    cap_records = [
        r for r in caplog.records
        if 'mRNA assembly' in r.getMessage()
        and '--mrna-max-constructs' in r.getMessage()
    ]
    assert cap_records, \
        "Expected an mRNA assembly top-k log line at the cap"
    assert all(r.levelno == logging.INFO for r in cap_records), \
        "Cap event must be INFO, not WARNING (top-k selection is the goal)"
    msg = cap_records[0].getMessage()
    # Selected 2 of 6 ranked antigens.
    assert "2" in msg
    assert "6" in msg
    # The next-best (not-selected) antigen is named for traceability.
    assert "GENE2" in msg


def test_mrna_min_antigen_length_warning(caplog):
    # mRNA also enforces min_antigen_length_aa via warning. Regression
    # for the field being declared but unused (#246 review).
    import logging
    pairs = [_variant_pair("KLQGH", gene_name='G')]  # 5 aa
    with caplog.at_level(logging.WARNING):
        assemble_mrna_constructs(
            pairs,
            options=RNAConstructConfig(
                signal_peptide=None, include_mitd=False,
                min_antigen_length_aa=15))
    assert any("below" in r.message and "mrna-min-antigen-length" in r.message
               for r in caplog.records), \
        "Expected a min-antigen-length warning for an undersized mRNA antigen"


def test_select_antigen_window_used_by_both_modalities():
    # The mutation-centered window helper lives in vaccine_library and
    # is used by both peptide.py and mrna.py. Catch a regression where
    # one modality forks its own copy of the logic.
    from vaxrank import mrna, peptide
    from vaxrank.vaccine_library import select_antigen_window
    assert mrna.select_antigen_window is select_antigen_window
    assert peptide.select_antigen_window is select_antigen_window


# ---- polyA + structured manifest + CSV (issue #252) -----------------------

def test_build_poly_a_default_120():
    from vaxrank.mrna import build_poly_a
    options = RNAConstructConfig()
    tail = build_poly_a(options)
    assert tail == "A" * 120


def test_build_poly_a_zero_length():
    from vaxrank.mrna import build_poly_a
    options = RNAConstructConfig(poly_a_length=0)
    assert build_poly_a(options) == ""


def test_build_poly_a_segmented_bnt162b2_pattern():
    """BNT162b2 architecture: A30 + GCATATGACT + A70 = 100 A's split
    by a 10-nt linker (Xia 2021, PMC8310186)."""
    from vaxrank.mrna import build_poly_a
    options = RNAConstructConfig(
        poly_a_length=100, poly_a_segmented=True,
        poly_a_first_segment=30,
        poly_a_segment_linker="GCATATGACT")
    tail = build_poly_a(options)
    assert tail == "A" * 30 + "GCATATGACT" + "A" * 70


def test_build_poly_a_negative_raises():
    from vaxrank.mrna import build_poly_a
    options = RNAConstructConfig(poly_a_length=-1)
    with pytest.raises(ValueError):
        build_poly_a(options)


def test_assemble_full_includes_polya_no_polya_does_not():
    """``c.full_nt`` must end with the polyA tail; ``c.no_polya_nt``
    must end with the 3' UTR (no polyA). This is the contract for
    the three-FASTA output."""
    pairs = [_variant_pair("KLQGHSAPVL")]
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        utr_3p='HBB', poly_a_length=120)
    [c] = assemble_mrna_constructs(pairs, options=options)
    assert c.no_polya_nt.endswith(UTR_3P_HBB)
    assert not c.no_polya_nt.endswith("A" * 120)
    assert c.full_nt.endswith("A" * 120)
    # full = no_polya + polyA (exact concatenation, nothing else)
    assert c.full_nt == c.no_polya_nt + ("A" * 120)
    # cds_nt is no-UTR no-polyA, ends with stop codon
    assert c.cds_nt.endswith("TAA")
    assert UTR_5P_HBB not in c.cds_nt


def test_assemble_structured_elements_carry_aa_and_nt():
    """Each element in c.elements must carry both AA (where applicable)
    and nt, with lengths matching."""
    pairs = [_variant_pair("KLQGHSAPVL", gene_name='GENEA'),
             _variant_pair("MNNVDEILGRWESPV", contig='2', start=200,
                           gene_name='GENEB')]
    options = RNAConstructConfig(
        signal_peptide='HLA_A', linker='G4S',  # short linker for clear math
        optimize_linkers=False,
        include_mitd=True, mitd='HLA_A',
        utr_3p='HBB', poly_a_length=10,
        antigens_per_construct=2,
        max_antigen_length_aa=20)
    [c] = assemble_mrna_constructs(pairs, options=options)

    el = c.elements
    # 5' UTR + 3' UTR are nt-only (no AA)
    assert el['utr_5p']['nt'] == UTR_5P_HBB
    assert el['utr_5p']['length_nt'] == len(UTR_5P_HBB)
    assert el['utr_3p']['nt'] == UTR_3P_HBB
    # Signal peptide has aa + nt; nt length = aa length × 3
    sp = el['signal_peptide']
    assert sp['aa'].startswith("M")
    assert sp['length_nt'] == sp['length_aa'] * 3
    assert len(sp['nt']) == sp['length_nt']
    # MITD has aa + nt
    m = el['mitd']
    assert m is not None
    assert m['length_nt'] == m['length_aa'] * 3
    # PolyA element
    assert el['poly_a']['nt'] == "A" * 10
    assert el['poly_a']['length_nt'] == 10
    assert el['poly_a']['segmented'] is False
    # Antigens have AA + nt with length math
    assert len(c.antigens) == 2
    for a in c.antigens:
        assert a['length_nt'] == a['length_aa'] * 3
        assert len(a['nt']) == a['length_nt']
    # The CDS view stitches signal + antigens + linkers + MITD
    assert c.cds_aa.startswith(sp['aa'])
    # Stop codon present at end of cds_nt
    assert c.cds_nt.endswith("TAA")


def test_write_mrna_three_fastas_have_matching_records():
    """Each of the three FASTAs has the same number of records (one
    per construct), with matching names, but different sequences."""
    pairs = [_variant_pair("KLQGHSAPVL", gene_name='GENEA'),
             _variant_pair("MNNVDEILGR", contig='2', start=200,
                           gene_name='GENEB')]
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        utr_3p='HBB', poly_a_length=50,
        antigens_per_construct=2)
    constructs = assemble_mrna_constructs(pairs, options=options)
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        write_mrna_outputs(constructs, out_dir)
        cds = open(os.path.join(out_dir, "cds.fasta")).read()
        no_polya = open(os.path.join(out_dir, "no_polyA.fasta")).read()
        full = open(os.path.join(out_dir, "full.fasta")).read()
    # All three have the mrna_001 header
    for txt in (cds, no_polya, full):
        assert ">mrna_001" in txt
    # full > no_polyA > cds in body length
    cds_body = "".join(line for line in cds.split("\n") if not line.startswith(">"))
    no_pa_body = "".join(line for line in no_polya.split("\n") if not line.startswith(">"))
    full_body = "".join(line for line in full.split("\n") if not line.startswith(">"))
    assert len(cds_body) < len(no_pa_body) < len(full_body)
    # full = no_polyA + 50 A's
    assert full_body == no_pa_body + ("A" * 50)


def test_write_mrna_csv_one_row_per_element():
    """CSV is long-format: one row per (construct, element). All key
    element kinds appear; nt strings are non-empty for nt-bearing rows."""
    import csv as _csv

    pairs = [_variant_pair("KLQGHSAPVL", gene_name='GENEA')]
    options = RNAConstructConfig(
        signal_peptide='HLA_A', linker='G4S',
        optimize_linkers=False,
        include_mitd=True, mitd='HLA_A',
        utr_3p='HBB', poly_a_length=20,
        antigens_per_construct=1,
        max_antigen_length_aa=20)
    constructs = assemble_mrna_constructs(pairs, options=options)
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        csv_path = os.path.join(tmp, "manifest.csv")
        write_mrna_outputs(constructs, out_dir, csv_path=csv_path)
        with open(csv_path) as f:
            rows = list(_csv.DictReader(f))
    kinds = {r['element_kind'] for r in rows}
    expected = {
        'utr_5p', 'signal_peptide', 'antigen', 'mitd',
        'stop_codon', 'utr_3p', 'poly_a',
        'cds', 'no_polyA', 'full',
    }
    assert expected <= kinds, "missing element kinds: %s" % (expected - kinds)
    # Per-construct: each row's nt column is non-empty for nt-bearing kinds.
    for r in rows:
        if r['element_kind'] in {'utr_5p', 'utr_3p', 'poly_a',
                                  'stop_codon', 'no_polyA', 'full',
                                  'antigen', 'signal_peptide', 'mitd', 'cds'}:
            # poly_a_length=20 + 0-length linkers shouldn't have empty
            # nt rows for these element kinds.
            assert r['nt'], (
                "Expected nt for %s row, got empty: %r"
                % (r['element_kind'], r))


def test_polya_segmented_appears_in_full_fasta_and_csv():
    """End-to-end: --mrna-poly-a-segmented round-trips through the
    construct, FASTA, and CSV."""
    import csv as _csv

    pairs = [_variant_pair("KLQGHSAPVL", gene_name='G')]
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        utr_3p='HBB',
        poly_a_length=100, poly_a_segmented=True,
        poly_a_first_segment=30,
        poly_a_segment_linker="GCATATGACT",
        antigens_per_construct=1)
    constructs = assemble_mrna_constructs(pairs, options=options)
    [c] = constructs
    expected_tail = "A" * 30 + "GCATATGACT" + "A" * 70
    assert c.poly_a_nt == expected_tail
    assert c.full_nt.endswith(expected_tail)

    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        csv_path = os.path.join(tmp, "manifest.csv")
        write_mrna_outputs(constructs, out_dir, csv_path=csv_path)
        with open(csv_path) as f:
            rows = list(_csv.DictReader(f))
    poly_rows = [r for r in rows if r['element_kind'] == 'poly_a']
    assert len(poly_rows) == 1
    assert poly_rows[0]['nt'] == expected_tail
    assert poly_rows[0]['note'] == 'segmented'


# ---- Review-fix coverage (PR #253 second-pass) ---------------------------

def test_output_dir_rejects_existing_file(tmp_path):
    """Old API took a FASTA file path; new API takes a directory.
    Pointing --output-dir at an existing file must raise loudly
    when --vaccine-type includes mrna."""
    pairs = [_variant_pair("KLQGHSAP")]
    constructs = assemble_mrna_constructs(
        pairs, options=RNAConstructConfig(
            signal_peptide=None, include_mitd=False,
            utr_3p='HBB', poly_a_length=0))
    existing_file = tmp_path / "out.fasta"
    existing_file.write_text("placeholder")
    import pytest as _pytest
    with _pytest.raises(ValueError, match="directory"):
        write_mrna_outputs(constructs, str(existing_file))


def test_output_dir_rejects_fasta_suffix(tmp_path):
    """Even if the path doesn't exist, .fasta / .fa suffixes are blocked
    so 'vaxrank --vaccine-type=mrna --output-dir out.fasta' fails
    clearly instead of silently creating a directory called out.fasta/."""
    pairs = [_variant_pair("KLQGHSAP")]
    constructs = assemble_mrna_constructs(
        pairs, options=RNAConstructConfig(
            signal_peptide=None, include_mitd=False,
            utr_3p='HBB', poly_a_length=0))
    import pytest as _pytest
    for bad in ("out.fasta", "out.fa", "OUT.FASTA"):
        path = tmp_path / bad
        with _pytest.raises(ValueError, match="directory"):
            write_mrna_outputs(constructs, str(path))


def test_output_dir_idempotent_makedirs(tmp_path):
    """Writing twice into the same dir must succeed; second write
    overwrites the FASTAs in-place (no error from existing dir)."""
    pairs = [_variant_pair("KLQGHSAP")]
    constructs = assemble_mrna_constructs(
        pairs, options=RNAConstructConfig(
            signal_peptide=None, include_mitd=False,
            utr_3p='HBB', poly_a_length=0))
    out_dir = str(tmp_path / "mrna")
    write_mrna_outputs(constructs, out_dir)
    write_mrna_outputs(constructs, out_dir)
    assert os.path.isfile(os.path.join(out_dir, "cds.fasta"))


def test_csv_pre_mitd_linker_uses_index_label_not_index(tmp_path):
    """Pre-MITD linker rows put 'mitd' in index_label, not in the
    integer-typed index column. Spreadsheet sorts on `index` stay
    numeric."""
    import csv as _csv
    pairs = [_variant_pair("KLQGHSAPVL", gene_name='GENEA'),
             _variant_pair("MNNVDEILGRWESPV", contig='2', start=200,
                           gene_name='GENEB')]
    options = RNAConstructConfig(
        signal_peptide=None, linker='G4S',
        optimize_linkers=False,
        include_mitd=True, mitd='HLA_A',
        utr_3p='HBB', poly_a_length=0,
        antigens_per_construct=2,
        max_antigen_length_aa=20)
    constructs = assemble_mrna_constructs(pairs, options=options)
    csv_path = tmp_path / "manifest.csv"
    out_dir = tmp_path / "out"
    write_mrna_outputs(constructs, str(out_dir), csv_path=str(csv_path))
    with open(csv_path) as f:
        rows = list(_csv.DictReader(f))
    linker_rows = [r for r in rows if r['element_kind'] == 'linker']
    assert linker_rows, "Expected linker rows in CSV"
    pre_mitd = [r for r in linker_rows if r['index_label'] == 'mitd']
    inter = [r for r in linker_rows if r['index_label'] != 'mitd']
    assert len(pre_mitd) == 1, (
        "Expected exactly one pre-MITD linker row, got %d" % len(pre_mitd))
    assert pre_mitd[0]['index'] == '', (
        "Pre-MITD linker must leave the integer index blank "
        "(label-only); got %r" % pre_mitd[0]['index'])
    assert pre_mitd[0]['note'] == 'pre_mitd'
    # Inter-antigen linkers have integer indices, not labels.
    for r in inter:
        assert r['index'].isdigit(), (
            "Inter-antigen linker index must be numeric, got %r"
            % r['index'])
        assert r['index_label'] == ''


def test_csv_lean_mode_omits_full_rows(tmp_path):
    """csv_include_full_rows=False suppresses the cds / no_polyA /
    full summary rows so the CSV's widest cells are per-element."""
    import csv as _csv
    pairs = [_variant_pair("KLQGHSAP")]
    constructs = assemble_mrna_constructs(
        pairs, options=RNAConstructConfig(
            signal_peptide=None, include_mitd=False,
            utr_3p='HBB', poly_a_length=10))
    csv_path = tmp_path / "manifest.csv"
    out_dir = tmp_path / "out"
    write_mrna_outputs(
        constructs, str(out_dir),
        csv_path=str(csv_path), csv_include_full_rows=False)
    with open(csv_path) as f:
        rows = list(_csv.DictReader(f))
    kinds = {r['element_kind'] for r in rows}
    # full-rows are suppressed
    assert 'cds' not in kinds
    assert 'no_polyA' not in kinds
    assert 'full' not in kinds
    # per-element rows remain
    assert 'utr_5p' in kinds
    assert 'utr_3p' in kinds
    assert 'poly_a' in kinds


def test_dropped_dead_linker_name_param():
    """_build_protein_with_segments no longer takes a linker_name arg.
    Pin the signature so a future regression that re-adds the dead
    parameter is caught immediately."""
    import inspect
    from vaxrank.mrna import _build_protein_with_segments
    sig = inspect.signature(_build_protein_with_segments)
    assert 'linker_name' not in sig.parameters, (
        "linker_name was a dead parameter and should remain removed; "
        "got params %s" % list(sig.parameters))
    # signal_peptide_name and mitd_name are *used* and must stay.
    assert 'signal_peptide_name' in sig.parameters
    assert 'mitd_name' in sig.parameters


def test_start_codon_segment_has_named_kind():
    """The pre-2.14 code emitted name=None for the prepended start
    codon segment. Pin that we now use 'start_codon' as the name."""
    from vaxrank.mrna import _build_protein_with_segments
    from vaxrank.vaccine_library import get_linker
    linker = get_linker("G4S")
    # No signal peptide and antigen doesn't start with M → assembler
    # prepends a start codon segment.
    _, _, segments = _build_protein_with_segments(
        antigen_aas=["KLQGH"], antigen_names=["A"],
        signal_peptide_aa="", signal_peptide_name=None,
        linker=linker, mitd_aa="", mitd_name=None)
    start = next(s for s in segments if s['kind'] == 'start_codon')
    assert start['name'] == 'start_codon', (
        "start_codon segment name must not be None; got %r"
        % start['name'])


def test_mitd_nt_slice_is_bounded(tmp_path):
    """MITD nt slice must use [start_aa*3:end_aa*3], not an open-ended
    [start_aa*3:] slice. Pin by checking that nt length = aa length × 3
    exactly (open-ended slice would equal the rest of coding_dna,
    which today happens to match but is fragile)."""
    pairs = [_variant_pair("KLQGHSAPVL")]
    options = RNAConstructConfig(
        signal_peptide=None, linker='G4S',
        optimize_linkers=False,
        include_mitd=True, mitd='HLA_A',
        utr_3p='HBB', poly_a_length=0,
        antigens_per_construct=1,
        max_antigen_length_aa=15)
    [c] = assemble_mrna_constructs(pairs, options=options)
    m = c.elements['mitd']
    assert m is not None
    # Strict invariant: nt length matches AA length × 3 exactly.
    assert m['length_nt'] == m['length_aa'] * 3
    assert len(m['nt']) == m['length_aa'] * 3


def test_junction_swap_meta_appears_in_elements_when_optimizer_runs():
    """When the junction-aware optimizer runs (predictor + alleles
    supplied), elements['junction_swap'] is populated. Without
    predictor, elements['junction_swap'] reflects the fallback."""
    from types import SimpleNamespace
    from vaxrank.junction_swap import junction_kmers
    from varcode import Variant

    # Two antigens; predictor scores a strong hit on (G4S)2 only,
    # so AAA should win.
    a1, a2 = "KLQGHSAPVL", "DVIVNCDESLLAS"
    g4s2_kmers = junction_kmers(a1, "GGGGSGGGGS", a2, k_lengths=(9,))
    rank_table = {(g4s2_kmers[0], "HLA-A*02:01"): 0.05}

    class StubPredictor:
        def __init__(self, table): self.table = table
        def predict_peptides(self, peptides):
            out = []
            for p in peptides:
                out.append(SimpleNamespace(
                    peptide=p, allele="HLA-A*02:01",
                    percentile_rank=self.table.get(
                        (p, "HLA-A*02:01"), 99.0)))
            return out

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
        optimize_linkers=True,
        junction_swap_candidates=("(G4S)2", "AAA"),
        junction_kmer_lengths=(9,),
        antigens_per_construct=2, max_constructs=1,
        max_antigen_length_aa=20,
        utr_3p='HBB', poly_a_length=0,
    )
    [c] = assemble_mrna_constructs(
        pairs, options=options,
        mhc_predictor=StubPredictor(rank_table),
        mhc_alleles=["HLA-A*02:01"])
    # The structured elements dict must carry junction_swap
    assert 'junction_swap' in c.elements
    assert c.elements['junction_swap']['enabled'] is True
    assert 'AAA' in c.elements['junction_swap']['chosen']
