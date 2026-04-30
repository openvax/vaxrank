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

"""Vaccine-design matrix: every combination of (vaccine_type,
antigen_content, antigens_per_construct).

The two design axes are orthogonal and shared across vaccine types:
- ``antigen_content`` ∈ {mutation_spanning, minimal_epitope}
- ``antigens_per_construct`` ∈ {1, N}

Combined with vaccine_type ∈ {peptide, mrna}, there are 8 design
points. This module pins each one end-to-end so a future refactor of
the dispatch / config plumbing can't silently drop a combination.

Eight distinct designs:

    +----------+-------------------+-----------------+-------------------------+
    | type     | antigen_content   | per_construct   | name                    |
    +----------+-------------------+-----------------+-------------------------+
    | peptide  | mutation_spanning | 1               | SLP (PGV-001 canonical) |
    | peptide  | mutation_spanning | N               | multi-SLP concatenated  |
    | peptide  | minimal_epitope   | 1               | minimal-ligand peptide  |
    | peptide  | minimal_epitope   | N               | concat-minimal-ligand   |
    | mrna     | mutation_spanning | N               | BioNTech FixVac / iNeST |
    | mrna     | mutation_spanning | 1               | single-antigen mRNA     |
    | mrna     | minimal_epitope   | N               | string-of-beads mRNA    |
    | mrna     | minimal_epitope   | 1               | single-ligand mRNA      |
    +----------+-------------------+-----------------+-------------------------+
"""

from types import SimpleNamespace

import pytest
from varcode import Variant

from vaxrank.epitope_prediction import EpitopePrediction
from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs


def _make_vaccine_peptide(
        amino_acids, gene_name, mut_start, mut_end, n_alt_reads,
        epitope_seqs_with_ic50):
    """Build a stub VaccinePeptide-shaped object that the assemblers
    can consume without going through the full ranking pipeline.

    ``epitope_seqs_with_ic50`` is a list of ``(peptide_sequence,
    ic50_nM)``. The first entry becomes the top mutant epitope.
    """
    fragment = SimpleNamespace(
        amino_acids=amino_acids,
        gene_name=gene_name,
        mutant_amino_acid_start_offset=mut_start,
        mutant_amino_acid_end_offset=mut_end,
        n_alt_reads=n_alt_reads,
    )
    epitope_predictions = [
        EpitopePrediction(
            allele="HLA-A*02:01",
            peptide_sequence=seq,
            wt_peptide_sequence="",
            ic50=ic50,
            wt_ic50=10000.0,
            percentile_rank=ic50 / 100.0,
            prediction_method_name="stub",
            overlaps_mutation=True,
            source_sequence=amino_acids,
            offset=0,
            occurs_in_reference=False,
        )
        for seq, ic50 in epitope_seqs_with_ic50
    ]
    return SimpleNamespace(
        mutant_protein_fragment=fragment,
        mutant_epitope_predictions=epitope_predictions,
        wildtype_epitope_predictions=[],
        epitope_predictions=epitope_predictions,
    )


def _ranked_three_variants():
    """Three variants × multiple ligand candidates each, suitable for
    every cell of the design matrix."""
    a1 = _make_vaccine_peptide(
        "AAKLQGHSAPVLDVIVNCD", gene_name="GENEA",
        mut_start=2, mut_end=12, n_alt_reads=10,
        epitope_seqs_with_ic50=[
            ("KLQGHSAPV", 30.0),    # top
            ("LQGHSAPVL", 80.0),    # second
            ("QGHSAPVLD", 220.0),   # third
        ])
    a2 = _make_vaccine_peptide(
        "MNNVDEILGRWESPVKLPK", gene_name="GENEB",
        mut_start=2, mut_end=12, n_alt_reads=8,
        epitope_seqs_with_ic50=[
            ("NVDEILGRW", 45.0),
            ("VDEILGRWE", 110.0),
        ])
    a3 = _make_vaccine_peptide(
        "AAVVVGADGVGKSALTIIQ", gene_name="GENEC",
        mut_start=4, mut_end=14, n_alt_reads=6,
        epitope_seqs_with_ic50=[
            ("VVVGADGVG", 65.0),
        ])
    return [
        (Variant("1", 100, "A", "T"), [a1]),
        (Variant("2", 200, "G", "C"), [a2]),
        (Variant("3", 300, "C", "T"), [a3]),
    ]


# ---------------------------------------------------------------------------
# Peptide design matrix (4 cells)
# ---------------------------------------------------------------------------

def test_peptide_slp():
    """SLP: mutation_spanning content, 1 antigen per construct.
    Canonical PGV-001 design — one synthetic long peptide per ranked
    variant."""
    options = PeptideConstructConfig(
        antigen_content="mutation_spanning",
        antigens_per_construct=1,
        max_antigen_length_aa=15,
        min_antigen_length_aa=5,
        max_constructs=10)
    constructs = assemble_peptide_constructs(_ranked_three_variants(), options)
    assert len(constructs) == 3, "SLP: one construct per variant"
    for c in constructs:
        assert len(c.antigen_names) == 1
        # Peptide is mutation-spanning (≤ 15 aa here)
        assert 5 <= len(c.sequence) <= 15
        assert c.components['antigen_content'] == "mutation_spanning"
        assert c.components['antigens_per_construct'] == 1


def test_peptide_multi_slp():
    """Multi-SLP: mutation_spanning, N per construct. Several SLPs
    concatenated with a linker into one longer peptide."""
    options = PeptideConstructConfig(
        mode="multi_epitope",  # legacy shorthand still works
        antigens_per_construct=3,
        linker="G4S3",
        max_antigen_length_aa=15,
        min_antigen_length_aa=5,
        max_constructs=10)
    constructs = assemble_peptide_constructs(_ranked_three_variants(), options)
    assert len(constructs) == 1
    [c] = constructs
    assert len(c.antigen_names) == 3
    # All three antigens linked
    assert c.sequence.count("GGGGS") >= 2  # at least 2 linkers between 3 antigens
    assert c.components['antigen_content'] == "mutation_spanning"
    assert c.components['antigens_per_construct'] == 3


def test_peptide_minimal_epitope():
    """Minimal-epitope peptide: single short MHC ligand per variant.
    Each construct is just the top mutant ligand (~9 aa)."""
    options = PeptideConstructConfig(
        antigen_content="minimal_epitope",
        antigens_per_construct=1,
        epitopes_per_antigen=1,
        min_antigen_length_aa=5,
        max_constructs=10)
    constructs = assemble_peptide_constructs(_ranked_three_variants(), options)
    assert len(constructs) == 3
    # Each construct's sequence is the top ligand from its variant
    seqs = sorted(c.sequence for c in constructs)
    assert seqs == sorted(["KLQGHSAPV", "NVDEILGRW", "VVVGADGVG"])
    for c in constructs:
        assert c.components['antigen_content'] == "minimal_epitope"
        # Antigen name disambiguator carries the _epitope suffix
        assert c.antigen_names[0].endswith("_epitope")


def test_peptide_concat_minimal_ligands():
    """Multi-epitope minimal: concatenated short MHC ligands. New
    combo unlocked by the orthogonal axes (was not expressible
    pre-PR)."""
    options = PeptideConstructConfig(
        antigen_content="minimal_epitope",
        antigens_per_construct=3,
        epitopes_per_antigen=1,
        linker="AAY",  # canonical CTL-spacer for concatenated minimal epitopes
        max_constructs=10,
        min_antigen_length_aa=5)
    constructs = assemble_peptide_constructs(_ranked_three_variants(), options)
    [c] = constructs
    assert len(c.antigen_names) == 3
    # Three 9-mers with AAY linkers: 9 + 3 + 9 + 3 + 9 = 33 aa
    assert len(c.sequence) == 33
    assert c.sequence.count("AAY") == 2  # 2 linkers between 3 antigens
    assert c.components['antigen_content'] == "minimal_epitope"


def test_peptide_minimal_with_multiple_top_ligands():
    """epitopes_per_antigen > 1: multiple top ligands from the same
    variant, each as its own antigen. Pins that we don't hardcode
    'just one top ligand matters'."""
    options = PeptideConstructConfig(
        antigen_content="minimal_epitope",
        antigens_per_construct=1,  # 1 construct per ligand
        epitopes_per_antigen=2,    # take top-2 from each VP
        min_antigen_length_aa=5,
        max_constructs=20)
    constructs = assemble_peptide_constructs(_ranked_three_variants(), options)
    # GENEA has 3 ligands → take 2; GENEB has 2 → take 2; GENEC has 1 → take 1
    # Total: 5 single-ligand constructs.
    assert len(constructs) == 5
    seqs = sorted(c.sequence for c in constructs)
    assert seqs == sorted([
        "KLQGHSAPV", "LQGHSAPVL",      # GENEA top-2
        "NVDEILGRW", "VDEILGRWE",      # GENEB top-2
        "VVVGADGVG",                    # GENEC only has 1
    ])
    # Names include the disambiguating suffix when >1 ligand from same VP
    suffixes = sorted(c.antigen_names[0].split("_")[-1] for c in constructs)
    assert "epitope1" in suffixes
    assert "epitope2" in suffixes


# ---------------------------------------------------------------------------
# mRNA design matrix (4 cells)
# ---------------------------------------------------------------------------

def test_mrna_biontech_canonical():
    """BioNTech FixVac / iNeST: mutation_spanning, N per construct,
    HLA-B signal peptide + MITD + tandem HBB UTRs. The default
    RNAConstructConfig — pinned here as a regression test."""
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,  # simplify slice math
        antigens_per_construct=3,
        max_antigen_length_aa=15,
        min_antigen_length_aa=5,
        utr_3p='HBB',  # not HBB_FI for slice simplicity
        poly_a_length=20,
        optimize_linkers=False,
        max_constructs=1)
    constructs = assemble_mrna_constructs(_ranked_three_variants(), options)
    assert len(constructs) == 1
    [c] = constructs
    assert len(c.antigens) == 3
    # Each antigen is a mutation-centered window (≤ 15 aa)
    for ant in c.antigens:
        assert ant['length_aa'] <= 15


def test_mrna_single_antigen():
    """Single-antigen mRNA: mutation_spanning, 1 per construct. One
    SLP per construct, multiple constructs spilling. Less common
    clinically but a valid design."""
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        antigens_per_construct=1,
        max_constructs=10,
        max_antigen_length_aa=15,
        min_antigen_length_aa=5,
        utr_3p='HBB',
        poly_a_length=20,
        optimize_linkers=False)
    constructs = assemble_mrna_constructs(_ranked_three_variants(), options)
    assert len(constructs) == 3
    for c in constructs:
        assert len(c.antigens) == 1


def test_mrna_string_of_beads():
    """String-of-beads mRNA: minimal_epitope, N per construct.
    Concatenated short MHC ligands in an mRNA backbone — analogous
    to Velten et al. 2018's polyepitope mRNA design."""
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        antigen_content="minimal_epitope",
        epitopes_per_antigen=1,
        antigens_per_construct=3,
        max_constructs=1,
        max_antigen_length_aa=15,
        min_antigen_length_aa=5,
        linker="AAY",  # CTL spacer for minimal-epitope mRNA
        utr_3p='HBB',
        poly_a_length=20,
        optimize_linkers=False)
    constructs = assemble_mrna_constructs(_ranked_three_variants(), options)
    [c] = constructs
    assert len(c.antigens) == 3
    # Each antigen is a 9-aa minimal ligand
    for ant in c.antigens:
        assert ant['length_aa'] == 9
        # The CDS contains the 9-aa ligand sequences
    # AA names carry the _epitope suffix from minimal_epitope mode
    for name in c.antigen_names:
        assert name.endswith("_epitope")


def test_mrna_single_ligand():
    """Single-ligand mRNA: minimal_epitope, 1 per construct. Unusual
    (mRNA payload usually wants more density) but a valid corner of
    the design matrix."""
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        antigen_content="minimal_epitope",
        epitopes_per_antigen=1,
        antigens_per_construct=1,
        max_constructs=10,
        max_antigen_length_aa=15,
        min_antigen_length_aa=5,
        utr_3p='HBB',
        poly_a_length=20,
        optimize_linkers=False)
    constructs = assemble_mrna_constructs(_ranked_three_variants(), options)
    assert len(constructs) == 3
    for c in constructs:
        assert len(c.antigens) == 1
        assert c.antigens[0]['length_aa'] == 9


def test_mrna_string_of_beads_multiple_top_ligands():
    """epitopes_per_antigen > 1 in mRNA: multiple top ligands from
    the same variant become separate antigens. Pins the
    "don't-assume-just-one-matters" semantics for the mRNA path too."""
    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        antigen_content="minimal_epitope",
        epitopes_per_antigen=3,  # all GENEA ligands; GENEB top-2; GENEC top-1
        antigens_per_construct=10,  # generous cap
        max_constructs=1,
        max_antigen_length_aa=15,
        min_antigen_length_aa=5,
        linker="AAY",
        utr_3p='HBB',
        poly_a_length=20,
        max_length_nt=4000,
        optimize_linkers=False)
    constructs = assemble_mrna_constructs(_ranked_three_variants(), options)
    [c] = constructs
    # 3 (GENEA) + 2 (GENEB) + 1 (GENEC) = 6 minimal-ligand antigens
    assert len(c.antigens) == 6
    # Each is a 9-aa ligand
    for ant in c.antigens:
        assert ant['length_aa'] == 9


# ---------------------------------------------------------------------------
# Validation: invalid axis values raise
# ---------------------------------------------------------------------------

def test_peptide_invalid_antigen_content_raises():
    with pytest.raises(ValueError, match="antigen_content"):
        PeptideConstructConfig(antigen_content="garbage")


def test_peptide_invalid_mode_raises():
    with pytest.raises(ValueError, match="mode must be one of"):
        PeptideConstructConfig(mode="not_a_mode")


def test_mrna_invalid_antigen_content_raises():
    """mRNA dispatches at antigen-extraction time; invalid content
    surfaces there with a clear message."""
    options = RNAConstructConfig(
        antigen_content="garbage", utr_3p='HBB', poly_a_length=0,
        optimize_linkers=False)
    with pytest.raises(ValueError, match="antigen_content"):
        assemble_mrna_constructs(_ranked_three_variants(), options)
