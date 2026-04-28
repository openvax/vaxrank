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

"""Tests for the shared vaccine_library and the mrna_library compatibility shim."""

import pytest

from vaxrank import mrna_library, vaccine_library


def test_mrna_library_exports_old_string_constants():
    # 2.10.0 publicly exposed plain-string linker constants. Ensure the
    # 2.11+ refactor (Linker dataclass in vaccine_library) preserves
    # those imports.
    assert mrna_library.LINKER_GS == "GGGGS"
    assert mrna_library.LINKER_GS3 == "GGGGSGGGGSGGGGS"
    assert mrna_library.LINKER_AAY == "AAY"
    assert mrna_library.LINKER_GPGPG == "GPGPG"


def test_mrna_library_linkers_dict_returns_strings():
    # Pre-refactor consumers used LINKERS["GS3"] expecting a string,
    # not a Linker object; preserve that shape in mrna_library.
    for name, value in mrna_library.LINKERS.items():
        assert isinstance(value, str), \
            "mrna_library.LINKERS['%s'] must remain a string for compat" % name


def test_vaccine_library_exposes_linker_objects():
    # The richer dataclass form lives in vaccine_library; the canonical
    # name for the (G4S)3 linker is now "G4S3", with "GS3" as a back-compat
    # alias resolved by get_linker.
    g4s3 = vaccine_library.LINKERS["G4S3"]
    assert isinstance(g4s3, vaccine_library.Linker)
    assert g4s3.amino_acids == "GGGGSGGGGSGGGGS"
    assert g4s3.dna is None
    assert g4s3.freeze_in_mrna is False
    # Aliases dereference to the canonical Linker object.
    assert vaccine_library.get_linker("GS3") is g4s3
    assert vaccine_library.get_linker("GS") is vaccine_library.LINKERS["G4S"]


def test_vaccine_library_p2a_carries_blessed_dna():
    p2a = vaccine_library.LINKERS["P2A"]
    assert p2a.amino_acids == "GSGATNFSLLKQAGDVEENPGP"
    # AA / DNA length consistency: 22 aa * 3 = 66 nt
    assert len(p2a.dna) == len(p2a.amino_acids) * 3
    assert p2a.freeze_in_mrna is True
    assert p2a.inert_in_peptide_mode is True


_STANDARD_CODON_TABLE = {
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
}


@pytest.mark.parametrize("name", ["P2A", "T2A", "F2A", "E2A"])
def test_2a_blessed_dna_translates_to_aa(name):
    """Each 2A entry's DNA must (a) be 3x its AA length and (b) translate
    back to the same amino-acid string. A typo in any blessed DNA would
    silently produce a non-functional 2A in mRNA constructs."""
    linker = vaccine_library.LINKERS[name]
    assert linker.dna is not None, "%s missing blessed DNA" % name
    assert len(linker.dna) == len(linker.amino_acids) * 3, \
        "%s DNA length mismatch" % name
    translated = ''.join(
        _STANDARD_CODON_TABLE[linker.dna[i:i + 3]]
        for i in range(0, len(linker.dna), 3)
    )
    assert translated == linker.amino_acids, \
        "%s DNA translates to '%s', expected '%s'" % (
            name, translated, linker.amino_acids)


@pytest.mark.parametrize("name", ["P2A", "T2A", "F2A", "E2A"])
def test_2a_linker_flags(name):
    linker = vaccine_library.LINKERS[name]
    assert linker.freeze_in_mrna is True
    assert linker.inert_in_peptide_mode is True


def test_get_linker_is_case_sensitive():
    # CLI normalizes case before reaching get_linker via type=str.upper, so
    # the underlying function expects the canonical (uppercase) form.
    with pytest.raises(ValueError):
        vaccine_library.get_linker("p2a")


def test_construct_name_format_consistent_across_modalities():
    """mRNA constructs use 'seq_NNN'; peptide constructs use 'peptide_NNN'.
    Same width and zero-padding so consumers can sort/parse uniformly.
    """
    import re
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.mrna import MRNAOptions, assemble_mrna_constructs
    from vaxrank.peptide import PeptideOptions, assemble_peptide_constructs

    fragment = SimpleNamespace(
        amino_acids="KLQGHSAPVLDVIVN", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=15)
    peptide = SimpleNamespace(
        mutant_protein_fragment=fragment, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    pairs = [(Variant('1', 1000, 'A', 'T'), [peptide])]

    [m] = assemble_mrna_constructs(pairs, options=MRNAOptions())
    [p] = assemble_peptide_constructs(pairs, options=PeptideOptions())

    # Both should match <prefix>_NNN with the same numeric width.
    pattern = re.compile(r"^(seq|peptide)_\d{3}$")
    assert pattern.match(m.name), \
        "mRNA construct name '%s' must match <prefix>_NNN" % m.name
    assert pattern.match(p.name), \
        "Peptide construct name '%s' must match <prefix>_NNN" % p.name


def test_vaccine_library_get_linker_unknown_raises():
    with pytest.raises(ValueError):
        vaccine_library.get_linker("notalinker")


def test_signal_peptides_include_tcr_additions():
    # CD8A and CD28 added from openvax/tcrsift's vetted set; confirm
    # they're discoverable through the existing SIGNAL_PEPTIDES dict.
    assert mrna_library.SIGNAL_PEPTIDES["CD8A"] == "MALPVTALLLPLALLLHAARP"
    assert mrna_library.SIGNAL_PEPTIDES["CD28"] == "MLRLLLALNLFPSIQVTG"


# ---- expanded linker family + aliases -------------------------------------

@pytest.mark.parametrize("name,aa", [
    ("G2S", "GGS"),
    ("G3S", "GGGS"),
    ("G4S", "GGGGS"),
    ("G5S", "GGGGGS"),
    ("G4S2", "GGGGSGGGGS"),
    ("G4S3", "GGGGSGGGGSGGGGS"),
    ("G4S4", "GGGGSGGGGSGGGGSGGGGS"),
])
def test_gs_linker_family(name, aa):
    linker = vaccine_library.LINKERS[name]
    assert linker.amino_acids == aa
    assert linker.dna is None  # not codon-frozen
    assert linker.freeze_in_mrna is False


@pytest.mark.parametrize("name,aa", [
    ("EAAAK", "EAAAK"),
    ("EAAAK2", "EAAAKEAAAK"),
    ("EAAAK3", "EAAAKEAAAKEAAAK"),
])
def test_eaaak_linker_family(name, aa):
    linker = vaccine_library.LINKERS[name]
    assert linker.amino_acids == aa
    assert "Arai" in linker.citation


@pytest.mark.parametrize("name,aa", [
    ("RKRR", "RKRR"),
    ("RVKR", "RVKR"),
    ("RKRKR", "RKRKR"),
])
def test_furin_linkers(name, aa):
    linker = vaccine_library.LINKERS[name]
    assert linker.amino_acids == aa
    assert "furin" in linker.citation.lower() or "Thomas" in linker.citation


def test_gs_aliases_resolve_to_canonical():
    # 2.10.0 used "GS"/"GS3"; preserve as aliases of the canonical names.
    assert vaccine_library.get_linker("GS").name == "G4S"
    assert vaccine_library.get_linker("GS3").name == "G4S3"


def test_all_linker_names_includes_aliases():
    names = vaccine_library.all_linker_names()
    assert "G4S3" in names
    assert "GS3" in names  # alias still listed
    assert "P2A" in names
    assert "EAAAK2" in names
    assert "RKRR" in names


# ---- codon species normalization ------------------------------------------

@pytest.mark.parametrize("user_input,expected", [
    ("human", "h_sapiens"),
    ("HUMAN", "h_sapiens"),
    ("Homo sapiens", "h_sapiens"),
    ("h. sapiens", "h_sapiens"),
    ("h_sapiens", "h_sapiens"),
    ("9606", "h_sapiens"),
    ("mouse", "m_musculus"),
    ("murine", "m_musculus"),
    ("Mus musculus", "m_musculus"),
    ("m. musculus", "m_musculus"),
    ("10090", "m_musculus"),
    ("yeast", "s_cerevisiae"),
    ("E. coli", "e_coli"),
])
def test_normalize_codon_species(user_input, expected):
    from vaxrank.mrna import normalize_codon_species
    assert normalize_codon_species(user_input) == expected


def test_normalize_codon_species_passes_through_unknown():
    # If we don't recognize the input, leave it for DnaChisel to validate.
    from vaxrank.mrna import normalize_codon_species
    assert normalize_codon_species("foo_bar") == "foo_bar"
