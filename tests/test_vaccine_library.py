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
    # The richer dataclass form lives in vaccine_library. Single-unit
    # linkers (G4S, AAY, EAAAK, ...) live in LINKERS; repeats are
    # parsed compositionally from (G4S)3 / G4Sx3 etc.
    g4s = vaccine_library.LINKERS["G4S"]
    assert isinstance(g4s, vaccine_library.Linker)
    assert g4s.amino_acids == "GGGGS"
    assert g4s.dna is None
    assert g4s.freeze_in_mrna is False
    # GS alias still resolves to the static G4S entry.
    assert vaccine_library.get_linker("GS") is g4s
    # GS3 alias forwards to the compositional form (G4S)3.
    g4s3 = vaccine_library.get_linker("GS3")
    assert g4s3.amino_acids == "GGGGSGGGGSGGGGS"


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


def test_get_linker_is_case_insensitive():
    # 2.13: get_linker uppercases the input itself so 'p2a', 'P2A',
    # 'g4sx2', 'G4Sx2' all resolve identically.
    assert vaccine_library.get_linker("p2a") is vaccine_library.LINKERS["P2A"]
    assert vaccine_library.get_linker("g4s") is vaccine_library.LINKERS["G4S"]
    assert vaccine_library.get_linker("g4sx2").amino_acids == "GGGGSGGGGS"


def test_construct_name_format_consistent_across_modalities():
    """mRNA constructs use 'seq_NNN'; peptide constructs use 'peptide_NNN'.
    Same width and zero-padding so consumers can sort/parse uniformly.
    """
    import re
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs

    fragment = SimpleNamespace(
        amino_acids="KLQGHSAPVLDVIVN", gene_name='G',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=15)
    peptide = SimpleNamespace(
        mutant_protein_fragment=fragment, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    pairs = [(Variant('1', 1000, 'A', 'T'), [peptide])]

    [m] = assemble_mrna_constructs(pairs, options=RNAConstructConfig())
    [p] = assemble_peptide_constructs(pairs, options=PeptideConstructConfig())

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
])
def test_gs_linker_family_static(name, aa):
    """Single-unit GnSm entries live in LINKERS. Repeats use the
    compositional grammar (see test_paren_repeat_grammar_g4s)."""
    linker = vaccine_library.LINKERS[name]
    assert linker.amino_acids == aa
    assert linker.dna is None
    assert linker.freeze_in_mrna is False


@pytest.mark.parametrize("name,aa", [
    ("(G4S)2", "GGGGSGGGGS"),
    ("(G4S)3", "GGGGSGGGGSGGGGS"),
    ("(G4S)4", "GGGGSGGGGSGGGGSGGGGS"),
    ("G4Sx2", "GGGGSGGGGS"),
    ("G4Sx3", "GGGGSGGGGSGGGGS"),
    ("(G4S)X3", "GGGGSGGGGSGGGGS"),  # caps X also accepted
])
def test_gs_linker_repeats_via_grammar(name, aa):
    linker = vaccine_library.get_linker(name)
    assert linker.amino_acids == aa


def test_gs_compact_form_now_literal():
    # 2.13 breaking change: bare 'G4S2' is no longer (G4S)2 — it parses
    # as the literal "GGGGSS" (4 glycines + 2 serines). This avoids
    # ambiguity with the (G4S)2 / G4Sx2 repeat forms.
    assert vaccine_library.get_linker("G4S2").amino_acids == "GGGGSS"
    assert vaccine_library.get_linker("G4S3").amino_acids == "GGGGSSS"


def test_eaaak_static_only():
    eaaak = vaccine_library.LINKERS["EAAAK"]
    assert eaaak.amino_acids == "EAAAK"
    assert "Arai" in eaaak.citation
    # Repeats use the compositional grammar.
    assert vaccine_library.get_linker("(EAAAK)3").amino_acids == "EAAAK" * 3
    assert vaccine_library.get_linker("EAAAKx2").amino_acids == "EAAAK" * 2


@pytest.mark.parametrize("name,aa", [
    ("RKRR", "RKRR"),
    ("RVKR", "RVKR"),
    ("RKRKR", "RKRKR"),
])
def test_furin_linkers(name, aa):
    linker = vaccine_library.LINKERS[name]
    assert linker.amino_acids == aa
    assert "furin" in linker.citation.lower() or "Thomas" in linker.citation


def test_gs_aliases_resolve():
    # 'GS' → static G4S; 'GS3' → compositional (G4S)3. Both 2.10.0 names
    # still resolve to a Linker with the original 2.10.0 amino-acid string.
    assert vaccine_library.get_linker("GS").amino_acids == "GGGGS"
    assert vaccine_library.get_linker("GS3").amino_acids == "GGGGSGGGGS" + "GGGGS"


def test_all_linker_names_lists_static_and_aliases():
    names = vaccine_library.all_linker_names()
    assert "G4S" in names
    assert "GS3" in names  # alias still listed
    assert "P2A" in names
    assert "RKRR" in names
    # G4S3, G4S2, G4S4, EAAAK2, EAAAK3 are NOT static — they're parsed
    # compositionally and don't appear in the listing.
    assert "G4S3" not in names
    assert "EAAAK2" not in names


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


def test_select_antigen_window_centers_on_mutation():
    # Regression for the bug where SLP truncation took [:max_length]
    # and could drop the mutation when it sat past the head of the
    # fragment. The shared helper must center on mutant_amino_acid_*_offset.
    from types import SimpleNamespace
    from vaxrank.vaccine_library import select_antigen_window
    fragment = SimpleNamespace(
        amino_acids="A" * 40 + "WHY" + "A" * 7,  # mutation at 40-43
        mutant_amino_acid_start_offset=40,
        mutant_amino_acid_end_offset=43,
    )
    window = select_antigen_window(fragment, "test", 30)
    assert len(window) == 30
    assert "WHY" in window


def test_select_antigen_window_short_fragment_passthrough():
    from types import SimpleNamespace
    from vaxrank.vaccine_library import select_antigen_window
    fragment = SimpleNamespace(
        amino_acids="KLQGHSAPVL",
        mutant_amino_acid_start_offset=0,
        mutant_amino_acid_end_offset=10,
    )
    assert select_antigen_window(fragment, "test", 30) == "KLQGHSAPVL"


def test_select_antigen_window_warns_when_mutation_exceeds_cap(caplog):
    import logging
    from types import SimpleNamespace
    from vaxrank.vaccine_library import select_antigen_window
    fragment = SimpleNamespace(
        amino_acids="A" * 5 + "M" * 25 + "A" * 10,  # 25 aa mutation
        mutant_amino_acid_start_offset=5,
        mutant_amino_acid_end_offset=30,
    )
    with caplog.at_level(logging.WARNING):
        result = select_antigen_window(fragment, "test", 20)
    # Untruncated since mutation > cap
    assert result == fragment.amino_acids
    assert any("longer than" in r.message for r in caplog.records)


# ---- end-to-end: both modalities from the same ranked list ----------------

def test_both_modalities_emit_compatible_manifests(tmp_path):
    """Smoke test: a single ranked list produces both peptide and mRNA
    constructs with manifest schemas that share the top-level keys
    (modality, name, length, length_unit, antigen_names, components,
    manufacturability). Catches divergence between the two writers.
    """
    import json
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.peptide import (
        PeptideConstructConfig, assemble_peptide_constructs,
        write_peptide_outputs)
    from vaxrank.mrna import (
        RNAConstructConfig, assemble_mrna_constructs, write_mrna_outputs)

    fragment = SimpleNamespace(
        amino_acids="KLQGHSAPVLDVIVN", gene_name='GENE',
        mutant_amino_acid_start_offset=5, mutant_amino_acid_end_offset=10)
    peptide = SimpleNamespace(
        mutant_protein_fragment=fragment, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    pairs = [(Variant('1', 1000, 'A', 'T'), [peptide])]

    p_constructs = assemble_peptide_constructs(
        pairs, options=PeptideConstructConfig())
    m_constructs = assemble_mrna_constructs(
        pairs, options=RNAConstructConfig())

    p_fasta = tmp_path / "peptide.fasta"
    p_manifest = tmp_path / "peptide.json"
    m_fasta = tmp_path / "mrna.fasta"
    m_manifest = tmp_path / "mrna.json"
    write_peptide_outputs(
        p_constructs, str(p_fasta), str(p_manifest))
    write_mrna_outputs(
        m_constructs, str(m_fasta), str(m_manifest))

    p_entry = json.loads(p_manifest.read_text())[0]
    m_entry = json.loads(m_manifest.read_text())[0]

    shared_keys = {
        'modality', 'name', 'length', 'length_unit',
        'antigen_names', 'components', 'manufacturability',
    }
    assert shared_keys <= set(p_entry), \
        "peptide manifest missing keys: %s" % (shared_keys - set(p_entry),)
    assert shared_keys <= set(m_entry), \
        "mRNA manifest missing keys: %s" % (shared_keys - set(m_entry),)
    assert p_entry['modality'] == 'peptide'
    assert p_entry['length_unit'] == 'aa'
    assert m_entry['modality'] == 'mrna'
    assert m_entry['length_unit'] == 'nt'


# ---- iter_named_antigens shared helper (#248) -----------------------------

def test_iter_named_antigens_naming_format():
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.vaccine_library import iter_named_antigens

    fragment = SimpleNamespace(
        amino_acids="KLQGH", gene_name='GENE',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=5)
    peptide = SimpleNamespace(
        mutant_protein_fragment=fragment, mutant_epitope_predictions=[])
    pairs = [(Variant('1', 1000, 'A', 'T'), [peptide])]
    [(name, frag, pep)] = list(iter_named_antigens(pairs))
    assert name == "GENE_1_1000_A_T"
    assert frag is fragment
    assert pep is peptide


def test_iter_named_antigens_alt_suffix():
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.vaccine_library import iter_named_antigens

    def _peptide(aa):
        f = SimpleNamespace(
            amino_acids=aa, gene_name='G',
            mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=len(aa))
        return SimpleNamespace(mutant_protein_fragment=f)

    pairs = [(Variant('1', 100, 'A', 'T'), [_peptide("ABC"), _peptide("DEF"), _peptide("GHI")])]
    names = [n for n, _, _ in iter_named_antigens(pairs, candidates_per_slot=3)]
    assert names == ["G_1_100_A_T", "G_1_100_A_T_alt1", "G_1_100_A_T_alt2"]


def test_iter_named_antigens_skips_empty_peptide_lists():
    from varcode import Variant
    from vaxrank.vaccine_library import iter_named_antigens
    pairs = [(Variant('1', 100, 'A', 'T'), [])]
    assert list(iter_named_antigens(pairs)) == []


def test_iter_named_antigens_handles_missing_gene_name():
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.vaccine_library import iter_named_antigens
    fragment = SimpleNamespace(
        amino_acids="KLQ", gene_name=None,
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=3)
    peptide = SimpleNamespace(mutant_protein_fragment=fragment)
    pairs = [(Variant('1', 100, 'A', 'T'), [peptide])]
    [(name, _, _)] = list(iter_named_antigens(pairs))
    assert name == "unknown_1_100_A_T"


def test_iter_named_antigens_caps_at_candidates_per_slot():
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.vaccine_library import iter_named_antigens

    def _peptide(aa):
        f = SimpleNamespace(
            amino_acids=aa, gene_name='G',
            mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=len(aa))
        return SimpleNamespace(mutant_protein_fragment=f)

    pairs = [(Variant('1', 100, 'A', 'T'),
              [_peptide("ABC"), _peptide("DEF"), _peptide("GHI")])]
    # candidates_per_slot=2 → only top 2 walked
    out = list(iter_named_antigens(pairs, candidates_per_slot=2))
    assert len(out) == 2


def test_peptide_and_mrna_names_match_for_same_input():
    """mRNA and peptide modes use the shared iter_named_antigens helper,
    so the antigen names in their respective manifests are identical
    for any given ranked input. Catches a regression where one modality
    forks its own naming convention.
    """
    from types import SimpleNamespace
    from varcode import Variant
    from vaxrank.peptide import (
        PeptideConstructConfig, assemble_peptide_constructs)
    from vaxrank.mrna import (
        RNAConstructConfig, assemble_mrna_constructs)

    fragment = SimpleNamespace(
        amino_acids="KLQGHSAPVLDVIVN", gene_name='GENE',
        mutant_amino_acid_start_offset=5, mutant_amino_acid_end_offset=10)
    peptide = SimpleNamespace(
        mutant_protein_fragment=fragment, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    pairs = [(Variant('1', 1000, 'A', 'T'), [peptide])]

    [p] = assemble_peptide_constructs(pairs, options=PeptideConstructConfig())
    [m] = assemble_mrna_constructs(pairs, options=RNAConstructConfig())
    assert p.antigen_names == m.antigen_names


# ---- compositional grammar (#247 prep) ------------------------------------

def test_paren_repeat_grammar_g2s():
    linker = vaccine_library.get_linker("(G2S)3")
    assert linker.amino_acids == "GGSGGSGGS"
    assert linker.name == "(G2S)3"
    # Inherits citation from base (G2S → GS family)
    assert "Huston" in linker.citation


def test_paren_repeat_grammar_aay():
    linker = vaccine_library.get_linker("(AAY)2")
    assert linker.amino_acids == "AAYAAY"
    assert "Velders" in linker.citation


def test_paren_repeat_rejects_2a():
    # 2A skipping is positional; repeating it doesn't add cleavage events.
    with pytest.raises(ValueError, match="codon-frozen"):
        vaccine_library.get_linker("(P2A)2")


def test_paren_repeat_caps_count():
    with pytest.raises(ValueError, match="repeat count"):
        vaccine_library.get_linker("(G4S)1000")


def test_gnsm_grammar_literal():
    # GnSm parses literally: n glycines + m serines, single unit.
    # Repeats require explicit parens or the 'x' separator.
    assert vaccine_library.get_linker("G6S").amino_acids == "GGGGGGS"
    assert vaccine_library.get_linker("G2S5").amino_acids == "GG" + "SSSSS"
    assert vaccine_library.get_linker("G4S3").amino_acids == "GGGG" + "SSS"


def test_any_grammar():
    # A2Y = canonical AAY (resolves to static entry first via its
    # canonical name rather than the regex; the regex fires for n != 2)
    a3 = vaccine_library.get_linker("A3Y")
    assert a3.amino_acids == "AAAY"
    assert a3.name == "A3Y"
    a4 = vaccine_library.get_linker("A4Y")
    assert a4.amino_acids == "AAAAY"
    a8 = vaccine_library.get_linker("A8Y")
    assert a8.amino_acids == "A" * 8 + "Y"


def test_any_citation_flags_extrapolation():
    # A3Y onwards: citation must explicitly note the extrapolation.
    a3 = vaccine_library.get_linker("A3Y")
    assert ("extrapolation" in a3.citation.lower()
            or "less standardized" in a3.citation.lower()
            or "without independent" in a3.citation.lower())


def test_paren_and_compact_forms_now_differ_for_g_family():
    # 2.13 breaking change: bare 'G2S3' is the LITERAL "GGSSS" (2 G's
    # + 3 S's), not the (G2S)3 repeat. The two forms now mean different
    # things by design — see module docstring.
    paren = vaccine_library.get_linker("(G2S)3").amino_acids
    compact = vaccine_library.get_linker("G2S3").amino_acids
    assert paren == "GGSGGSGGS"
    assert compact == "GGSSS"
    assert paren != compact


def test_unknown_compositional_form_raises():
    # Random gibberish doesn't match any of the regexes.
    with pytest.raises(ValueError, match="Unknown linker"):
        vaccine_library.get_linker("ZZZ123")


def test_grammar_rejects_zero_count():
    with pytest.raises(ValueError):
        vaccine_library.get_linker("(G4S)0")
    with pytest.raises(ValueError):
        vaccine_library.get_linker("A0Y")
