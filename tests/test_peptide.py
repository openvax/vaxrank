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

"""Tests for peptide vaccine construct assembly."""

import csv
import json
import os
import tempfile
from types import SimpleNamespace

import pytest
from varcode import Variant

from vaxrank.peptide import (
    PeptideConstructConfig,
    assemble_peptide_constructs,
    write_peptide_outputs,
)


def _peptide_stub(amino_acids, gene_name='GENE',
                  mutant_epitope_predictions=None,
                  manufacturability_scores=None,
                  mut_start=None, mut_end=None):
    fragment = SimpleNamespace(
        amino_acids=amino_acids,
        gene_name=gene_name,
        mutant_amino_acid_start_offset=(0 if mut_start is None else mut_start),
        mutant_amino_acid_end_offset=(
            len(amino_acids) if mut_end is None else mut_end),
    )
    return SimpleNamespace(
        mutant_protein_fragment=fragment,
        mutant_epitope_predictions=mutant_epitope_predictions or [],
        manufacturability_scores=manufacturability_scores,
    )


def _variant_pair(amino_acids, contig='1', start=1000, gene_name='GENE',
                  epitopes=None, manufacturability=None,
                  mut_start=None, mut_end=None):
    variant = Variant(contig, start, 'A', 'T')
    peptide = _peptide_stub(
        amino_acids, gene_name=gene_name,
        mutant_epitope_predictions=epitopes,
        manufacturability_scores=manufacturability,
        mut_start=mut_start, mut_end=mut_end)
    return (variant, [peptide])


def test_slp_mode_one_construct_per_peptide():
    pairs = [
        _variant_pair("KLQGHSAPVLDVIVN", start=100, gene_name='GENEA'),
        _variant_pair("MNNVDEILGRWESPV", start=200, gene_name='GENEB'),
    ]
    constructs = assemble_peptide_constructs(pairs, options=PeptideConstructConfig())
    assert len(constructs) == 2
    assert constructs[0].sequence == "KLQGHSAPVLDVIVN"
    assert constructs[0].name == "peptide_001"
    assert constructs[0].antigen_names == ["GENEA_1_100_A_T"]
    assert constructs[1].sequence == "MNNVDEILGRWESPV"
    assert constructs[1].name == "peptide_002"
    assert constructs[0].components['mode'] == 'slp'


def test_slp_truncates_oversize_peptides():
    long_aa = "A" * 50
    [c] = assemble_peptide_constructs(
        [_variant_pair(long_aa, mut_start=24, mut_end=26)],
        options=PeptideConstructConfig(max_antigen_length_aa=25))
    assert c.sequence == "A" * 25


def test_slp_truncation_preserves_mutation():
    # 50-aa fragment with the mutation at offset 40:43; a naive [:30]
    # truncation would drop the mutation entirely. The assembler must
    # center the window on the mutation instead.
    aa = "A" * 40 + "WHY" + "A" * 7  # mutation residues at 40-43
    [c] = assemble_peptide_constructs(
        [_variant_pair(aa, mut_start=40, mut_end=43)],
        options=PeptideConstructConfig(max_antigen_length_aa=30))
    assert len(c.sequence) == 30
    assert "WHY" in c.sequence


def test_slp_emits_full_fragment_when_mutation_exceeds_cap():
    # Inframe insertion longer than --peptide-max-length-aa: refuse
    # to truncate (the mutation can't be retained otherwise).
    aa = "A" * 5 + "M" * 25 + "A" * 10  # 40 aa, mutation spans 25 aa
    [c] = assemble_peptide_constructs(
        [_variant_pair(aa, mut_start=5, mut_end=30)],
        options=PeptideConstructConfig(max_antigen_length_aa=20))
    assert c.sequence == aa


def test_minimal_epitope_uses_top_prediction():
    epitope = SimpleNamespace(peptide_sequence="KLAGHSPVL")
    pairs = [_variant_pair(
        "KLQGHSAPVLDVIVNCDESLLAS", gene_name='GENEA',
        epitopes=[epitope])]
    [c] = assemble_peptide_constructs(
        pairs, options=PeptideConstructConfig(mode='minimal_epitope'))
    assert c.sequence == "KLAGHSPVL"
    assert c.antigen_names == ["GENEA_1_1000_A_T_epitope"]


def test_minimal_epitope_skips_peptides_without_predictions():
    pairs = [
        _variant_pair("KLQGHSAPVLDVIVN", gene_name='GENEA', epitopes=[]),
        _variant_pair(
            "MNNVDEILGRWESPV", start=200, gene_name='GENEB',
            epitopes=[SimpleNamespace(peptide_sequence="MNNVDEILG")]),
    ]
    constructs = assemble_peptide_constructs(
        pairs, options=PeptideConstructConfig(mode='minimal_epitope'))
    assert len(constructs) == 1
    assert constructs[0].sequence == "MNNVDEILG"


def test_multi_epitope_concatenates_with_linker():
    pairs = [
        _variant_pair("KLQGH", start=100, gene_name='GENEA'),
        _variant_pair("MNNVD", start=200, gene_name='GENEB'),
    ]
    [c] = assemble_peptide_constructs(
        pairs,
        options=PeptideConstructConfig(mode='multi_epitope', linker='AAY',
                               max_antigen_length_aa=100,
                               antigens_per_construct=2))
    assert c.sequence == "KLQGHAAYMNNVD"
    assert c.components['linker'] == 'AAY'
    assert c.components['linker_inert'] is False


def test_multi_epitope_2a_linker_marked_inert():
    pairs = [
        _variant_pair("KLQGH", start=100, gene_name='GENEA'),
        _variant_pair("MNNVD", start=200, gene_name='GENEB'),
    ]
    [c] = assemble_peptide_constructs(
        pairs,
        options=PeptideConstructConfig(mode='multi_epitope', linker='P2A',
                               max_antigen_length_aa=200,
                               antigens_per_construct=2))
    assert "GSGATNFSLLKQAGDVEENPGP" in c.sequence
    assert c.components['linker_inert'] is True


def test_multi_epitope_splits_on_max_length():
    # 4 antigens at 20 aa each. With antigens_per_construct=2 and
    # max_antigen_length_aa=20, each construct caps at ~45 aa (2*20 + 1*5).
    # We'd need 2 constructs to fit all 4.
    pairs = [_variant_pair("A" * 20, start=i) for i in range(4)]
    constructs = assemble_peptide_constructs(
        pairs,
        options=PeptideConstructConfig(mode='multi_epitope', linker='GS',
                               max_antigen_length_aa=20,
                               antigens_per_construct=2,
                               max_constructs=10))
    assert len(constructs) >= 2


def test_slp_manufacturability_recomputed_after_truncation():
    # Source fragment ends in C (would have cysteine_count=1); after
    # truncation to a window that excludes the C, the construct's
    # manufacturability must reflect the emitted sequence, not the
    # source.
    aa = "A" * 30 + "WHY" + "A" * 6 + "C"  # 40-aa, mutation at 30-33
    [c] = assemble_peptide_constructs(
        [_variant_pair(aa, mut_start=30, mut_end=33)],
        options=PeptideConstructConfig(max_antigen_length_aa=15))
    assert "C" not in c.sequence
    assert c.manufacturability['cysteine_count'] == 0


def test_minimal_epitope_manufacturability_matches_emitted_sequence():
    # The emitted sequence is the predicted epitope, not the full
    # source vaccine peptide; manufacturability must follow the emitted
    # sequence.
    epitope = SimpleNamespace(peptide_sequence="KAAAAAA")  # no cysteines
    pairs = [_variant_pair(
        "KAAAAAACCC",  # source has cysteines; emitted does not
        epitopes=[epitope])]
    [c] = assemble_peptide_constructs(
        pairs, options=PeptideConstructConfig(mode='minimal_epitope'))
    assert c.manufacturability['cysteine_count'] == 0


def test_multi_epitope_manufacturability_populated():
    pairs = [
        _variant_pair("KLQGH", start=100, gene_name='GENEA'),
        _variant_pair("MNNCD", start=200, gene_name='GENEB'),
    ]
    [c] = assemble_peptide_constructs(
        pairs,
        options=PeptideConstructConfig(mode='multi_epitope', linker='AAY',
                               max_antigen_length_aa=100,
                               antigens_per_construct=2))
    assert c.manufacturability  # not empty
    assert c.manufacturability['cysteine_count'] == 1


def test_unknown_mode_raises():
    with pytest.raises(ValueError):
        assemble_peptide_constructs(
            [_variant_pair("KLQ")],
            options=PeptideConstructConfig(mode='nonsense'))


def test_no_antigens_returns_empty():
    assert assemble_peptide_constructs([], options=PeptideConstructConfig()) == []


def test_write_peptide_outputs_fasta_manifest_orderform():
    pairs = [_variant_pair(
        "KLQGHSAPVLDVIVN",
        manufacturability=SimpleNamespace(
            cterm_7mer_gravy_score=0.1, max_7mer_gravy_score=0.2,
            difficult_n_terminal_residue=False, c_terminal_cysteine=False,
            c_terminal_proline=False, cysteine_count=0,
            n_terminal_asparagine=False, n_terminal_methionine=False,
            aspartate_proline_bond_count=0,
        ))]
    options = PeptideConstructConfig(n_terminal_acetylation=True,
                             c_terminal_amidation=True)
    constructs = assemble_peptide_constructs(pairs, options=options)
    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "out.fasta")
        manifest_path = os.path.join(tmp, "out.json")
        order_path = os.path.join(tmp, "order.csv")
        write_peptide_outputs(
            constructs,
            fasta_path=fasta_path,
            manifest_path=manifest_path,
            order_form_path=order_path,
            options=options,
        )
        with open(fasta_path) as f:
            text = f.read()
        assert ">peptide_001" in text
        assert "KLQGHSAPVLDVIVN" in text

        with open(manifest_path) as f:
            manifest = json.load(f)
        entry = manifest[0]
        assert entry['modality'] == 'peptide'
        assert entry['length_unit'] == 'aa'
        assert entry['length'] == len("KLQGHSAPVLDVIVN")
        assert 'cterm_7mer_gravy_score' in entry['manufacturability']

        with open(order_path) as f:
            rows = list(csv.DictReader(f))
        row = rows[0]
        assert row['name'] == 'peptide_001'
        assert row['n_terminal_modification'] == 'Acetyl'
        assert row['c_terminal_modification'] == 'Amide'
        assert row['displayed_sequence'] == 'Ac-KLQGHSAPVLDVIVN-NH2'


def test_order_form_omits_displayed_sequence_when_no_modifications():
    """Without N-/C-terminal modifications, the displayed_sequence column
    would be a verbatim duplicate of `sequence` — so it should be dropped
    from the order form."""
    pairs = [_variant_pair("KLQGHSAPVLDVIVN")]
    options = PeptideConstructConfig()  # no mods
    constructs = assemble_peptide_constructs(pairs, options=options)
    with tempfile.TemporaryDirectory() as tmp:
        order_path = os.path.join(tmp, "order.csv")
        write_peptide_outputs(
            constructs,
            fasta_path=os.path.join(tmp, "out.fasta"),
            order_form_path=order_path,
            options=options,
        )
        with open(order_path) as f:
            rows = list(csv.DictReader(f))
        assert 'displayed_sequence' not in rows[0]
        assert rows[0]['n_terminal_modification'] == 'Free'
        assert rows[0]['c_terminal_modification'] == 'Free'


def test_candidates_per_slot_emits_alternates():
    # When candidates_per_slot=2, each variant with multiple ranked
    # peptides emits an "_alt1" alternate construct alongside the top.
    fragment_a = SimpleNamespace(
        amino_acids="KLQGHSAPVLD", gene_name='GENEA',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=11)
    fragment_b = SimpleNamespace(
        amino_acids="LQGHSAPVLDV", gene_name='GENEA',  # alternate window
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=11)
    peptide_a = SimpleNamespace(
        mutant_protein_fragment=fragment_a, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    peptide_b = SimpleNamespace(
        mutant_protein_fragment=fragment_b, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    from varcode import Variant
    pairs = [(Variant('1', 100, 'A', 'T'), [peptide_a, peptide_b])]
    constructs = assemble_peptide_constructs(
        pairs,
        options=PeptideConstructConfig(candidates_per_slot=2))
    assert len(constructs) == 2
    assert constructs[0].sequence == "KLQGHSAPVLD"
    assert constructs[1].sequence == "LQGHSAPVLDV"
    assert "_alt1" in constructs[1].antigen_names[0]


def test_candidates_per_slot_default_is_one():
    fragment_a = SimpleNamespace(
        amino_acids="KLQGHSAPVLD", gene_name='GENE',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=11)
    fragment_b = SimpleNamespace(
        amino_acids="LQGHSAPVLDV", gene_name='GENE',
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=11)
    pa = SimpleNamespace(
        mutant_protein_fragment=fragment_a, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    pb = SimpleNamespace(
        mutant_protein_fragment=fragment_b, mutant_epitope_predictions=[],
        manufacturability_scores=None)
    from varcode import Variant
    pairs = [(Variant('1', 100, 'A', 'T'), [pa, pb])]
    [c] = assemble_peptide_constructs(pairs, options=PeptideConstructConfig())
    assert c.sequence == "KLQGHSAPVLD"


def test_max_constructs_caps_peptide_pool():
    # Default max_constructs=20: anything beyond drops with a warning.
    pairs = []
    from varcode import Variant
    for i in range(25):
        fragment = SimpleNamespace(
            amino_acids="KLQGHSAPVLD", gene_name='G%d' % i,
            mutant_amino_acid_start_offset=0,
            mutant_amino_acid_end_offset=11)
        peptide = SimpleNamespace(
            mutant_protein_fragment=fragment, mutant_epitope_predictions=[],
            manufacturability_scores=None)
        pairs.append((Variant('1', 100 + i, 'A', 'T'), [peptide]))
    constructs = assemble_peptide_constructs(pairs, options=PeptideConstructConfig())
    assert len(constructs) == 20  # default


# ---- review-gap regression tests ------------------------------------------

def test_min_antigen_length_warning_when_emitted_below_floor(caplog):
    # Source fragment is 10 aa; min_antigen_length_aa=15 means the
    # emitted SLP will be too short. Vaxrank must warn rather than
    # silently shipping an undersized peptide.
    import logging
    pairs = [_variant_pair("KLQGHSAPVL", gene_name='G')]  # 10 aa
    with caplog.at_level(logging.WARNING):
        constructs = assemble_peptide_constructs(
            pairs, options=PeptideConstructConfig(min_antigen_length_aa=15))
    assert len(constructs) == 1
    assert any("below" in r.message and "min-antigen-length" in r.message
               for r in caplog.records), \
        "Expected a min-antigen-length warning for an undersized SLP"


def test_multi_epitope_warns_on_dropped_antigens(caplog):
    # 4 antigens, antigens_per_construct=2, max_constructs=1 → only
    # the first 2 antigens fit; the last 2 must produce a drop warning
    # so the user can see they're being silently dropped.
    import logging
    pairs = [_variant_pair("A" * 20, start=i) for i in range(4)]
    with caplog.at_level(logging.WARNING):
        constructs = assemble_peptide_constructs(
            pairs,
            options=PeptideConstructConfig(
                mode='multi_epitope', linker='GS', antigens_per_construct=2,
                max_constructs=1, max_antigen_length_aa=20))
    assert len(constructs) == 1
    assert len(constructs[0].antigen_names) == 2
    drop_warnings = [
        r for r in caplog.records
        if "max-constructs" in r.message
        and ("dropping" in r.message.lower() or "dropped" in r.message.lower())
    ]
    assert drop_warnings, \
        "Expected a max-constructs drop warning. Got: %s" % (
            [r.message for r in caplog.records],)


def test_order_form_includes_scale_purity_counterion():
    pairs = [_variant_pair("KLQGHSAPVLDVIVN")]
    options = PeptideConstructConfig(
        scale_mg=10.0, purity_percent=98.0, counterion='acetate')
    constructs = assemble_peptide_constructs(pairs, options=options)
    with tempfile.TemporaryDirectory() as tmp:
        order_path = os.path.join(tmp, "order.csv")
        write_peptide_outputs(
            constructs, fasta_path=os.path.join(tmp, "out.fasta"),
            order_form_path=order_path, options=options)
        with open(order_path) as f:
            rows = list(csv.DictReader(f))
        row = rows[0]
        assert float(row['scale_mg']) == 10.0
        assert float(row['purity_percent']) == 98.0
        assert row['counterion'] == 'acetate'


def test_order_form_default_scale_purity_counterion():
    pairs = [_variant_pair("KLQGHSAPVLDVIVN")]
    constructs = assemble_peptide_constructs(
        pairs, options=PeptideConstructConfig())
    with tempfile.TemporaryDirectory() as tmp:
        order_path = os.path.join(tmp, "order.csv")
        write_peptide_outputs(
            constructs, fasta_path=os.path.join(tmp, "out.fasta"),
            order_form_path=order_path)
        with open(order_path) as f:
            rows = list(csv.DictReader(f))
        row = rows[0]
        assert float(row['scale_mg']) == 5.0
        assert float(row['purity_percent']) == 95.0
        assert row['counterion'] == 'TFA'
