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
    MRNAOptions,
    assemble_mrna_constructs,
    codon_optimize,
    write_mrna_outputs,
)
from vaxrank.mrna_library import (
    MITD_HLA_A,
    SIGNAL_PEPTIDE_TPA,
    UTR_3P_HBB,
    UTR_5P_HBB,
)
from vaxrank.vaccine_library import LINKER_GS3


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
    options = MRNAOptions(signal_peptide='tPA', linker='GS3',
                          include_mitd=False)
    constructs = assemble_mrna_constructs(pairs, options=options)
    assert len(constructs) == 1
    c = constructs[0]
    assert c.antigen_names == ['GENEA_1_100_A_T', 'GENEB_2_200_A_T']
    assert c.sequence.startswith(UTR_5P_HBB)
    assert c.sequence.endswith(UTR_3P_HBB)
    # CDS sits between UTRs and ends with a stop codon and starts with ATG
    cds = c.sequence[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    assert cds[:3] == 'ATG'
    assert cds[-3:] == 'TAA'
    protein = _translate(cds[:-3])
    expected = SIGNAL_PEPTIDE_TPA + LINKER_GS3.amino_acids.join(['KLQGHSAPVLDVIVN', 'MNNVDEILGRWESPV'])
    assert protein == expected


def test_assembly_prepends_atg_when_no_signal_peptide():
    # Antigens that don't naturally start with M should still produce a
    # translatable CDS — assembler prepends an M.
    pairs = [_variant_pair("KLQGHSAPVLDVIVN", gene_name='GENE')]
    options = MRNAOptions(signal_peptide=None, linker='GS3', include_mitd=False)
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.sequence[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    assert cds[:3] == 'ATG'
    protein = _translate(cds[:-3])
    assert protein == "MKLQGHSAPVLDVIVN"


def test_assembly_does_not_double_prepend_when_antigen_starts_with_m():
    pairs = [_variant_pair("MNNVDEILGRWESPV", gene_name='GENE')]
    options = MRNAOptions(signal_peptide=None, linker='GS3', include_mitd=False)
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.sequence[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    protein = _translate(cds[:-3])
    assert protein == "MNNVDEILGRWESPV"


def test_assembly_with_mitd_appends_trafficking_domain():
    pairs = [_variant_pair("KLQGHSAPVLDVIVN")]
    options = MRNAOptions(signal_peptide=None, linker='GS3',
                          include_mitd=True, mitd='HLA_A')
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.sequence[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
    protein = _translate(cds[:-3])
    # No signal peptide selected and antigen doesn't start with M, so the
    # assembler prepends one to keep the CDS translatable.
    expected = "M" + "KLQGHSAPVLDVIVN" + LINKER_GS3.amino_acids + MITD_HLA_A
    assert protein == expected


def test_assembly_splits_on_max_antigens():
    pairs = [_variant_pair("AAAA", start=i) for i in range(5)]
    options = MRNAOptions(signal_peptide=None, linker='GS',
                          include_mitd=False, max_antigens_per_construct=2)
    constructs = assemble_mrna_constructs(pairs, options=options)
    assert [len(c.antigen_names) for c in constructs] == [2, 2, 1]


def test_assembly_splits_on_max_length():
    long_aa = "A" * 200
    pairs = [_variant_pair(long_aa, start=i) for i in range(4)]
    options = MRNAOptions(signal_peptide=None, linker='GS',
                          include_mitd=False, max_antigens_per_construct=10,
                          max_length_nt=len(UTR_5P_HBB) + len(UTR_3P_HBB) + 3 + 600 * 2)
    constructs = assemble_mrna_constructs(pairs, options=options)
    assert len(constructs) > 1


def test_assembly_no_antigens_returns_empty():
    assert assemble_mrna_constructs([], options=MRNAOptions()) == []


def test_unknown_signal_peptide_raises():
    pairs = [_variant_pair("KLQ")]
    with pytest.raises(ValueError):
        assemble_mrna_constructs(pairs, options=MRNAOptions(signal_peptide='nope'))


def test_assembly_p2a_linker_preserves_blessed_dna():
    # When a 2A linker is selected, the blessed DNA from the literature
    # must appear verbatim in the optimized CDS.
    from vaxrank.vaccine_library import LINKER_P2A
    pairs = [
        _variant_pair("KLQGHSAPVLDVIVN", contig='1', start=100, gene_name='GENEA'),
        _variant_pair("MNNVDEILGRWESPV", contig='2', start=200, gene_name='GENEB'),
    ]
    options = MRNAOptions(signal_peptide='tPA', linker='P2A',
                          include_mitd=False)
    [c] = assemble_mrna_constructs(pairs, options=options)
    cds = c.sequence[len(UTR_5P_HBB):-len(UTR_3P_HBB)]
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
        options=MRNAOptions(signal_peptide=None, include_mitd=False))
    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "out.fasta")
        manifest_path = os.path.join(tmp, "out.json")
        write_mrna_outputs(constructs, fasta_path, manifest_path)
        with open(fasta_path) as f:
            text = f.read()
        assert text.startswith(">seq_001")
        assert UTR_5P_HBB in text.replace('\n', '')
        with open(manifest_path) as f:
            manifest = json.load(f)
        entry = manifest[0]
        assert entry['modality'] == 'mrna'
        assert entry['length_unit'] == 'nt'
        assert entry['name'] == 'seq_001'
        assert entry['antigen_names'] == ['GENE_1_1000_A_T']
        assert entry['length'] == len(constructs[0].sequence)
        assert 'components' in entry
        assert 'manufacturability' in entry
