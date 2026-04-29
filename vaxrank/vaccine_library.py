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

"""Shared assembly utilities used by both mRNA and peptide construct modes.

Contains the linker library and the mutation-centered window selector
that both ``vaxrank.mrna`` and ``vaxrank.peptide`` need.

A linker is anything that goes between concatenated antigens. Some
linkers are pure amino-acid spacers (GS-family, AAY, GPGPG, EAAAK,
furin sites); 2A self-cleaving peptides also live here because the
*vocabulary* is shared even though the modality semantics differ:

- In mRNA constructs, 2A linkers cause ribosomal skipping; their
  blessed DNA is preserved during codon optimization
  (``freeze_in_mrna=True``).
- In synthesized peptides, 2A motifs do not cleave (the mechanism
  is co-translational); they are included only for vocabulary
  symmetry and the manifest annotates them as inert
  (``inert_in_peptide_mode=True``).

Every entry carries a citation in the dataclass ``citation`` field.
Do not introduce new entries without one.

Aliases (``ALIASES``) map informal historical names to the canonical
form so that ``--peptide-linker GS3`` and ``--peptide-linker G4S3``
both resolve to the same Linker object.
"""

import logging
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Linker:
    """A linker entry for vaccine construct assembly."""
    name: str
    amino_acids: str
    dna: Optional[str] = None
    freeze_in_mrna: bool = False
    inert_in_peptide_mode: bool = False
    citation: str = ""


# -- GnSm flexible linker family --------------------------------------------
#
# Canonical form is GnSm — n glycines followed by m serines. The (G4S)k
# repeated form is the dominant choice in vaccine and antibody literature
# (Huston et al., PNAS 85:5879, 1988, doi:10.1073/pnas.85.16.5879;
# reviewed in Chen et al., Adv Drug Deliv Rev 65:1357, 2013,
# doi:10.1016/j.addr.2012.09.039). Shorter / longer Gn-S variants give
# different stiffness/length trade-offs.

_GS_CITATION = (
    "Huston et al., PNAS 85:5879, 1988 (doi:10.1073/pnas.85.16.5879) — "
    "(Gly4Ser)n linker family, originally connecting VH-VL in scFv. "
    "Length variants (G2S/G3S/G5S, repeats up to (G4S)4) reviewed in "
    "Chen et al., Adv Drug Deliv Rev 65:1357, 2013 "
    "(doi:10.1016/j.addr.2012.09.039)."
)

LINKER_G2S = Linker(name="G2S", amino_acids="GGS", citation=_GS_CITATION)
LINKER_G3S = Linker(name="G3S", amino_acids="GGGS", citation=_GS_CITATION)
LINKER_G4S = Linker(name="G4S", amino_acids="GGGGS", citation=_GS_CITATION)
LINKER_G5S = Linker(name="G5S", amino_acids="GGGGGS", citation=_GS_CITATION)
LINKER_G4S2 = Linker(name="G4S2", amino_acids="GGGGS" * 2, citation=_GS_CITATION)
LINKER_G4S3 = Linker(name="G4S3", amino_acids="GGGGS" * 3, citation=_GS_CITATION)
LINKER_G4S4 = Linker(name="G4S4", amino_acids="GGGGS" * 4, citation=_GS_CITATION)


# -- Rigid alpha-helical linkers --------------------------------------------

_EAAAK_CITATION = (
    "Arai et al., Protein Eng 14:529, 2001 (doi:10.1093/protein/14.8.529) — "
    "(Glu-Ala-Ala-Ala-Lys)n alpha-helical linker; rigid spacer used when "
    "antigen domains need separation rather than flex."
)

LINKER_EAAAK = Linker(name="EAAAK", amino_acids="EAAAK", citation=_EAAAK_CITATION)
LINKER_EAAAK2 = Linker(
    name="EAAAK2", amino_acids="EAAAK" * 2, citation=_EAAAK_CITATION)
LINKER_EAAAK3 = Linker(
    name="EAAAK3", amino_acids="EAAAK" * 3, citation=_EAAAK_CITATION)


# -- Furin cleavage sites ----------------------------------------------------
#
# Cleaved in the trans-Golgi network by the proprotein convertase furin;
# used in some multi-antigen mRNA designs as an alternative to 2A
# ribosomal skipping. Consensus motif is R-X-(K/R)-R; the three sites
# below are well-characterized minimal forms.

_FURIN_CITATION = (
    "Thomas G., Nat Rev Mol Cell Biol 3:753, 2002 (doi:10.1038/nrm934) — "
    "Furin cleavage at R-X-(K/R)-R motifs in the trans-Golgi network. "
    "The minimal R-X-K-R motifs (RKRR / RVKR / RKRKR) are commonly used "
    "to liberate concatenated antigens in multi-cistronic constructs."
)

LINKER_FURIN_RKRR = Linker(
    name="RKRR", amino_acids="RKRR", citation=_FURIN_CITATION)
LINKER_FURIN_RVKR = Linker(
    name="RVKR", amino_acids="RVKR", citation=_FURIN_CITATION)
LINKER_FURIN_RKRKR = Linker(
    name="RKRKR", amino_acids="RKRKR", citation=_FURIN_CITATION)


# -- MHC-epitope spacers ----------------------------------------------------

LINKER_AAY = Linker(
    name="AAY",
    amino_acids="AAY",
    citation=(
        "Velders et al., J Immunol 166:5366, 2001 "
        "(doi:10.4049/jimmunol.166.9.5366) — established that AAY supports "
        "proper proteasomal/TAP processing of joined CTL epitope strings "
        "in DNA-vaccine constructs."),
)

LINKER_GPGPG = Linker(
    name="GPGPG",
    amino_acids="GPGPG",
    citation=(
        "Livingston et al., J Immunol 168:5499, 2002 "
        "(doi:10.4049/jimmunol.168.11.5499) — demonstrated that GPGPG "
        "breaks junctional HLA-DR epitopes while preserving each helper "
        "epitope when assembling multi-epitope constructs."),
)


# -- 2A self-cleaving peptides ----------------------------------------------
#
# 2A skipping is co-translational and depends on the C-terminal
# D(V/I)E(E/S)NPGP motif (Doronina et al., MCB 28:4227, 2008,
# doi:10.1128/MCB.00421-08; Sharma et al., NAR 40:3143, 2012,
# doi:10.1093/nar/gkr1176). Empirically, codon-optimized 2A variants
# work in mammalian cells (Kim et al., PLoS ONE 6:e18556, 2011,
# doi:10.1371/journal.pone.0018556; Liu et al., Sci Rep 7:2193, 2017,
# doi:10.1038/s41598-017-02460-2). Vaxrank still ships these with
# freeze_in_mrna=True as a conservative default — published constructs
# all use canonical codons.
#
# All four entries include the GSG N-terminal spacer recommended by
# Szymczak et al., Nat Biotechnol 22:589, 2004 (doi:10.1038/nbt957) and
# shown to improve cleavage efficiency by Wang et al., Sci Rep 5:16273,
# 2015 (doi:10.1038/srep16273).

LINKER_P2A = Linker(
    name="P2A",
    amino_acids="GSGATNFSLLKQAGDVEENPGP",
    dna=(
        "GGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGA"
        "GACGTGGAGGAGAACCCTGGACCT"
    ),
    freeze_in_mrna=True,
    inert_in_peptide_mode=True,
    citation=(
        "AA: Donnelly et al., J Gen Virol 82:1027, 2001 "
        "(doi:10.1099/0022-1317-82-5-1027). DNA: Kim et al., "
        "PLoS ONE 6:e18556, 2011 (doi:10.1371/journal.pone.0018556). "
        "GSG spacer: Szymczak et al., Nat Biotechnol 22:589, 2004."),
)

LINKER_T2A = Linker(
    name="T2A",
    amino_acids="GSGEGRGSLLTCGDVEENPGP",
    dna=(
        "GGAAGCGGAGAGGGCAGAGGAAGTCTGCTAACATGCGGT"
        "GACGTCGAGGAGAATCCTGGACCT"
    ),
    freeze_in_mrna=True,
    inert_in_peptide_mode=True,
    citation=(
        "AA: Donnelly et al., J Gen Virol 82:1013, 2001 "
        "(doi:10.1099/0022-1317-82-5-1013). DNA: Kim et al., "
        "PLoS ONE 6:e18556, 2011."),
)

LINKER_F2A = Linker(
    name="F2A",
    amino_acids="GSGVKQTLNFDLLKLAGDVESNPGP",
    dna=(
        "GGAAGCGGAGTGAAACAGACTTTGAATTTTGACCTTCTCAAG"
        "TTGGCAGGAGACGTTGAGTCCAACCCTGGGCCT"
    ),
    freeze_in_mrna=True,
    inert_in_peptide_mode=True,
    citation=(
        "AA: Ryan & Drew, EMBO J 13:928, 1994 "
        "(doi:10.1002/j.1460-2075.1994.tb06337.x); Donnelly et al. 2001. "
        "DNA: Kim et al., PLoS ONE 6:e18556, 2011."),
)

LINKER_E2A = Linker(
    name="E2A",
    amino_acids="GSGQCTNYALLKLAGDVESNPGP",
    dna=(
        "GGAAGCGGACAGTGTACTAATTATGCTCTCTTGAAATTG"
        "GCTGGAGATGTTGAGAGCAACCCTGGACCT"
    ),
    freeze_in_mrna=True,
    inert_in_peptide_mode=True,
    citation=(
        "AA: Donnelly et al., J Gen Virol 82:1013, 2001. "
        "DNA: Kim et al., PLoS ONE 6:e18556, 2011."),
)


# Canonical name → Linker. CLI normalizes input to upper-case and
# dereferences ALIASES first, so users can use either the canonical or
# the historical name.
LINKERS = {
    LINKER_G2S.name: LINKER_G2S,
    LINKER_G3S.name: LINKER_G3S,
    LINKER_G4S.name: LINKER_G4S,
    LINKER_G5S.name: LINKER_G5S,
    LINKER_G4S2.name: LINKER_G4S2,
    LINKER_G4S3.name: LINKER_G4S3,
    LINKER_G4S4.name: LINKER_G4S4,
    LINKER_EAAAK.name: LINKER_EAAAK,
    LINKER_EAAAK2.name: LINKER_EAAAK2,
    LINKER_EAAAK3.name: LINKER_EAAAK3,
    LINKER_FURIN_RKRR.name: LINKER_FURIN_RKRR,
    LINKER_FURIN_RVKR.name: LINKER_FURIN_RVKR,
    LINKER_FURIN_RKRKR.name: LINKER_FURIN_RKRKR,
    LINKER_AAY.name: LINKER_AAY,
    LINKER_GPGPG.name: LINKER_GPGPG,
    LINKER_P2A.name: LINKER_P2A,
    LINKER_T2A.name: LINKER_T2A,
    LINKER_F2A.name: LINKER_F2A,
    LINKER_E2A.name: LINKER_E2A,
}

# Historical or convenience names → canonical name.
ALIASES = {
    "GS": "G4S",         # 2.10.0 default; matches G4S in standard nomenclature
    "GS3": "G4S3",       # 2.10.0 default; matches (G4S)3
}


def all_linker_names():
    """Canonical names + aliases. Use as argparse ``choices=`` value."""
    return sorted(set(LINKERS) | set(ALIASES))


def get_linker(name):
    """Resolve a linker by name (after de-aliasing); raise ValueError on miss."""
    canonical = ALIASES.get(name, name)
    if canonical not in LINKERS:
        raise ValueError(
            "Unknown linker '%s'. Available: %s" % (
                name, ', '.join(all_linker_names())))
    return LINKERS[canonical]


def select_antigen_window(fragment, base_name, max_length_aa):
    """Return a sub-window of the vaccine peptide that preserves the mutation.

    A naive ``amino_acids[:max_length_aa]`` truncation can drop the
    mutated residues entirely when the mutation sits past the head of
    the fragment. This helper centers a window of ``max_length_aa`` on
    the mutation (using ``mutant_amino_acid_*_offset`` from the fragment).
    If the mutation itself exceeds the cap (rare; long inframe
    insertions), the full fragment is emitted unchanged with a warning.

    Used by both ``vaxrank.peptide`` (SLP mode) and ``vaxrank.mrna``
    (per-antigen window in concatenated constructs).
    """
    aa = fragment.amino_acids
    if len(aa) <= max_length_aa:
        return aa
    mut_start = getattr(fragment, 'mutant_amino_acid_start_offset', 0)
    mut_end = getattr(fragment, 'mutant_amino_acid_end_offset', len(aa))
    mut_len = mut_end - mut_start
    if mut_len > max_length_aa:
        logger.warning(
            "Mutation in %s spans %d aa, longer than the antigen-length "
            "cap (%d); emitting full fragment without truncation.",
            base_name, mut_len, max_length_aa)
        return aa
    midpoint = (mut_start + mut_end) // 2
    half = max_length_aa // 2
    start = max(0, midpoint - half)
    end = min(len(aa), start + max_length_aa)
    start = max(0, end - max_length_aa)
    return aa[start:end]
