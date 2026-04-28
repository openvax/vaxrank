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

"""Shared linker library for mRNA and peptide vaccine constructs.

A linker is anything that goes between concatenated antigens. Some
linkers are pure amino-acid spacers (GS-family, AAY, GPGPG); 2A
self-cleaving peptides also live here because the *vocabulary* is
shared even though the modality semantics differ:

- In mRNA constructs, 2A linkers cause ribosomal skipping; their
  blessed DNA is preserved during codon optimization
  (``freeze_in_mrna=True``).
- In synthesized peptides, 2A motifs do not cleave (the mechanism
  is co-translational); they are included only for vocabulary
  symmetry and the manifest annotates them as inert
  (``inert_in_peptide_mode=True``).

Every entry carries a citation in its docstring or in the dataclass
``citation`` field. Do not introduce new entries without one.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Linker:
    """A linker entry for vaccine construct assembly."""
    name: str
    amino_acids: str
    dna: Optional[str] = None
    freeze_in_mrna: bool = False
    inert_in_peptide_mode: bool = False
    citation: str = ""


# -- Pure amino-acid linkers (back-translated freely in mRNA) ----------------

LINKER_GS = Linker(
    name="GS",
    amino_acids="GGGGS",
    citation=(
        "Huston et al., PNAS 85:5879, 1988 (doi:10.1073/pnas.85.16.5879) — "
        "first description of the (Gly4Ser)n linker family connecting VH-VL "
        "in single-chain Fv. Reviewed in Chen et al., Adv Drug Deliv Rev "
        "65:1357, 2013 (doi:10.1016/j.addr.2012.09.039)."),
)

LINKER_GS3 = Linker(
    name="GS3",
    amino_acids="GGGGSGGGGSGGGGS",
    citation=LINKER_GS.citation,
)

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


LINKERS = {
    LINKER_GS.name: LINKER_GS,
    LINKER_GS3.name: LINKER_GS3,
    LINKER_AAY.name: LINKER_AAY,
    LINKER_GPGPG.name: LINKER_GPGPG,
    LINKER_P2A.name: LINKER_P2A,
    LINKER_T2A.name: LINKER_T2A,
    LINKER_F2A.name: LINKER_F2A,
    LINKER_E2A.name: LINKER_E2A,
}


def get_linker(name):
    """Resolve a linker by name; raise ValueError on miss."""
    if name not in LINKERS:
        raise ValueError(
            "Unknown linker '%s'. Available: %s" % (
                name, ', '.join(sorted(LINKERS))))
    return LINKERS[name]
