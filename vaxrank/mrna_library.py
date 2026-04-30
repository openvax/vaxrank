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

"""Reference sequences for mRNA vaccine construct assembly.

Every sequence in this module is sourced from a primary database
(RefSeq, UniProt) or peer-reviewed publication. The accession or DOI
is recorded in the docstring next to each constant. Do not introduce
new entries without a citation.

Conventions
-----------
- ``UTR_5P_*`` / ``UTR_3P_*`` constants are RNA in DNA spelling
  (T not U), matching how varcode/isovar carry sequences.
- ``SIGNAL_PEPTIDE_*``, ``MITD_*``, ``LINKER_*`` are amino-acid
  strings — they get back-translated and codon-optimized at assembly time.
"""

# -- 5' / 3' UTRs ------------------------------------------------------------

UTR_5P_HBB = "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC"
"""Human β-globin (HBB) 5' UTR, 50 nt.

Source: NCBI RefSeq NM_000518.5, bases 1-50 (CDS starts at base 51).
Used in BioNTech FixVac and many published mRNA vaccine constructs
(Sahin et al., Nature 547:222, 2017, doi:10.1038/nature23003;
Holtkamp et al., Blood 108:4009, 2006, doi:10.1182/blood-2006-04-015024).
"""

UTR_3P_HBB = (
    "GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTAC"
    "TAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATT"
    "TATTTTCATTGCAA"
)
"""Human β-globin (HBB) 3' UTR, 134 nt.

Source: NCBI RefSeq NM_000518.5, bases 495-628. Original sequence:
Proudfoot, Cell 10:559, 1977 (PMID 67897).
"""

UTR_3P_HBB_FI = UTR_3P_HBB + UTR_3P_HBB
"""Human β-globin tandem 2× 3' UTR ("2hBg" / FI element), 268 nt.

This is the 3' UTR architecture used in BioNTech IVAC MUTANOME /
FixVac mRNAs as published in Sahin et al., Nature 547:222, 2017
(doi:10.1038/nature23003) and originally optimized in Holtkamp
et al., Blood 108:4009, 2006 (doi:10.1182/blood-2006-04-015024).
Two head-to-tail copies of the canonical HBB 3' UTR.

Note: BioNTech's later platform (BNT162b2 COVID, post-Orlandini
von Niessen et al., Mol Ther 27:824, 2019, doi:10.1016/j.ymthe.2018.12.011)
migrated to an AES + mtRNR1 combination that outperformed 2hBg in
their cellular library screen. AES + mtRNR1 sequences are not yet
bundled here.
"""


# -- Signal peptides ---------------------------------------------------------

SIGNAL_PEPTIDE_HLA_A = "MAVMAPRTLLLLLSGALALTQTWA"
"""HLA-A signal peptide, 24 aa.

Source: UniProt P04439 (HLA class I A alpha chain), signal peptide
annotation residues 1-24. Originally cited as the BioNTech FixVac
N-terminal signal (per Kreiter et al., J Immunol 180:309, 2008,
doi:10.4049/jimmunol.180.1.309), but Kreiter 2008 itself does not
disclose the allele. Sahin et al., Nature 547:222, 2017
(doi:10.1038/nature23003) reports an HLA-B-family leader paired
with the HLA-B-family MITD; see ``SIGNAL_PEPTIDE_HLA_B``.
"""

SIGNAL_PEPTIDE_HLA_B = "MLVMAPRTVLLLLSAALALTETWA"
"""HLA-B signal peptide, 24 aa.

Source: UniProt P01889 (HLA class I B alpha chain), signal peptide
annotation residues 1-24. This is the leader BioNTech FixVac /
IVAC MUTANOME pairs with the HLA-B-family MITD per Sahin et al.,
Nature 547:222, 2017 (doi:10.1038/nature23003). Slight residue
differences from the Sahin 2017 paper's quoted sequence are likely
HLA-B subtype variation; this is the canonical UniProt P01889 form.
"""

SIGNAL_PEPTIDE_TPA = "MDAMKRGLCCVLLLCGAVFVSPS"
"""Tissue plasminogen activator (tPA) signal peptide, 23 aa.

Source: UniProt P00750 (PLAT_HUMAN), signal peptide annotation
residues 1-23. Widely used as a heterologous secretion signal in
DNA/mRNA vaccines (Luke et al., Vaccine 27:6911, 2009,
doi:10.1016/j.vaccine.2009.09.005).
"""

SIGNAL_PEPTIDE_IGK = "METDTLLLWVLLLWVPGSTG"
"""Mouse Ig kappa light-chain signal peptide, 20 aa.

Source: Coloma et al., J Immunol Methods 152:89, 1992 (PMID 1640114);
GenBank V00777 (mouse Ig kappa V-J2-C). Distributed in pSecTag /
pVAC1 expression vectors and used in many DNA-vaccine constructs
(Kou et al., Sci Rep 7:42408, 2017, doi:10.1038/srep42408).
"""

SIGNAL_PEPTIDE_CD8A = "MALPVTALLLPLALLLHAARP"
"""Human CD8A signal peptide, 21 aa.

Source: UniProt P01732 (CD8A_HUMAN), signal peptide annotation
residues 1-21. Common choice in TCR-engineering and cell-surface
expression constructs.
"""

SIGNAL_PEPTIDE_CD28 = "MLRLLLALNLFPSIQVTG"
"""Human CD28 signal peptide, 18 aa.

Source: UniProt P10747 (CD28_HUMAN), signal peptide annotation
residues 1-18. Pairs with CD8A in dual-chain TCR constructs to
distinguish alpha/beta secretion signals.
"""


# -- MHC-I trafficking domain -----------------------------------------------

MITD_HLA_A = "IVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSLTACKV"
"""MHC class I trafficking domain (MITD) from HLA-A, 58 aa.

Source: UniProt P04439 (HLA class I A alpha chain), residues 308-365
(transmembrane helix 309-332 + cytoplasmic domain 333-365).
Mechanistically described in Kreiter et al., J Immunol 180:309, 2008
(doi:10.4049/jimmunol.180.1.309) for routing antigens to the
endolysosomal compartment to enhance MHC class I and class II
presentation. Note: Kreiter 2008 does not disclose the specific HLA
allele; ``MITD_HLA_B`` is what BioNTech FixVac (Sahin et al., Nature
547:222, 2017) actually uses, paired with the HLA-B signal peptide.
"""

MITD_HLA_B = "IVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSLTA"
"""MHC class I trafficking domain (MITD) from HLA-B, 54 aa.

Source: UniProt P01889 (HLA class I B alpha chain), residues 309-362
(transmembrane helix 309-332 + cytoplasmic domain 333-362). This is
the MITD BioNTech FixVac / IVAC MUTANOME pairs with the HLA-B signal
peptide per Sahin et al., Nature 547:222, 2017
(doi:10.1038/nature23003). Diagnostic differences from HLA-A:
``IVGIV`` (B) vs ``IVGII`` (A); ``CRRKSSGGK`` (B) vs ``WRRKSSDRK``
(A); HLA-B's tail ends ``...VSLTA`` (no CKV).

Sahin 2017's quoted sequence has minor residue differences from the
canonical UniProt P01889 — those are likely HLA-B subtype variation.
This entry uses the canonical UniProt form.
"""


# -- Lookup tables -----------------------------------------------------------

UTRS_5P = {
    "HBB": UTR_5P_HBB,
}

UTRS_3P = {
    "HBB": UTR_3P_HBB,
    "HBB_FI": UTR_3P_HBB_FI,  # tandem 2× HBB; BioNTech FixVac canonical
}

SIGNAL_PEPTIDES = {
    "HLA_A": SIGNAL_PEPTIDE_HLA_A,
    "HLA_B": SIGNAL_PEPTIDE_HLA_B,  # BioNTech FixVac canonical
    "tPA": SIGNAL_PEPTIDE_TPA,
    "IgK": SIGNAL_PEPTIDE_IGK,
    "CD8A": SIGNAL_PEPTIDE_CD8A,
    "CD28": SIGNAL_PEPTIDE_CD28,
}

MITDS = {
    "HLA_A": MITD_HLA_A,
    "HLA_B": MITD_HLA_B,  # BioNTech FixVac canonical
}


# -- Linkers (string compatibility shim) -------------------------------------
#
# 2.10.0 exposed linker entries as plain amino-acid strings via this
# module: ``from vaxrank.mrna_library import LINKER_GS3``. The richer
# ``Linker`` dataclass (with optional blessed DNA, freeze flags, and
# citations) lives in ``vaxrank.vaccine_library``; for callers that
# imported the string form, we keep the originals here unchanged.

LINKER_GS = "GGGGS"
LINKER_GS3 = "GGGGSGGGGSGGGGS"
LINKER_AAY = "AAY"
LINKER_GPGPG = "GPGPG"

LINKERS = {
    "GS": LINKER_GS,
    "GS3": LINKER_GS3,
    "AAY": LINKER_AAY,
    "GPGPG": LINKER_GPGPG,
}
