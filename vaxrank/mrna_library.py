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

Note: the BioNTech FixVac "FI" element uses a tandem 2x HBB 3' UTR
variant — see Orlandini von Niessen et al., Mol Ther 27:824, 2019
(doi:10.1016/j.ymthe.2018.12.011). Not yet provided here; users who
want it can pass a custom 3' UTR string.
"""


# -- Signal peptides ---------------------------------------------------------

SIGNAL_PEPTIDE_HLA_A = "MAVMAPRTLLLLLSGALALTQTWA"
"""HLA-A signal peptide, 24 aa.

Source: UniProt P04439 (HLA class I A alpha chain), signal peptide
annotation residues 1-24. This is the canonical N-terminal signal
paired with MITD in BioNTech antigen-MITD constructs (Kreiter et al.,
J Immunol 180:309, 2008, doi:10.4049/jimmunol.180.1.309;
Sahin et al., Nature 547:222, 2017).
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
Originally described in Kreiter et al., J Immunol 180:309, 2008
(doi:10.4049/jimmunol.180.1.309) for routing antigens to the
endolysosomal compartment to enhance MHC class I and class II
presentation. Used in BioNTech IVAC MUTANOME / FixVac trials
(Sahin et al., Nature 547:222, 2017, doi:10.1038/nature23003).
"""


# -- Lookup tables -----------------------------------------------------------

UTRS_5P = {
    "HBB": UTR_5P_HBB,
}

UTRS_3P = {
    "HBB": UTR_3P_HBB,
}

SIGNAL_PEPTIDES = {
    "HLA_A": SIGNAL_PEPTIDE_HLA_A,
    "tPA": SIGNAL_PEPTIDE_TPA,
    "IgK": SIGNAL_PEPTIDE_IGK,
    "CD8A": SIGNAL_PEPTIDE_CD8A,
    "CD28": SIGNAL_PEPTIDE_CD28,
}

MITDS = {
    "HLA_A": MITD_HLA_A,
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
