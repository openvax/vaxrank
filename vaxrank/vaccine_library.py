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

## Compositional grammar

In addition to named entries, ``get_linker`` accepts compositional
forms parsed at lookup time. **Repeats require explicit parens or
an ``x`` separator** — bare digit suffixes are parsed literally to
avoid ambiguity (e.g. ``G4S2`` = "GGGGSS", *not* "(G4S)2").

- **Repeat forms** — ``(BASE)N``, ``(BASE)xN``, ``BASExN``. Examples:
  ``(G4S)2``, ``(G4S)x2``, and ``G4Sx2`` all → "GGGGSGGGGS".
  2A entries (codon-frozen) are rejected — repeating a 2A linker
  would not produce additional ribosomal-skipping events.
- ``GnSm`` — n glycines followed by m serines, **as a literal single
  unit, not a repeat**. Examples: ``G4S`` = "GGGGS"; ``G4S2`` =
  "GGGGSS"; ``G6S`` = "GGGGGGS". Use ``(G4S)2`` for the (Gly4Ser)2
  repeat from Huston 1988 / Chen 2013.
- ``AnY`` — n alanines followed by tyrosine. Example: ``A3Y`` →
  "AAAY". AAY (= A2Y) is the conventional form (Wang 2004 / Velders
  2001); n>2 has no primary-literature footing — pure mechanistic
  extrapolation. The synthesized Linker's ``citation`` flags this.
- ``An`` — n alanines, no tyrosine. Example: ``A4`` → "AAAA". AAA
  is the empirically-tested form (Aguilar-Gurrieri 2023, the only
  published alanine-spacer bake-off); longer An are length
  extrapolations.
- ``Gn`` — n glycines, no serine. Example: ``G4`` → "GGGG". Pure
  polyglycine has biophysical characterization (Klement 2018) as the
  most-flexible linker but NO published vaccine empirical use. The
  citation flags this — prefer GS-family forms (G4S, (G4S)2, etc.)
  for any expression-dependent platform.

Repeat counts are capped at 100 to prevent accidental megasequences.
"""

import logging
import re
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
    "Huston et al., PNAS 85:5879, 1988 (doi:10.1073/pnas.85.16.5879) "
    "established the (Gly4Ser)n family for connecting VH-VL in scFv. "
    "Length variants (G2S/G3S/G5S and shorter/longer repeats) appear "
    "throughout subsequent scFv / fusion-protein literature; reviewed "
    "in Chen et al., Adv Drug Deliv Rev 65:1357, 2013 "
    "(doi:10.1016/j.addr.2012.09.039). Klement et al., Biochemistry "
    "57:1378, 2018 (doi:10.1021/acs.biochem.7b00902) showed persistence "
    "length scales with glycine fraction (more Gly = more flexible). "
    "Clinical mRNA neoantigen vaccine use: BioNTech IVAC MUTANOME / "
    "FixVac (Sahin et al., Nature 547:222, 2017, "
    "doi:10.1038/nature23003) and iNeST / autogene cevumeran (Rojas "
    "et al., Nature 618:144, 2023, doi:10.1038/s41586-023-06063-y) use "
    "a 10-aa Gly/Ser linker — i.e. (G4S)2 — between 25-27mer minigenes. "
    "This is the only linker family with published clinical mRNA "
    "vaccine use as of 2025. Caveats: empirical vaccine data on the "
    "GS family is thin — Yang 2015 (PMC4514284) showed GGGS worked "
    "while AAY failed in HIV DNA vaccines; Aguilar-Gurrieri 2023 "
    "(doi:10.1007/s00262-023-03409-3) showed GGGS performed WORST for "
    "MHC-I presentation among AAA / AAL / ADL / single-A / GGGS. No "
    "published vaccine study has varied (G4S)n length and measured "
    "immunogenicity; length effects are inherited from antibody / "
    "scFv engineering."
)

_POLYG_CITATION = (
    "Pure polyglycine (Gn, no serine) has biophysical characterization "
    "as the most-flexible linker family (lowest persistence length; "
    "Klement et al., Biochemistry 57:1378, 2018, "
    "doi:10.1021/acs.biochem.7b00902 measured persistence length 4.5 Å "
    "for full-Gly vs 6.2 Å for full-Ser). However, polyglycine has NO "
    "published primary use as an inter-epitope spacer in any vaccine "
    "modality (peptide, DNA, mRNA, viral vector) with an immunogenicity "
    "readout. The choice trade-off vs (G4S)n: more flexibility but "
    "lower solubility (no polar residue). Use at your own design risk; "
    "consider AAA (Aguilar-Gurrieri 2023, the only empirically-tested "
    "alanine-only spacer for vaccines) as an alternative."
)

# Single-unit GnSm entries only. Repeats use the compositional grammar:
# (G4S)2, (G4S)x2, or G4Sx2 — never bare 'G4S2' (which now parses as
# the literal "GGGGSS"; see module docstring).
LINKER_G2S = Linker(name="G2S", amino_acids="GGS", citation=_GS_CITATION)
LINKER_G3S = Linker(name="G3S", amino_acids="GGGS", citation=_GS_CITATION)
LINKER_G4S = Linker(name="G4S", amino_acids="GGGGS", citation=_GS_CITATION)
LINKER_G5S = Linker(name="G5S", amino_acids="GGGGGS", citation=_GS_CITATION)


# -- Rigid alpha-helical linkers --------------------------------------------

_EAAAK_CITATION = (
    "Arai et al., Protein Eng 14:529, 2001 (doi:10.1093/protein/14.8.529) — "
    "(Glu-Ala-Ala-Ala-Lys)n alpha-helical linker; rigid spacer used when "
    "antigen domains need separation rather than flex."
)

# Single-unit only; repeats use (EAAAK)n / (EAAAK)xn / EAAAKxn.
LINKER_EAAAK = Linker(name="EAAAK", amino_acids="EAAAK", citation=_EAAAK_CITATION)


# -- Furin cleavage sites ----------------------------------------------------
#
# Cleaved in the trans-Golgi network by the proprotein convertase furin;
# used in some multi-antigen mRNA designs as an alternative to 2A
# ribosomal skipping. Consensus motif is R-X-(K/R)-R; the three sites
# below are well-characterized minimal forms.

_FURIN_CITATION = (
    "Furin cleaves at the consensus motif R-X-(K/R)-R in the trans-Golgi "
    "network. Primary establishment of the consensus: Hosaka et al., "
    "J Biol Chem 266:12127, 1991 (PMID 1905715) — 'Arg-X-Lys/Arg-Arg "
    "motif as a signal for precursor cleavage catalyzed by furin within "
    "the constitutive secretory pathway'. Reviewed in Thomas G., "
    "Nat Rev Mol Cell Biol 3:753, 2002 (doi:10.1038/nrm934). "
    "The minimal motifs (RKRR / RVKR / RKRKR) appear in multi-cistronic "
    "expression and DNA-vaccine constructs as alternatives to 2A "
    "ribosomal skipping. Cho & Celis, Cancer Immunol Immunother "
    "61:343, 2012 (PMC4019994) tested RVKR vs furin-resistant VRVV in "
    "a preclinical melanoma DNA vaccine and observed equivalent CD8 "
    "responses. As of 2025 NO published vaccine clinical trial of any "
    "modality (mRNA, DNA, peptide, viral vector) uses an engineered "
    "furin site as an inter-antigen linker; native furin sites in "
    "single antigens (HIV gp160→gp120/gp41; SARS-CoV-2 spike S1/S2; "
    "influenza HA0) are not the same thing. Furin-T2A cassettes are "
    "common in CAR-T cell therapies but those are cell therapies, "
    "not vaccines."
)

LINKER_FURIN_RKRR = Linker(
    name="RKRR", amino_acids="RKRR", citation=_FURIN_CITATION)
LINKER_FURIN_RVKR = Linker(
    name="RVKR", amino_acids="RVKR", citation=_FURIN_CITATION)
LINKER_FURIN_RKRKR = Linker(
    name="RKRKR", amino_acids="RKRKR", citation=_FURIN_CITATION)


# -- MHC-epitope spacers ----------------------------------------------------

_AAY_CITATION = (
    "Empirical mechanistic foundation: Livingston et al., Vaccine "
    "19:4652-4660, 2001 (PMID 11535313, doi:10.1016/S0264-410X(01)00233-X) "
    "— the Sette/Epimmune group built a 94-epitope flanking-residue "
    "database and demonstrated that small/amide/basic residues at C+1 "
    "(immediately after the CTL epitope) modulate processing and "
    "immunogenicity; Ala and Tyr at C+1 gave the best CTL responses in "
    "HLA-transgenic mice. This is the foundational basis for AAY (Ala "
    "at P-2/P-1, Tyr at C+1). Velders et al., J Immunol 166:5366, 2001 "
    "(doi:10.4049/jimmunol.166.9.5366) showed spacer-vs-no-spacer in HPV16 "
    "DNA vaccines but did not compare AAY against other linkers. Wang et "
    "al., Vaccine 22:3622, 2004 (PMID 15320877) used AAY in a TB DNA "
    "vaccine. Subsequent immunoinformatics designs adopted AAY by "
    "convention rather than head-to-head comparison. Cho & Celis, "
    "Cancer Immunol Immunother 61:343, 2012 (PMC4019994) explicitly "
    "noted AAY remained empirically unvalidated against alternatives. "
    "Important caveat: Yang et al., Hum Vaccin Immunother 11:795, 2015 "
    "(PMC4514284) tested AAY vs GGGS in HIV-1 multi-epitope DNA vaccines "
    "(MEG3 vs MEG2; one construct each, n=1 per linker, identical "
    "epitope order). MEG3 (AAY) produced no detectable Western blot "
    "signal in 293T, no IgG response, and ~20× lower ELISpot than the "
    "GGGS constructs. Authors propose two mechanisms: (a) AAY changed "
    "predicted structure (less flex, more alpha-helix) impairing "
    "expression, (b) proteasomal cleavage at AAY happened too rapidly "
    "during translation. Hypothesis (b) is paradoxical — fast proteasomal "
    "release is precisely what AAY is designed for and should ENHANCE "
    "MHC presentation, so (a) is more consistent with the observed low "
    "ELISpot. Either mechanism would apply equally to mRNA constructs "
    "(same translation, same proteasome, same MHC pathway as DNA "
    "vaccines); only synthesized peptide vaccines (no expression "
    "dependency) sidestep these failure modes. No equivalent AAY-vs-GS "
    "head-to-head has been published in an mRNA-vaccine context, so "
    "Yang 2015's result should be treated as a real warning about AAY "
    "in any expression-dependent platform — not just DNA. The strongest "
    "empirical alanine-spacer data is "
    "Aguilar-Gurrieri et al., Cancer Immunol Immunother 72:2113, 2023 "
    "(doi:10.1007/s00262-023-03409-3), which tested AAA / AAL / ADL / "
    "single-A / GGGS for MHC-I presentation in a polypeptide context — "
    "AAA won, AAY itself was NOT tested. The mechanistic claim that Y "
    "is THE P1 proteasome anchor is oversold: Toes et al., J Exp Med "
    "194:1, 2001 (PMC2193442) shows L > F at P1 with Y in the same "
    "hydrophobic cluster but not specifically preferred."
)

_ANY_CITATION = (
    "AnY variants beyond A2Y (= AAY) have NO primary-literature support. "
    "AAAY / A4Y / A5Y appear in immunoinformatics design papers but no "
    "primary publication uses these as deliberately-chosen spacers with "
    "empirical validation. They are mechanistic extrapolations combining "
    "(a) Aguilar-Gurrieri et al. 2023 (doi:10.1007/s00262-023-03409-3) "
    "showing polyalanine spacers outperform GGGS for MHC-I presentation, "
    "untested at length n>3, with (b) the conventional but oversold claim "
    "that Y is the canonical P1 proteasome anchor (Toes et al. 2001, "
    "PMC2193442 actually shows L > F > Y). For an empirical alanine "
    "linker with published data, use AAA (LINKER_AAA in this library)."
)

_ALANINE_CITATION = (
    "Aguilar-Gurrieri et al., Cancer Immunol Immunother 72:2113, 2023 "
    "(doi:10.1007/s00262-023-03409-3, PMC10264286) — tested AAA / AAL / "
    "ADL / single-A / GGGS as inter-epitope spacers in a 13-neoantigen + "
    "SIINFEKL polypeptide construct. AAA gave the strongest H-2Kb/SIINFEKL "
    "surface presentation by 25-D1.16 staining; GGGS was worst. Sole "
    "published empirical bake-off of alanine-based spacers as of 2025. "
    "Did NOT test AAY, AAF, AAW, or longer An (AAAA, AAAAA). Single-A "
    "underperformed AAA, suggesting a length floor; no upper-bound data."
)

LINKER_AAY = Linker(
    name="AAY",
    amino_acids="AAY",
    citation=_AAY_CITATION,
)

LINKER_AAA = Linker(
    name="AAA",
    amino_acids="AAA",
    citation=_ALANINE_CITATION,
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
    LINKER_EAAAK.name: LINKER_EAAAK,
    LINKER_FURIN_RKRR.name: LINKER_FURIN_RKRR,
    LINKER_FURIN_RVKR.name: LINKER_FURIN_RVKR,
    LINKER_FURIN_RKRKR.name: LINKER_FURIN_RKRKR,
    LINKER_AAY.name: LINKER_AAY,
    LINKER_AAA.name: LINKER_AAA,
    LINKER_GPGPG.name: LINKER_GPGPG,
    LINKER_P2A.name: LINKER_P2A,
    LINKER_T2A.name: LINKER_T2A,
    LINKER_F2A.name: LINKER_F2A,
    LINKER_E2A.name: LINKER_E2A,
}

# Historical or convenience names → canonical form. The right-hand side
# is fed back through get_linker, so it can be a static name or a
# compositional form.
ALIASES = {
    "GS": "G4S",          # 2.10.0 default
    "GS3": "(G4S)3",      # 2.10.0 default; was the (G4S)3 repeat
}


# Default candidate set for the per-junction linker optimizer
# (see vaxrank.junction_swap). All five resolve via the compositional
# grammar; mechanistically diverse (length variants of GS-family +
# alanine spacer) without straying into 2A / EAAAK / furin which have
# different mechanisms or no clinical track record.
JUNCTION_SWAP_CANDIDATES = ("G3S", "G4S", "(G3S)2", "(G4S)2", "AAA")


# Compositional-name parsers (see module docstring for grammar).
# Repeats require explicit parens or 'x' separator; bare digit suffixes
# parse literally (G4S2 = "GGGGSS", not "(G4S)2").
_PAREN_REPEAT_RE = re.compile(r"^\((?P<base>[A-Z][A-Z0-9]*)\)X?(?P<count>\d+)$")
_X_REPEAT_RE = re.compile(r"^(?P<base>[A-Z][A-Z0-9]*)X(?P<count>\d+)$")
_GNSM_RE = re.compile(r"^G(?P<g>\d+)S(?P<s>\d+)?$")
_GN_RE = re.compile(r"^G(?P<n>\d+)$")
_ANY_RE = re.compile(r"^A(?P<n>\d+)Y$")
_AN_RE = re.compile(r"^A(?P<n>\d+)$")

_MAX_REPEAT = 100  # safety cap so '(G4S)1000000' can't materialize a megasequence


def all_linker_names():
    """Canonical names + aliases. Use as argparse ``choices=`` value."""
    return sorted(set(LINKERS) | set(ALIASES))


def _check_repeat(count, what):
    if count < 1 or count > _MAX_REPEAT:
        raise ValueError(
            "Linker %s must be between 1 and %d (got %d)" % (
                what, _MAX_REPEAT, count))


def _make_repeat_of(base_name, count, canonical_name):
    """Build a synthetic Linker by repeating ``base_name`` ``count`` times.

    Refuses to repeat 2A linkers (codon-frozen, positional skipping
    mechanism — repeating wouldn't add cleavage events) or any other
    Linker carrying a blessed DNA sequence.
    """
    base = get_linker(base_name)
    if base.dna:
        raise ValueError(
            "Linker '%s' has a codon-frozen DNA sequence (e.g. 2A "
            "skipping is positional); repeating it would not produce "
            "additional functional events. Use the base linker once "
            "or pick a different family." % base_name)
    return Linker(
        name=canonical_name,
        amino_acids=base.amino_acids * count,
        freeze_in_mrna=False,
        inert_in_peptide_mode=base.inert_in_peptide_mode,
        citation=base.citation,
    )


def get_linker(name):
    """Resolve a linker by name; supports the compositional grammar in
    the module docstring.

    Resolution order: aliases → static LINKERS → ``(BASE)N`` /
    ``(BASE)xN`` → ``BASExN`` → ``GnSm`` (literal) → ``AnY`` → ``An``
    → ``Gn``.
    Names are uppercased before lookup so ``g4s2``, ``G4Sx2``, and
    ``g4sx2`` all resolve identically. ValueError on miss.
    """
    name = name.upper()
    canonical = ALIASES.get(name, name)
    if canonical in LINKERS:
        return LINKERS[canonical]

    # (BASE)N or (BASE)xN — explicit repeat with parens
    m = _PAREN_REPEAT_RE.match(canonical)
    if m:
        count = int(m.group("count"))
        _check_repeat(count, "repeat count")
        return _make_repeat_of(m.group("base"), count, canonical)

    # BASExN — explicit repeat with x separator
    m = _X_REPEAT_RE.match(canonical)
    if m:
        count = int(m.group("count"))
        _check_repeat(count, "repeat count")
        return _make_repeat_of(m.group("base"), count, canonical)

    # GnSm — n glycines + m serines, literal single unit
    m = _GNSM_RE.match(canonical)
    if m:
        g = int(m.group("g"))
        s = int(m.group("s")) if m.group("s") else 1
        _check_repeat(g, "glycine count")
        _check_repeat(s, "serine count")
        return Linker(
            name=canonical,
            amino_acids="G" * g + "S" * s,
            citation=_GS_CITATION,
        )

    # AnY — n alanines + tyrosine
    m = _ANY_RE.match(canonical)
    if m:
        n = int(m.group("n"))
        _check_repeat(n, "alanine count")
        # Use AAY's citation for the canonical n=2 case; for longer
        # variants, the synthesized citation explains the extrapolation.
        citation = _AAY_CITATION if n == 2 else _ANY_CITATION
        return Linker(
            name=canonical,
            amino_acids="A" * n + "Y",
            citation=citation,
        )

    # An — polyalanine, no Y. AAA is the empirically-tested form
    # (Aguilar-Gurrieri 2023); longer variants (AAAA, AAAAA) are
    # length extrapolations.
    m = _AN_RE.match(canonical)
    if m:
        n = int(m.group("n"))
        _check_repeat(n, "alanine count")
        return Linker(
            name=canonical,
            amino_acids="A" * n,
            citation=_ALANINE_CITATION,
        )

    # Gn — pure polyglycine, no serine. Biophysically the most
    # flexible linker family but unstudied as an inter-epitope spacer
    # in vaccines (see _POLYG_CITATION).
    m = _GN_RE.match(canonical)
    if m:
        n = int(m.group("n"))
        _check_repeat(n, "glycine count")
        return Linker(
            name=canonical,
            amino_acids="G" * n,
            citation=_POLYG_CITATION,
        )

    raise ValueError(
        "Unknown linker '%s'. Available named entries: %s. "
        "Compositional forms: (BASE)N or BASExN for repeats, GnSm for "
        "literal n-glycines + m-serines, AnY for n-alanines + Y." % (
            name, ', '.join(all_linker_names())))


def iter_named_antigens(ranked_vaccine_peptides, candidates_per_slot=1):
    """Yield ``(name, fragment, vaccine_peptide)`` per ranked candidate.

    Both peptide and mRNA assembly walk the same ranked list and need
    the same per-candidate name (``gene_chr_pos_ref_alt`` with an
    ``_alt<k>`` suffix on alternates). Centralizing the loop here keeps
    that naming consistent across modalities.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(varcode.Variant, list[VaccinePeptide])]
    candidates_per_slot : int
        How many ranked alternates to walk per variant. The first
        candidate gets the bare name; subsequent candidates get
        ``_alt1``, ``_alt2``, etc.
    """
    for variant, peptides in ranked_vaccine_peptides:
        if not peptides:
            continue
        for idx, peptide in enumerate(peptides[:max(1, candidates_per_slot)]):
            fragment = peptide.mutant_protein_fragment
            base_name = "%s_%s_%s_%s_%s" % (
                getattr(fragment, 'gene_name', None) or 'unknown',
                variant.contig,
                variant.start,
                variant.ref or '.',
                variant.alt or '.',
            )
            if idx > 0:
                base_name = "%s_alt%d" % (base_name, idx)
            yield base_name, fragment, peptide


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


def target_epitopes_sorted(vaccine_peptide):
    """All mutant ``CandidateEpitope`` records on a VaccinePeptide, score-sorted
    (best first).

    Returns the full list (possibly empty). The pipeline already sorts
    ``target_epitopes`` by score, so this is a thin wrapper that
    exposes the slice point to callers.

    Used by both peptide and mRNA vaccine assembly when
    ``antigen_content='minimal_epitope'`` — the antigen is one or more
    short MHC ligands instead of a mutation-spanning long window. The
    caller chooses how many to take (typically via
    ``epitopes_per_antigen``); designs that pack multiple top ligands
    from the same variant are first-class, not an afterthought.
    """
    return list(
        getattr(vaccine_peptide, 'target_epitopes', None) or [])


_STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def truncate_at_stop_codon(amino_acids):
    """Truncate an amino-acid string at the first ``*`` (stop codon).

    Translation halts at the stop, so any AA after a ``*`` doesn't
    exist in the cell. Stop-loss / readthrough variants are the
    rare exception where downstream residues become real, but for
    those the upstream loader / Isovar should emit a fragment that
    stops at the new C-terminus rather than embedding the original
    stop in the middle of the sequence. So: truncate at the first
    ``*`` unconditionally.

    Returns the substring up to (but not including) the first ``*``.
    Strings without ``*`` are returned unchanged.
    """
    if not amino_acids:
        return amino_acids
    idx = amino_acids.find('*')
    if idx < 0:
        return amino_acids
    return amino_acids[:idx]


def has_only_standard_amino_acids(amino_acids):
    """True iff every character of ``amino_acids`` is one of the 20
    standard amino acid letters. Used as a sanity gate before passing
    sequences to manufacturability / hydropathy / codon-optimization
    code that key off a 20-letter alphabet.
    """
    if not amino_acids:
        return True
    return all(aa in _STANDARD_AMINO_ACIDS for aa in amino_acids)


def top_target_epitopes(vaccine_peptide, n=1):
    """Top ``n`` score-sorted mutant epitopes on a VaccinePeptide.

    Returns a list of length ``min(n, len(predictions))`` — possibly
    empty when no predictions exist. ``n=1`` (the default) gives the
    legacy "single top ligand" semantics; pass a larger ``n`` to get
    multiple top ligands from the same variant.
    """
    if n <= 0:
        return []
    return target_epitopes_sorted(vaccine_peptide)[:n]
