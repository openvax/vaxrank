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

"""Assemble mRNA vaccine constructs from ranked vaccine peptides.

A construct is laid out as:

    5' UTR ─ [signal peptide]? ─ antigen1 ─ linker ─ ... ─ antigenN
           ─ [linker ─ MITD]? ─ stop ─ 3' UTR

The CDS must begin with an ATG. Canonical signal peptides (HLA-A,
tPA, IgK) all start with methionine, so when a signal peptide is
selected the start codon is the M at its N-terminus. When the user
opts out of a signal peptide and the first antigen does not start
with M, the assembler prepends one to keep the CDS translatable.

Each amino-acid block is back-translated and codon-optimized for the
target species via DnaChisel (see ``codon_optimize``). When the
combined antigen string would exceed ``max_length_nt`` or
``antigens_per_construct``, antigens spill into additional
constructs returned in the same order.
"""

import json
import logging
from dataclasses import dataclass, field
from typing import Optional

from dnachisel import (
    AvoidChanges,
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceTranslation,
    Location,
)
from dnachisel.biotools import reverse_translate

from .mrna_library import (
    MITDS,
    SIGNAL_PEPTIDES,
    UTRS_3P,
    UTRS_5P,
)
from .vaccine_library import get_linker, select_antigen_window

logger = logging.getLogger(__name__)


STOP_CODON = "TAA"


# Friendly → DnaChisel canonical mapping for codon-table species. Keys are
# lowercased-stripped aliases; values are the species keys DnaChisel passes
# to python_codon_tables.
_SPECIES_ALIASES = {
    # Human
    "human": "h_sapiens",
    "homo sapiens": "h_sapiens",
    "h. sapiens": "h_sapiens",
    "h.sapiens": "h_sapiens",
    "h_sapiens": "h_sapiens",
    "9606": "h_sapiens",
    # Mouse
    "mouse": "m_musculus",
    "murine": "m_musculus",
    "mus musculus": "m_musculus",
    "m. musculus": "m_musculus",
    "m.musculus": "m_musculus",
    "m_musculus": "m_musculus",
    "10090": "m_musculus",
    # Fly
    "fly": "d_melanogaster",
    "drosophila": "d_melanogaster",
    "drosophila melanogaster": "d_melanogaster",
    "d. melanogaster": "d_melanogaster",
    "d_melanogaster": "d_melanogaster",
    # Yeast
    "yeast": "s_cerevisiae",
    "s. cerevisiae": "s_cerevisiae",
    "saccharomyces cerevisiae": "s_cerevisiae",
    "s_cerevisiae": "s_cerevisiae",
    # Worm
    "worm": "c_elegans",
    "c. elegans": "c_elegans",
    "caenorhabditis elegans": "c_elegans",
    "c_elegans": "c_elegans",
    # E. coli
    "e. coli": "e_coli",
    "e.coli": "e_coli",
    "escherichia coli": "e_coli",
    "e_coli": "e_coli",
    # Chicken
    "chicken": "g_gallus",
    "g. gallus": "g_gallus",
    "gallus gallus": "g_gallus",
    "g_gallus": "g_gallus",
    # B. subtilis
    "b. subtilis": "b_subtilis",
    "bacillus subtilis": "b_subtilis",
    "b_subtilis": "b_subtilis",
}


def normalize_codon_species(name):
    """Map a user-supplied species name to DnaChisel's canonical form.

    Accepts common names (human, mouse), genus-species ("homo sapiens",
    "mus musculus"), abbreviations ("h. sapiens", "m. musculus"), the
    underscore-prefixed DnaChisel form ("h_sapiens"), and a few NCBI
    taxids. Unknown names pass through unchanged (DnaChisel will raise
    if it can't resolve them).
    """
    if not name:
        return name
    key = name.strip().lower()
    return _SPECIES_ALIASES.get(key, name)


@dataclass
class RNAConstructConfig:
    """User-configurable mRNA construct parameters.

    Defaults reflect a typical neoantigen mRNA design: 5 antigen
    windows of 15-20 aa each in a single construct, with a tPA signal
    peptide and HBB UTRs. ``max_length_nt`` is a belt-and-suspenders
    cap that triggers spillover into additional constructs only on
    pathologically long inputs.
    """
    signal_peptide: Optional[str] = "tPA"
    linker: str = "G4S3"
    include_mitd: bool = False
    mitd: str = "HLA_A"
    utr_5p: str = "HBB"
    utr_3p: str = "HBB"
    codon_species: str = "h_sapiens"
    codon_method: str = "use_best_codon"
    min_antigen_length_aa: int = 15
    max_antigen_length_aa: int = 20
    antigens_per_construct: int = 5
    max_constructs: int = 1
    candidates_per_slot: int = 1
    max_length_nt: int = 4000
    avoid_patterns: tuple = ()


@dataclass
class RNAConstruct:
    """A single assembled mRNA construct."""
    name: str
    antigen_names: list
    sequence: str
    components: dict = field(default_factory=dict)


def _resolve_named(table, name, kind):
    if name not in table:
        raise ValueError(
            "Unknown %s '%s'. Available: %s" % (
                kind, name, ', '.join(sorted(table))))
    return table[name]


def codon_optimize(amino_acids, species="h_sapiens", method="use_best_codon",
                   avoid_patterns=(), frozen_segments=()):
    """Back-translate and codon-optimize an amino-acid sequence to DNA.

    Parameters
    ----------
    amino_acids : str
    species : str
        DnaChisel species key (e.g. 'h_sapiens'); see python_codon_tables.
    method : str
        DnaChisel CodonOptimize method ('use_best_codon', 'match_codon_usage',
        or 'harmonize_rca').
    avoid_patterns : iterable of str
        Patterns / restriction sites to avoid (passed to DnaChisel
        AvoidPattern constraints).
    frozen_segments : iterable of (start_aa, end_aa, blessed_dna)
        Amino-acid positions [start_aa, end_aa) whose DNA should match
        ``blessed_dna`` exactly and not be touched by codon optimization.
        Used to preserve published codon usage for 2A self-cleaving
        peptides (see vaccine_library.LINKER_P2A et al.).

    Returns
    -------
    str
        Optimized DNA sequence (same protein when translated).
    """
    if not amino_acids:
        return ""
    species = normalize_codon_species(species)
    frozen_segments = sorted(frozen_segments)
    # Build the initial DNA by substituting blessed nt segments where
    # provided, otherwise back-translating freely. Track nt ranges to
    # later pin via AvoidChanges.
    initial_parts = []
    frozen_nt_ranges = []
    cursor = 0
    nt_offset = 0
    for start_aa, end_aa, blessed_dna in frozen_segments:
        if start_aa < cursor:
            raise ValueError("frozen_segments must not overlap")
        if (end_aa - start_aa) * 3 != len(blessed_dna):
            raise ValueError(
                "blessed DNA length %d does not match AA window %d:%d"
                % (len(blessed_dna), start_aa, end_aa))
        prefix_aa = amino_acids[cursor:start_aa]
        if prefix_aa:
            prefix_dna = reverse_translate(prefix_aa)
            initial_parts.append(prefix_dna)
            nt_offset += len(prefix_dna)
        initial_parts.append(blessed_dna)
        frozen_nt_ranges.append((nt_offset, nt_offset + len(blessed_dna)))
        nt_offset += len(blessed_dna)
        cursor = end_aa
    suffix_aa = amino_acids[cursor:]
    if suffix_aa:
        initial_parts.append(reverse_translate(suffix_aa))
    initial = "".join(initial_parts)

    constraints = [EnforceTranslation()]
    for pattern in avoid_patterns:
        constraints.append(AvoidPattern(pattern))
    for start_nt, end_nt in frozen_nt_ranges:
        constraints.append(AvoidChanges(location=Location(start_nt, end_nt)))
    problem = DnaOptimizationProblem(
        sequence=initial,
        constraints=constraints,
        objectives=[CodonOptimize(species=species, method=method)],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence


def _antigen_aa_sequences(ranked_vaccine_peptides, max_antigen_length_aa,
                          min_antigen_length_aa=0, candidates_per_slot=1):
    """Yield ``(name, amino_acid_string)`` per antigen, mutation-centered
    and clipped to ``max_antigen_length_aa``.

    Warns when the emitted window falls below ``min_antigen_length_aa`` —
    typically a sign that the underlying fragment is shorter than the
    user's configured floor (e.g. a stop-loss extension that produced
    only a few mutant residues).

    ``candidates_per_slot`` walks through additional ranked candidates
    per variant; alternates get a ``_alt<k>`` suffix on their name so
    different constructs can be traced back.
    """
    for variant, peptides in ranked_vaccine_peptides:
        if not peptides:
            continue
        for idx, peptide in enumerate(peptides[:max(1, candidates_per_slot)]):
            fragment = peptide.mutant_protein_fragment
            name = "%s_%s_%s_%s_%s" % (
                fragment.gene_name or 'unknown',
                variant.contig,
                variant.start,
                variant.ref or '.',
                variant.alt or '.',
            )
            if idx > 0:
                name = "%s_alt%d" % (name, idx)
            window = select_antigen_window(
                fragment, name, max_antigen_length_aa)
            if len(window) < min_antigen_length_aa:
                logger.warning(
                    "Antigen %s emitted at %d aa, below "
                    "--mrna-min-antigen-length-aa (%d).",
                    name, len(window), min_antigen_length_aa)
            yield name, window


def _build_protein_with_segments(antigen_aas, signal_peptide_aa, linker,
                                 mitd_aa):
    """Concatenate signal peptide + antigens + linker/MITD into one protein.

    Returns ``(protein_str, frozen_segments)`` where ``frozen_segments``
    is a list of ``(start_aa, end_aa, blessed_dna)`` tuples for linker
    occurrences that should be codon-frozen during mRNA optimization
    (i.e. ``Linker.freeze_in_mrna`` and ``Linker.dna`` set).

    Ensures the result begins with M so the back-translated CDS has a
    start codon. When a signal peptide is supplied, its N-terminal M
    serves as the start; otherwise an M is prepended in front of the
    first antigen if it doesn't already start with one.
    """
    linker_aa = linker.amino_acids
    freeze = bool(linker.freeze_in_mrna and linker.dna)
    parts = []
    frozen = []
    aa_offset = 0
    if signal_peptide_aa:
        parts.append(signal_peptide_aa)
        aa_offset += len(signal_peptide_aa)
    if not signal_peptide_aa and antigen_aas and not antigen_aas[0].startswith("M"):
        parts.append("M")
        aa_offset += 1
    for i, aa in enumerate(antigen_aas):
        if i > 0:
            if freeze:
                frozen.append((aa_offset, aa_offset + len(linker_aa), linker.dna))
            parts.append(linker_aa)
            aa_offset += len(linker_aa)
        parts.append(aa)
        aa_offset += len(aa)
    if mitd_aa:
        if freeze:
            frozen.append((aa_offset, aa_offset + len(linker_aa), linker.dna))
        parts.append(linker_aa)
        aa_offset += len(linker_aa)
        parts.append(mitd_aa)
        aa_offset += len(mitd_aa)
    return "".join(parts), frozen


def _pack_constructs(antigen_pairs, options, signal_peptide_aa, linker,
                     mitd_aa, utr_5p_dna, utr_3p_dna):
    """Greedy bin-packing of antigens into constructs honoring the caps."""
    linker_aa = linker.amino_acids
    # When there is no signal peptide, the assembler prepends an ATG to the
    # CDS body if the first antigen doesn't already start with M (see
    # _build_protein_with_segments). Reserve 3 nt up-front to keep packing
    # honest.
    start_codon_nt = 3 if not signal_peptide_aa else 0
    fixed_overhead_nt = (
        len(utr_5p_dna)
        + (len(signal_peptide_aa) * 3 if signal_peptide_aa else 0)
        + start_codon_nt
        + ((len(linker_aa) + len(mitd_aa)) * 3 if mitd_aa else 0)
        + len(STOP_CODON)
        + len(utr_3p_dna)
    )
    constructs = []
    current = []
    current_aa_nt = 0
    for name, aa in antigen_pairs:
        antigen_nt = len(aa) * 3
        linker_nt = len(linker_aa) * 3 if current else 0
        projected = fixed_overhead_nt + current_aa_nt + linker_nt + antigen_nt
        cap_hit = (
            projected > options.max_length_nt
            or len(current) >= options.antigens_per_construct
        )
        if cap_hit and current:
            constructs.append(current)
            current = []
            current_aa_nt = 0
            antigen_nt = len(aa) * 3
            linker_nt = 0
            if len(constructs) >= options.max_constructs:
                logger.warning(
                    "Reached --mrna-max-constructs (%d); dropping "
                    "remaining antigens including %s.",
                    options.max_constructs, name)
                return constructs
        # A single antigen larger than the length cap is still emitted in
        # its own construct — splitting one antigen would corrupt the
        # epitope. Warn the user so they can adjust caps or drop it.
        if not current and (fixed_overhead_nt + antigen_nt) > options.max_length_nt:
            logger.warning(
                "Antigen %s exceeds --mrna-max-length-nt (%d > %d) on its "
                "own; emitting in a single-antigen construct anyway.",
                name, fixed_overhead_nt + antigen_nt, options.max_length_nt)
        current.append((name, aa))
        current_aa_nt += linker_nt + antigen_nt
    if current and len(constructs) < options.max_constructs:
        constructs.append(current)
    elif current:
        logger.warning(
            "Dropped final %d antigen(s); --mrna-max-constructs (%d) reached.",
            len(current), options.max_constructs)
    return constructs


def assemble_mrna_constructs(ranked_vaccine_peptides, options=None):
    """Assemble mRNA constructs from ranked vaccine peptides.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(varcode.Variant, list[VaccinePeptide])]
    options : RNAConstructConfig or None

    Returns
    -------
    list[RNAConstruct]
    """
    options = options or RNAConstructConfig()
    signal_peptide_aa = (
        _resolve_named(SIGNAL_PEPTIDES, options.signal_peptide,
                       'signal peptide')
        if options.signal_peptide else ""
    )
    # Natural secretory signal peptides start with the initial M (it's the
    # start codon's product). All bundled SIGNAL_PEPTIDES satisfy this; guard
    # against future additions that don't, since otherwise the CDS would lack
    # an ATG and the ribosome wouldn't initiate.
    if signal_peptide_aa and not signal_peptide_aa.startswith("M"):
        raise ValueError(
            "Signal peptide '%s' does not start with M; the resulting CDS "
            "would have no start codon." % options.signal_peptide)
    linker = get_linker(options.linker)
    mitd_aa = _resolve_named(MITDS, options.mitd, 'MITD') if options.include_mitd else ""
    utr_5p_dna = _resolve_named(UTRS_5P, options.utr_5p, "5' UTR")
    utr_3p_dna = _resolve_named(UTRS_3P, options.utr_3p, "3' UTR")

    antigens = list(_antigen_aa_sequences(
        ranked_vaccine_peptides,
        max_antigen_length_aa=options.max_antigen_length_aa,
        min_antigen_length_aa=options.min_antigen_length_aa,
        candidates_per_slot=options.candidates_per_slot))
    if not antigens:
        return []

    packed = _pack_constructs(
        antigens, options, signal_peptide_aa, linker, mitd_aa,
        utr_5p_dna, utr_3p_dna)

    constructs = []
    for i, antigen_group in enumerate(packed):
        names = [n for n, _ in antigen_group]
        antigen_aas = [aa for _, aa in antigen_group]
        protein, frozen_segments = _build_protein_with_segments(
            antigen_aas, signal_peptide_aa, linker, mitd_aa)
        coding_dna = codon_optimize(
            protein,
            species=options.codon_species,
            method=options.codon_method,
            avoid_patterns=options.avoid_patterns,
            frozen_segments=frozen_segments,
        )
        sequence = utr_5p_dna + coding_dna + STOP_CODON + utr_3p_dna
        constructs.append(RNAConstruct(
            name="seq_%03d" % (i + 1),
            antigen_names=names,
            sequence=sequence,
            components={
                'utr_5p': options.utr_5p,
                'signal_peptide': options.signal_peptide,
                'linker': options.linker,
                'mitd': options.mitd if options.include_mitd else None,
                'utr_3p': options.utr_3p,
                'codon_species': options.codon_species,
                'codon_method': options.codon_method,
                'protein': protein,
            },
        ))
    return constructs


def write_mrna_outputs(constructs, fasta_path, manifest_path=None):
    """Write FASTA + optional JSON manifest describing each construct.

    Manifest entries follow a shared schema (``modality``, ``name``,
    ``length``, ``length_unit``, ``antigen_names``, ``components``,
    ``manufacturability``) so the peptide-mode writer can produce a
    structurally compatible manifest.
    """
    with open(fasta_path, 'w') as f:
        for c in constructs:
            f.write(">%s antigens=%s length=%d\n" % (
                c.name, ','.join(c.antigen_names), len(c.sequence)))
            for i in range(0, len(c.sequence), 80):
                f.write(c.sequence[i:i + 80] + "\n")

    if manifest_path:
        manifest = [
            {
                'modality': 'mrna',
                'name': c.name,
                'length': len(c.sequence),
                'length_unit': 'nt',
                'antigen_names': c.antigen_names,
                'components': c.components,
                'manufacturability': {},  # nt-level metrics: see openvax/vaxrank#245
            }
            for c in constructs
        ]
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)
