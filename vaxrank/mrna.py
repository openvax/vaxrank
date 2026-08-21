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

import csv
import json
import logging
import os
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
from .vaccine_library import (
    gene_names_from_antigen_names,
    get_linker,
    iter_named_antigens,
    select_antigen_window,
    top_target_epitopes,
)

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

    Defaults reproduce the BioNTech IVAC MUTANOME / FixVac design as
    published in Sahin et al., Nature 547:222, 2017 (doi:10.1038/nature23003)
    and originally engineered in Kreiter et al., J Immunol 180:309, 2008
    (doi:10.4049/jimmunol.180.1.309):

    - HLA-B signal peptide + HLA-B MITD (consistent self-pairing)
    - 5 antigens per construct, 25-aa windows
    - (G4S)2 inter-antigen linker
    - HBB 5' UTR + tandem 2× HBB 3' UTR (FI element / "2hBg")
    - Single construct per vaccine

    ``max_length_nt`` is a belt-and-suspenders cap that triggers
    spillover into additional constructs only on pathologically long
    inputs.
    """
    signal_peptide: Optional[str] = "HLA_B"
    linker: str = "(G4S)2"
    include_mitd: bool = True
    mitd: str = "HLA_B"
    utr_5p: str = "HBB"
    utr_3p: str = "HBB_FI"  # tandem 2× HBB
    codon_species: str = "h_sapiens"
    codon_method: str = "use_best_codon"
    # Antigen-design axes — shared semantics with PeptideConstructConfig.
    # ``antigen_content`` ∈ {'mutation_spanning', 'minimal_epitope'};
    # the latter packs top MHC ligands as antigens instead of
    # mutation-centered windows. ``epitopes_per_antigen`` only matters
    # for ``minimal_epitope`` content (legacy 1-per-VP semantics by
    # default; >1 packs multiple top ligands from the same variant).
    antigen_content: str = "mutation_spanning"
    epitopes_per_antigen: int = 1
    min_antigen_length_aa: int = 15
    max_antigen_length_aa: int = 25  # Sahin 2017 used 27mers; 25 is typical
    antigens_per_construct: int = 5
    # BioNTech FixVac canonical: 2× pentatope = 10 antigens at
    # 5 antigens/construct (Sahin 2017). Set to 1 for a single
    # construct; raise for broader coverage.
    max_constructs: int = 2
    candidates_per_slot: int = 1
    max_length_nt: int = 4000
    avoid_patterns: tuple = ()
    # Per-junction linker optimization (issue #247). On by default.
    # When enabled, the ``linker`` field above is the *fallback* if no
    # candidate outperforms it; otherwise each junction picks its
    # own linker from ``junction_swap_candidates`` to minimize
    # predicted MHC presentation of chimeric k-mers spanning the
    # junction. If the caller doesn't pass an mhc_predictor + alleles
    # to ``assemble_mrna_constructs``, the optimizer is skipped with
    # a warning and the shared linker is used at every junction.
    optimize_linkers: bool = True
    junction_swap_candidates: tuple = ()  # empty → use library default
    junction_kmer_lengths: tuple = (8, 9, 10, 11)
    junction_rank_strong: float = 0.5
    junction_rank_mild: float = 2.0
    # PolyA tail (issue #252). Default A120 matches BioNTech IVAC
    # MUTANOME / FixVac (Sahin 2017). When ``poly_a_segmented`` is
    # True, the tail is split as A_first + linker + A_rest, mirroring
    # BNT162b2 (A30 + GCATATGACT + A70 per Xia 2021, PMC8310186).
    poly_a_length: int = 120
    poly_a_segmented: bool = False
    poly_a_first_segment: int = 30
    poly_a_segment_linker: str = "GCATATGACT"


@dataclass
class RNAConstruct:
    """A single assembled mRNA construct.

    ``sequence`` is the full nt sequence including polyA tail (kept as
    the canonical "the construct" string for back-compat callers).
    The structured per-element manifest is exposed via ``elements``,
    ``antigens``, and the explicit ``cds_nt`` / ``no_polya_nt`` /
    ``full_nt`` views.
    """
    name: str
    antigen_names: list
    sequence: str
    components: dict = field(default_factory=dict)
    cds_aa: str = ""
    cds_nt: str = ""        # protein-coding DNA *including* stop codon
    no_polya_nt: str = ""   # 5' UTR + cds_nt + 3' UTR
    full_nt: str = ""       # no_polya_nt + polyA tail
    poly_a_nt: str = ""
    elements: dict = field(default_factory=dict)
    antigens: list = field(default_factory=list)


def _normalize_lookup_key(name):
    """Strip case + dash/underscore differences for table lookups.
    Lets users write ``HLA-B`` / ``hla_b`` / ``HLAB`` / ``HLA_B`` for
    a key the library spells ``HLA_B``."""
    return str(name).replace('-', '').replace('_', '').lower()


def _resolve_named(table, name, kind):
    """Resolve ``name`` against ``table`` with case + separator
    tolerance. Common biological identifiers (HLA-B, tPA, IgK) are
    written with mixed case and dashes in the literature; the table
    keys use a consistent ``HLA_B`` / ``tPA`` form, but lookups
    accept any sensible spelling."""
    if name in table:
        return table[name]
    needle = _normalize_lookup_key(name)
    for k, v in table.items():
        if _normalize_lookup_key(k) == needle:
            return v
    raise ValueError(
        "Unknown %s '%s'. Available: %s" % (
            kind, name, ', '.join(sorted(table))))


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


def summarize_mrna_ranking_decisions(
        ranked_variants_with_vaccine_peptides, options):
    """Produce a "which antigens land in the mRNA vaccine?" summary
    derived from the ranked-variants list and the
    :class:`RNAConstructConfig` caps. Used by report writers (#270)
    to surface the mRNA ranking-decision section of the ASCII / HTML
    / PDF report next to the existing per-variant tables.

    The mRNA assembler iterates the ranked variants in score order
    and packs antigens into constructs until
    ``antigens_per_construct × max_constructs`` is reached — so the
    first ``cap`` variants in ``ranked_variants_with_vaccine_peptides``
    are the ones that land in the vaccine, and the rest fall past the
    cap. (Edge case: a single antigen too long to fit is still
    emitted, which can shift the cut by one — the writer logs the
    edge-case warning separately.)

    Parameters
    ----------
    ranked_variants_with_vaccine_peptides : list[tuple]
        ``[(source, [VaccinePeptide])]`` — vaxrank's canonical
        ranked-output shape.
    options : RNAConstructConfig
        Construct-assembly knobs; the cap comes from
        ``options.antigens_per_construct * options.max_constructs``.

    Returns
    -------
    dict
        Shape:

        ``{
            'antigens_per_construct': int,
            'max_constructs': int,
            'cap': int,        # antigens_per_construct × max_constructs
            'total_ranked': int,
            'selected': [{'rank': 1, 'description': 'GENE_X_...',
                          'gene_name': 'GENE',
                          'combined_score': float, ...}, ...],
            'dropped': [...],   # past-cap, same shape
          }``

        ``selected`` and ``dropped`` are sorted by rank (ascending
        index in the ranked list).
    """
    cap = options.antigens_per_construct * options.max_constructs
    from .coverage import summarize_construction_decisions
    construction_summary = summarize_construction_decisions(
        ranked_variants_with_vaccine_peptides,
        cap=cap,
        target_alleles=[],
    )

    return {
        'antigens_per_construct': options.antigens_per_construct,
        'max_constructs': options.max_constructs,
        'cap': cap,
        'total_ranked': len(ranked_variants_with_vaccine_peptides),
        'selected': [
            {
                key: row[key]
                for key in ("rank", "gene_name", "description", "combined_score")
            }
            for row in construction_summary['selected']
        ],
        'dropped': [
            {
                key: row[key]
                for key in ("rank", "gene_name", "description", "combined_score")
            }
            for row in construction_summary['dropped']
        ],
    }


def _antigen_aa_sequences(ranked_vaccine_peptides, max_antigen_length_aa,
                          min_antigen_length_aa=0, candidates_per_slot=1,
                          antigen_content="mutation_spanning",
                          epitopes_per_antigen=1):
    """Yield ``(name, amino_acid_string)`` per antigen.

    Dispatches on ``antigen_content`` (shared semantics with
    ``peptide._antigen_records``):

    - ``'mutation_spanning'``: legacy configuration name for a window that
      preserves the antigen's targetable span, up to
      ``max_antigen_length_aa`` per VaccinePeptide (canonical mRNA antigen —
      BioNTech FixVac-style 25mer for mutation antigens).
    - ``'minimal_epitope'``: emit the top ``epitopes_per_antigen``
      MHC ligands per VaccinePeptide as separate (short) antigens.
      Concatenated minimal-epitope mRNA constructs (the
      "string-of-beads" design — Whitton et al., J Virol 67:348,
      1993, PMID 7677954) become first-class.

    Warns when a mutation_spanning window falls below
    ``min_antigen_length_aa`` (typically a stop-loss / frameshift
    fragment shorter than the floor); for minimal_epitope content,
    skip-with-info when no epitope predictions are available.

    Naming + alt-suffix logic comes from ``iter_named_antigens`` so
    the antigen names match the peptide-modality output.
    """
    if antigen_content not in ('mutation_spanning', 'minimal_epitope'):
        raise ValueError(
            "antigen_content must be 'mutation_spanning' or "
            "'minimal_epitope'; got %r" % antigen_content)
    for name, antigen, peptide in iter_named_antigens(
            ranked_vaccine_peptides, candidates_per_slot=candidates_per_slot):
        if antigen_content == "minimal_epitope":
            tops = top_target_epitopes(peptide, n=epitopes_per_antigen)
            if not tops:
                logger.info(
                    "Skipping %s in minimal_epitope mode: no mutant "
                    "epitope predictions available.", name)
                continue
            for k, ep in enumerate(tops):
                suffix = "_epitope" if len(tops) == 1 else "_epitope%d" % (k + 1)
                yield name + suffix, ep.sequence
        else:
            window = select_antigen_window(
                antigen, name, max_antigen_length_aa)
            if len(window) < min_antigen_length_aa:
                logger.warning(
                    "Antigen %s emitted at %d aa, below "
                    "--mrna-min-antigen-length-aa (%d).",
                    name, len(window), min_antigen_length_aa)
            yield name, window


def _build_protein_with_segments(antigen_aas, antigen_names, signal_peptide_aa,
                                 signal_peptide_name, linker, mitd_aa,
                                 mitd_name, per_junction_linkers=None,
                                 pre_mitd_linker=None):
    """Concatenate signal peptide + antigens + linker/MITD into one protein.

    Returns ``(protein_str, frozen_segments, aa_segments)`` where:
      - ``frozen_segments`` is a list of ``(start_aa, end_aa, blessed_dna)``
        for linkers that should be codon-frozen during mRNA optimization
        (``Linker.freeze_in_mrna`` and ``Linker.dna`` set).
      - ``aa_segments`` is a list of ``dict``s, one per CDS element, with
        keys ``kind`` ('signal_peptide' / 'antigen' / 'linker' / 'mitd' /
        'start_codon'), ``name``, ``start_aa``, ``end_aa``, ``aa``. The
        caller slices the back-translated DNA at ``[start*3:end*3]`` to
        recover per-element nt.

    Ensures the result begins with M so the back-translated CDS has a
    start codon. When a signal peptide is supplied, its N-terminal M
    serves as the start; otherwise an M is prepended in front of the
    first antigen if it doesn't already start with one.

    When ``per_junction_linkers`` is provided (length must equal
    ``len(antigen_aas) - 1``), each inter-antigen junction uses its
    own Linker instead of the shared ``linker``. ``pre_mitd_linker``
    similarly overrides the shared linker for the last-antigen ↔ MITD
    junction. Both come from the per-junction linker optimizer (#247).
    """
    if per_junction_linkers is not None and antigen_aas:
        expected = max(0, len(antigen_aas) - 1)
        if len(per_junction_linkers) != expected:
            raise ValueError(
                "per_junction_linkers length %d != antigens-1 (%d)" % (
                    len(per_junction_linkers), expected))
    if antigen_names is None:
        antigen_names = ["antigen_%d" % i for i in range(len(antigen_aas))]

    parts = []
    frozen = []
    segments = []
    aa_offset = 0

    def _emit(kind, name, aa, junction_index=None):
        nonlocal aa_offset
        seg = {
            'kind': kind,
            'name': name,
            'aa': aa,
            'start_aa': aa_offset,
            'end_aa': aa_offset + len(aa),
        }
        if junction_index is not None:
            seg['junction_index'] = junction_index
        segments.append(seg)
        parts.append(aa)
        aa_offset += len(aa)

    if signal_peptide_aa:
        _emit('signal_peptide', signal_peptide_name, signal_peptide_aa)
    if not signal_peptide_aa and antigen_aas and not antigen_aas[0].startswith("M"):
        _emit('start_codon', 'start_codon', "M")
    for i, aa in enumerate(antigen_aas):
        if i > 0:
            if per_junction_linkers is not None:
                this_linker = per_junction_linkers[i - 1]
            else:
                this_linker = linker
            this_aa = this_linker.amino_acids
            if this_linker.freeze_in_mrna and this_linker.dna:
                frozen.append(
                    (aa_offset, aa_offset + len(this_aa), this_linker.dna))
            _emit('linker', this_linker.name, this_aa, junction_index=i - 1)
        _emit('antigen', antigen_names[i], aa)
    if mitd_aa:
        # Pre-MITD linker: use override if supplied (per-junction
        # optimizer path), else fall back to the shared `linker`.
        mitd_link = pre_mitd_linker if pre_mitd_linker is not None else linker
        linker_aa = mitd_link.amino_acids
        if mitd_link.freeze_in_mrna and mitd_link.dna:
            frozen.append(
                (aa_offset, aa_offset + len(linker_aa), mitd_link.dna))
        # Tag the pre-MITD linker with junction_index = "mitd" sentinel so
        # the manifest can distinguish it from inter-antigen junctions.
        _emit('linker', mitd_link.name, linker_aa, junction_index='mitd')
        _emit('mitd', mitd_name, mitd_aa)
    return "".join(parts), frozen, segments


def build_poly_a(options):
    """Build the polyA tail nt string from RNAConstructConfig.

    Non-segmented: ``A`` * poly_a_length (matches Sahin 2017 BNT122).
    Segmented: ``A`` * first_segment + segment_linker + ``A`` * rest,
    where rest = max(0, poly_a_length - poly_a_first_segment). The
    segment linker's bases are NOT counted toward poly_a_length —
    the configured length is the total adenosine count, mirroring
    Xia 2021's description of BNT162b2 (A30 + GCAUAUGACU + A70 =
    100 A's split by 10-nt linker).
    """
    if options.poly_a_length < 0:
        raise ValueError(
            "poly_a_length must be >= 0 (got %d)" % options.poly_a_length)
    if not options.poly_a_segmented:
        return "A" * options.poly_a_length
    first = max(0, min(options.poly_a_first_segment, options.poly_a_length))
    rest = options.poly_a_length - first
    return "A" * first + options.poly_a_segment_linker + "A" * rest


def _packing_linker_aa(options, linker):
    """Linker AA string the packer should bill against, sized to the
    longest possible per-junction substitution.

    When ``optimize_linkers`` is on, the per-junction optimizer can
    swap in any candidate from ``junction_swap_candidates`` (or the
    library default ``JUNCTION_SWAP_CANDIDATES`` when the option tuple
    is empty). Use the longest of {shared linker, candidates} so the
    bin-packer is conservative against later substitution.
    """
    shared = linker.amino_acids
    if not options.optimize_linkers:
        return shared
    from .junction_swap import JUNCTION_SWAP_CANDIDATES
    cand_names = (
        tuple(options.junction_swap_candidates) or JUNCTION_SWAP_CANDIDATES
    )
    longest = shared
    for n in cand_names:
        cand = get_linker(n).amino_acids
        if len(cand) > len(longest):
            longest = cand
    return longest


def _pack_constructs(antigen_pairs, options, signal_peptide_aa, linker,
                     mitd_aa, utr_5p_dna, utr_3p_dna):
    """Greedy bin-packing of antigens into constructs honoring the caps.

    Packing uses ``max(shared_linker, longest_candidate)`` for size
    estimation when the per-junction optimizer (#247) is on, so the
    swap step at assembly time can never exceed the cap by picking a
    candidate longer than the shared linker.
    """
    linker_aa = _packing_linker_aa(options, linker)
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
    total_antigens = (
        len(antigen_pairs) if hasattr(antigen_pairs, '__len__') else None)
    for i, (name, aa) in enumerate(antigen_pairs):
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
                # Top-k selection log. Ranking N candidates and
                # keeping the top k (the chosen 5 antigens × 1
                # construct here) is *the whole point* — not a
                # warning, not "spill" / "drop" panic language.
                # Operators who want a different cut control it via
                # ``--mrna-max-constructs`` /
                # ``--mrna-antigens-per-construct``. ``total_antigens``
                # is None on a generator caller; the count line is
                # the only thing that requires it.
                kept = sum(len(c) for c in constructs)
                if total_antigens is not None:
                    logger.info(
                        "mRNA assembly: selected top %d / %d "
                        "ranked antigen(s) into %d construct(s) at "
                        "--mrna-max-constructs=%d. Lower-ranked "
                        "antigens past the cap (e.g. %s) are "
                        "available in --output-csv / "
                        "--output-json-file for downstream review.",
                        kept, total_antigens, len(constructs),
                        options.max_constructs, name)
                else:
                    logger.info(
                        "mRNA assembly: selected top %d ranked "
                        "antigen(s) into %d construct(s) at "
                        "--mrna-max-constructs=%d (next-best: %s).",
                        kept, len(constructs),
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


def assemble_mrna_constructs(ranked_vaccine_peptides, options=None,
                             mhc_predictor=None, mhc_alleles=None,
                             reference_proteome=None, *,
                             target_alleles=None):
    """Assemble mRNA constructs from ranked vaccine peptides.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(source, list[VaccinePeptide])]
    options : RNAConstructConfig or None
    mhc_predictor : optional, mhctools.BasePredictor
        Required when ``options.optimize_linkers`` is True.
        Used by the per-junction linker optimizer to score chimeric k-mers.
    mhc_alleles : optional, list[str]
        Patient HLA allele names. Required when
        ``optimize_linkers`` is True.
    reference_proteome : optional
        Container that answers ``kmer in reference_proteome``. Junction
        k-mers found in the reference proteome are filtered out before
        scoring (already-tolerated, not new presentation).
    target_alleles : optional, list[str]
        Patient MHC alleles for coverage-aware antigen selection.
        When non-empty, the selector reorders the ranked list to
        maximize coverage before bin-packing into constructs — see
        :func:`vaxrank.coverage.select_antigens_for_coverage`. When
        empty / None, today's pure-score order is preserved.
        Typically the same list as ``mhc_alleles`` (the linker
        optimizer's input); kept separate so the two can diverge if
        a future caller wants to optimize against a population
        cohort while ranking against the patient.

    Returns
    -------
    list[RNAConstruct]
    """
    options = options or RNAConstructConfig()
    if target_alleles:
        from .coverage import select_antigens_for_coverage
        cap = options.antigens_per_construct * options.max_constructs
        head = select_antigens_for_coverage(
            ranked_vaccine_peptides, target_alleles, cap)
        head_keys = {id(item) for item in head}
        tail = [item for item in ranked_vaccine_peptides
                if id(item) not in head_keys]
        ranked_vaccine_peptides = head + tail
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
        candidates_per_slot=options.candidates_per_slot,
        antigen_content=options.antigen_content,
        epitopes_per_antigen=options.epitopes_per_antigen))
    if not antigens:
        return []

    packed = _pack_constructs(
        antigens, options, signal_peptide_aa, linker, mitd_aa,
        utr_5p_dna, utr_3p_dna)

    constructs = []
    for i, antigen_group in enumerate(packed):
        names = [n for n, _ in antigen_group]
        antigen_aas = [aa for _, aa in antigen_group]

        per_junction_linkers = None
        pre_mitd_linker = None  # overrides shared linker for the MITD junction
        junction_swap_meta = None
        n_junctions = max(0, len(antigen_aas) - 1)
        if options.optimize_linkers and (n_junctions > 0 or mitd_aa):
            if mhc_predictor is None or not mhc_alleles:
                logger.warning(
                    "optimize_linkers=True but no mhc_predictor / "
                    "mhc_alleles supplied; falling back to the shared "
                    "linker at every junction. Pass mhc_predictor + "
                    "mhc_alleles to assemble_mrna_constructs (or run "
                    "via the CLI which auto-wires them) to enable the "
                    "optimizer.")
                junction_swap_meta = {
                    'enabled': False,
                    'note': "no predictor / alleles available",
                }
            else:
                from .junction_swap import (
                    JUNCTION_SWAP_CANDIDATES,
                    optimize_linkers,
                )
                cand_names = (
                    tuple(options.junction_swap_candidates)
                    or JUNCTION_SWAP_CANDIDATES
                )
                swap = optimize_linkers(
                    antigen_aas=antigen_aas,
                    alleles=mhc_alleles,
                    predictor=mhc_predictor,
                    candidate_names=cand_names,
                    k_lengths=tuple(options.junction_kmer_lengths),
                    reference_proteome=reference_proteome,
                    rank_strong=options.junction_rank_strong,
                    rank_mild=options.junction_rank_mild,
                    mitd_aa=mitd_aa or None,
                    default_linker_name=options.linker,
                )
                # `swap.default_*_burden` is always set when
                # default_linker_name is supplied (which it always is on
                # this code path); we asserted above that the optimizer
                # ran with a default. Read directly.
                default_strong = swap.default_strong_burden
                default_mild = swap.default_burden
                # Canonical name for the shared/default linker so the
                # manifest's `chosen` field is uniform across branches.
                default_name = linker.name
                n_total_junctions = n_junctions + (1 if mitd_aa else 0)
                if ((swap.strong_burden, swap.burden)
                        < (default_strong, default_mild)):
                    chosen_list = swap.chosen_linker_per_junction
                    # Last entry corresponds to MITD junction when mitd_aa
                    # was provided; split it off.
                    if mitd_aa and len(chosen_list) == n_junctions + 1:
                        per_junction_linkers = chosen_list[:n_junctions]
                        pre_mitd_linker = chosen_list[-1]
                    else:
                        per_junction_linkers = chosen_list
                    junction_swap_meta = {
                        'enabled': True,
                        'chosen': swap.linker_names(),
                        'burden_strong': swap.strong_burden,
                        'burden_mild': swap.burden,
                        'default_burden_strong': default_strong,
                        'default_burden_mild': default_mild,
                    }
                    logger.info(
                        "Construct %d: per-junction linker optimizer "
                        "reduced strong burden %d → %d, mild burden %d → %d",
                        i + 1, default_strong, swap.strong_burden,
                        default_mild, swap.burden)
                else:
                    junction_swap_meta = {
                        'enabled': True,
                        'chosen': [default_name] * n_total_junctions,
                        'burden_strong': default_strong,
                        'burden_mild': default_mild,
                        'note': "default linker beat or tied all candidates",
                    }

        protein, frozen_segments, aa_segments = _build_protein_with_segments(
            antigen_aas, names, signal_peptide_aa, options.signal_peptide,
            linker, mitd_aa,
            options.mitd if options.include_mitd else None,
            per_junction_linkers=per_junction_linkers,
            pre_mitd_linker=pre_mitd_linker)
        coding_dna = codon_optimize(
            protein,
            species=options.codon_species,
            method=options.codon_method,
            avoid_patterns=options.avoid_patterns,
            frozen_segments=frozen_segments,
        )
        # Three sequence views: cds_nt = back-translated protein + STOP;
        # no_polya_nt = 5' UTR + cds_nt + 3' UTR; full_nt = no_polya_nt
        # + polyA tail. ``sequence`` keeps full_nt for back-compat with
        # callers that read it as the canonical "the construct" string.
        cds_nt = coding_dna + STOP_CODON
        no_polya_nt = utr_5p_dna + cds_nt + utr_3p_dna
        poly_a_nt = build_poly_a(options)
        full_nt = no_polya_nt + poly_a_nt

        # Per-element manifest with both AA and nt forms. Antigens are
        # listed separately (one entry per packed antigen). Linkers per
        # junction are also listed separately so the manifest preserves
        # per-junction-swap information.
        antigens_meta = []
        linker_segments = []
        for seg in aa_segments:
            nt_slice = coding_dna[seg['start_aa'] * 3:seg['end_aa'] * 3]
            seg_record = {
                'name': seg['name'],
                'aa': seg['aa'],
                'nt': nt_slice,
                'length_aa': len(seg['aa']),
                'length_nt': len(nt_slice),
                'start_aa': seg['start_aa'],
                'end_aa': seg['end_aa'],
            }
            if seg['kind'] == 'antigen':
                antigens_meta.append(seg_record)
            elif seg['kind'] == 'linker':
                ls = dict(seg_record)
                ls['junction_index'] = seg.get('junction_index')
                linker_segments.append(ls)

        # Helper: build the per-element record for a given segment kind.
        def _seg_record(kind, declared_name):
            seg = next((s for s in aa_segments if s['kind'] == kind), None)
            if seg is None:
                return None
            nt = coding_dna[seg['start_aa'] * 3:seg['end_aa'] * 3]
            return {
                'name': declared_name,
                'aa': seg['aa'],
                'nt': nt,
                'length_aa': len(seg['aa']),
                'length_nt': len(nt),
            }

        elements = {}
        elements['utr_5p'] = {
            'name': options.utr_5p, 'nt': utr_5p_dna,
            'length_nt': len(utr_5p_dna),
        }
        elements['signal_peptide'] = (
            _seg_record('signal_peptide', options.signal_peptide)
            if signal_peptide_aa else None
        )
        elements['linkers_per_junction'] = linker_segments
        elements['mitd'] = (
            _seg_record('mitd', options.mitd) if mitd_aa else None
        )
        elements['stop_codon'] = STOP_CODON
        elements['utr_3p'] = {
            'name': options.utr_3p, 'nt': utr_3p_dna,
            'length_nt': len(utr_3p_dna),
        }
        elements['poly_a'] = {
            'length_nt': options.poly_a_length,
            'segmented': options.poly_a_segmented,
            'segment_linker': (
                options.poly_a_segment_linker
                if options.poly_a_segmented else None),
            'first_segment_nt': (
                options.poly_a_first_segment
                if options.poly_a_segmented else None),
            'nt': poly_a_nt,
        }
        elements['codon_species'] = options.codon_species
        elements['codon_method'] = options.codon_method
        if junction_swap_meta is not None:
            elements['junction_swap'] = junction_swap_meta

        # Back-compat ``components`` dict: same fields the 2.12.x
        # manifest carried, so existing readers don't break. The
        # richer per-element view lives on RNAConstruct.elements.
        components = {
            'utr_5p': options.utr_5p,
            'signal_peptide': options.signal_peptide,
            'linker': options.linker,
            'mitd': options.mitd if options.include_mitd else None,
            'utr_3p': options.utr_3p,
            'codon_species': options.codon_species,
            'codon_method': options.codon_method,
            'protein': protein,
        }
        if junction_swap_meta is not None:
            components['junction_swap'] = junction_swap_meta

        constructs.append(RNAConstruct(
            # Modality-stamped namespace (was bare ``seq_NNN``):
            # users mix peptide + mRNA outputs in the same downstream
            # pipeline and a bare ``seq_001`` collides with peptide
            # construct names. Mirrors the peptide writer's
            # ``peptide_NNN`` scheme.
            name="mrna_%03d" % (i + 1),
            antigen_names=names,
            sequence=full_nt,
            components=components,
            cds_aa=protein,
            cds_nt=cds_nt,
            no_polya_nt=no_polya_nt,
            full_nt=full_nt,
            poly_a_nt=poly_a_nt,
            elements=elements,
            antigens=antigens_meta,
        ))
    return constructs


def _wrap_fasta(seq, line_width=80):
    return "\n".join(seq[i:i + line_width] for i in range(0, len(seq), line_width))


_FASTA_LIKE_SUFFIXES = (".fasta", ".fa", ".fna", ".ffn", ".fas")


def _validate_output_dir(output_dir):
    """Reject paths that look like the pre-2.14 single-FASTA target.

    The old API took a FASTA file path here; the new API takes a
    directory. Silently creating ``out.fasta/`` when the user meant a
    file is a sharp footgun, so block it loudly.
    """
    if os.path.isfile(output_dir):
        raise ValueError(
            "--output-dir for --vaccine-type=mrna is a *directory* (writes cds.fasta / "
            "no_polyA.fasta / full.fasta into it); got an existing file "
            "%r. Pass a directory path instead." % output_dir)
    if any(output_dir.lower().endswith(s) for s in _FASTA_LIKE_SUFFIXES):
        raise ValueError(
            "--output-dir for --vaccine-type=mrna is a *directory* (writes cds.fasta / "
            "no_polyA.fasta / full.fasta into it); got %r which looks "
            "like a FASTA file path. Pass a directory path instead "
            "(e.g. drop the .fasta suffix)." % output_dir)


def write_mrna_outputs(constructs, output_dir, manifest_path=None,
                       csv_path=None, csv_include_full_rows=True):
    """Write three FASTAs (cds / full / no_polyA) + optional manifest + CSV.

    Output layout (issue #252):
      <output_dir>/cds.fasta        — coding DNA (start codon → stop codon)
      <output_dir>/no_polyA.fasta   — 5' UTR + CDS + 3' UTR
      <output_dir>/full.fasta       — no_polyA + polyA tail

    All three contain one record per construct, matched by ``name``.
    When ``poly_a_length=0``, full.fasta and no_polyA.fasta are
    identical by construction; an info-level log line notes this.

    ``manifest_path`` (optional): JSON manifest with the structured
    per-element view (``elements``, ``antigens``, all sequence variants).
    ``csv_path`` (optional): long-format CSV — one row per (construct,
    element) — exposing AA + nt for every layer for spreadsheet
    inspection.
    ``csv_include_full_rows``: emit summary rows for cds / no_polyA /
    full (each carrying the full-length nt). Default True. Set False
    when the per-element rows are all you want — keeps cells narrow
    enough for spreadsheet column views.

    Manifest schema notes
    ---------------------
    Some ``elements`` keys are **None when the corresponding component
    is disabled**: ``signal_peptide`` is None when no signal peptide is
    selected, ``mitd`` is None when ``include_mitd=False``. Downstream
    readers must null-check both. ``elements['linkers_per_junction']``
    is always a list (possibly empty for a single-antigen construct).
    """
    _validate_output_dir(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    fasta_specs = [
        ('cds.fasta', 'cds_nt'),
        ('no_polyA.fasta', 'no_polya_nt'),
        ('full.fasta', 'full_nt'),
    ]
    for filename, attr in fasta_specs:
        path = os.path.join(output_dir, filename)
        with open(path, 'w') as f:
            for c in constructs:
                seq = getattr(c, attr)
                # Use ';' as the antigen separator: GenBank
                # convention for compound description fields, and
                # avoids the trap where downstream tools split the
                # FASTA description on ',' (treating each antigen
                # name as a separate token).
                genes = gene_names_from_antigen_names(c.antigen_names)
                f.write(">%s antigens=%s length=%d genes=%s\n" % (
                    c.name, ';'.join(c.antigen_names), len(seq),
                    ';'.join(genes)))
                f.write(_wrap_fasta(seq))
                f.write("\n")

    if constructs and constructs[0].full_nt == constructs[0].no_polya_nt:
        logger.info(
            "poly_a_length=0: full.fasta and no_polyA.fasta are identical "
            "(no polyA tail appended).")

    if manifest_path:
        manifest = [
            {
                'modality': 'mrna',
                'name': c.name,
                'length': len(c.full_nt),  # canonical "the construct" length
                'length_unit': 'nt',
                'antigen_names': c.antigen_names,
                'lengths': {
                    'cds_aa': len(c.cds_aa),
                    'cds_nt': len(c.cds_nt),
                    'no_polya_nt': len(c.no_polya_nt),
                    'full_nt': len(c.full_nt),
                    'poly_a_nt': len(c.poly_a_nt),
                },
                'cds': {'aa': c.cds_aa, 'nt': c.cds_nt},
                'no_polya_nt': c.no_polya_nt,
                'full_nt': c.full_nt,
                'antigens': c.antigens,
                'elements': c.elements,
                'components': c.components,  # legacy 2.12 schema for back-compat
                'manufacturability': {},
            }
            for c in constructs
        ]
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)

    if csv_path:
        # Long format: one row per (construct, element). The `index`
        # column is integer-only (junction position for linkers,
        # antigen position for antigens, blank otherwise). Free-form
        # qualifiers (e.g. 'mitd', 'pre_mitd', 'segmented') go in the
        # `index_label` / `note` columns so spreadsheet sorts on
        # `index` stay numeric.
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow([
                'construct', 'element_kind', 'index', 'index_label', 'name',
                'aa', 'nt', 'length_aa', 'length_nt', 'note',
            ])
            for c in constructs:
                el = c.elements
                # 5' UTR
                u5 = el.get('utr_5p') or {}
                w.writerow([c.name, 'utr_5p', '', '', u5.get('name', ''),
                            '', u5.get('nt', ''), '', u5.get('length_nt', ''),
                            ''])
                # signal peptide
                sp = el.get('signal_peptide')
                if sp:
                    w.writerow([c.name, 'signal_peptide', '', '', sp['name'],
                                sp['aa'], sp['nt'],
                                sp['length_aa'], sp['length_nt'], ''])
                # antigens
                for j, a in enumerate(c.antigens):
                    w.writerow([c.name, 'antigen', j, '', a['name'],
                                a['aa'], a['nt'],
                                a['length_aa'], a['length_nt'], ''])
                # linkers per junction; the pre-MITD linker uses an
                # empty integer index + 'mitd' label so spreadsheet
                # sorts on 'index' stay numeric.
                for ls in el.get('linkers_per_junction', []):
                    ji = ls.get('junction_index')
                    if ji == 'mitd':
                        idx_int = ''
                        idx_label = 'mitd'
                        note = 'pre_mitd'
                    else:
                        idx_int = ji if isinstance(ji, int) else ''
                        idx_label = ''
                        note = ''
                    w.writerow([c.name, 'linker', idx_int, idx_label,
                                ls['name'], ls['aa'], ls['nt'],
                                ls['length_aa'], ls['length_nt'], note])
                # MITD
                m = el.get('mitd')
                if m:
                    w.writerow([c.name, 'mitd', '', '', m['name'],
                                m['aa'], m['nt'],
                                m['length_aa'], m['length_nt'], ''])
                # stop codon
                w.writerow([c.name, 'stop_codon', '', '', '', '',
                            el.get('stop_codon', ''), '', 3, ''])
                # 3' UTR
                u3 = el.get('utr_3p') or {}
                w.writerow([c.name, 'utr_3p', '', '', u3.get('name', ''),
                            '', u3.get('nt', ''), '', u3.get('length_nt', ''),
                            ''])
                # polyA
                pa = el.get('poly_a') or {}
                pa_note = ('segmented' if pa.get('segmented') else 'unsegmented')
                w.writerow([c.name, 'poly_a', '', '', '', '',
                            pa.get('nt', ''), '', pa.get('length_nt', ''),
                            pa_note])
                # full assembled views — opt-out via csv_include_full_rows.
                if csv_include_full_rows:
                    w.writerow([c.name, 'cds', '', '', '', c.cds_aa, c.cds_nt,
                                len(c.cds_aa), len(c.cds_nt), ''])
                    w.writerow([c.name, 'no_polyA', '', '', '', '',
                                c.no_polya_nt,
                                '', len(c.no_polya_nt), ''])
                    w.writerow([c.name, 'full', '', '', '', '', c.full_nt,
                                '', len(c.full_nt), ''])
