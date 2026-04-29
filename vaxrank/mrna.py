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
from .vaccine_library import (
    get_linker,
    iter_named_antigens,
    select_antigen_window,
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
    min_antigen_length_aa: int = 15
    max_antigen_length_aa: int = 25  # Sahin 2017 used 27mers; 25 is typical
    antigens_per_construct: int = 5
    max_constructs: int = 1
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
    optimize_junction_linkers: bool = True
    junction_swap_candidates: tuple = ()  # empty → use library default
    junction_kmer_lengths: tuple = (8, 9, 10, 11)
    junction_rank_strong: float = 0.5
    junction_rank_mild: float = 2.0


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

    Naming + alt-suffix logic comes from ``iter_named_antigens`` so the
    antigen names match the peptide-modality output.
    """
    for name, fragment, _peptide in iter_named_antigens(
            ranked_vaccine_peptides, candidates_per_slot=candidates_per_slot):
        window = select_antigen_window(fragment, name, max_antigen_length_aa)
        if len(window) < min_antigen_length_aa:
            logger.warning(
                "Antigen %s emitted at %d aa, below "
                "--mrna-min-antigen-length-aa (%d).",
                name, len(window), min_antigen_length_aa)
        yield name, window


def _build_protein_with_segments(antigen_aas, signal_peptide_aa, linker,
                                 mitd_aa, per_junction_linkers=None,
                                 pre_mitd_linker=None):
    """Concatenate signal peptide + antigens + linker/MITD into one protein.

    Returns ``(protein_str, frozen_segments)`` where ``frozen_segments``
    is a list of ``(start_aa, end_aa, blessed_dna)`` tuples for linker
    occurrences that should be codon-frozen during mRNA optimization
    (i.e. ``Linker.freeze_in_mrna`` and ``Linker.dna`` set).

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
            if per_junction_linkers is not None:
                this_linker = per_junction_linkers[i - 1]
            else:
                this_linker = linker
            this_aa = this_linker.amino_acids
            if this_linker.freeze_in_mrna and this_linker.dna:
                frozen.append((aa_offset, aa_offset + len(this_aa), this_linker.dna))
            parts.append(this_aa)
            aa_offset += len(this_aa)
        parts.append(aa)
        aa_offset += len(aa)
    if mitd_aa:
        # Pre-MITD linker: use override if supplied (per-junction
        # optimizer path), else fall back to the shared `linker`.
        mitd_link = pre_mitd_linker if pre_mitd_linker is not None else linker
        linker_aa = mitd_link.amino_acids
        if mitd_link.freeze_in_mrna and mitd_link.dna:
            frozen.append((aa_offset, aa_offset + len(linker_aa), mitd_link.dna))
        parts.append(linker_aa)
        aa_offset += len(linker_aa)
        parts.append(mitd_aa)
        aa_offset += len(mitd_aa)
    return "".join(parts), frozen


def _pack_constructs(antigen_pairs, options, signal_peptide_aa, linker,
                     mitd_aa, utr_5p_dna, utr_3p_dna):
    """Greedy bin-packing of antigens into constructs honoring the caps.

    Packing uses the *shared* linker length to estimate construct size.
    The per-junction linker optimizer (#247) can pick different linkers
    per junction at assembly time; for the candidate set
    JUNCTION_SWAP_CANDIDATES (3-10 aa range) the per-junction swing is
    at most ~7 aa = 21 nt per junction relative to the (G4S)2 default,
    well within the headroom of a 4000-nt cap. Reach for a tighter
    estimate only if max_length_nt becomes binding.
    """
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


def assemble_mrna_constructs(ranked_vaccine_peptides, options=None,
                             mhc_predictor=None, mhc_alleles=None,
                             reference_proteome=None):
    """Assemble mRNA constructs from ranked vaccine peptides.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(varcode.Variant, list[VaccinePeptide])]
    options : RNAConstructConfig or None
    mhc_predictor : optional, mhctools.BasePredictor
        Required when ``options.optimize_junction_linkers`` is True.
        Used by the per-junction linker optimizer to score chimeric k-mers.
    mhc_alleles : optional, list[str]
        Patient HLA allele names. Required when
        ``optimize_junction_linkers`` is True.
    reference_proteome : optional
        Container that answers ``kmer in reference_proteome``. Junction
        k-mers found in the reference proteome are filtered out before
        scoring (already-tolerated, not new presentation).

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

        per_junction_linkers = None
        pre_mitd_linker = None  # overrides shared linker for the MITD junction
        junction_swap_meta = None
        n_junctions = max(0, len(antigen_aas) - 1)
        if options.optimize_junction_linkers and (n_junctions > 0 or mitd_aa):
            if mhc_predictor is None or not mhc_alleles:
                logger.warning(
                    "optimize_junction_linkers=True but no mhc_predictor / "
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
                    optimize_junction_linkers,
                )
                cand_names = (
                    tuple(options.junction_swap_candidates)
                    or JUNCTION_SWAP_CANDIDATES
                )
                swap = optimize_junction_linkers(
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
                # `swap.default_*_burden` is set when default_linker_name
                # is supplied; only swap when a candidate strictly beats
                # the default.
                default_strong = getattr(swap, 'default_strong_burden', None)
                default_mild = getattr(swap, 'default_burden', None)
                if (default_strong is not None
                        and (swap.strong_burden, swap.burden)
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
                    n_total_junctions = n_junctions + (1 if mitd_aa else 0)
                    junction_swap_meta = {
                        'enabled': True,
                        'chosen': [options.linker] * n_total_junctions,
                        'burden_strong': default_strong or 0,
                        'burden_mild': default_mild or 0,
                        'note': "default linker beat or tied all candidates",
                    }

        protein, frozen_segments = _build_protein_with_segments(
            antigen_aas, signal_peptide_aa, linker, mitd_aa,
            per_junction_linkers=per_junction_linkers,
            pre_mitd_linker=pre_mitd_linker)
        coding_dna = codon_optimize(
            protein,
            species=options.codon_species,
            method=options.codon_method,
            avoid_patterns=options.avoid_patterns,
            frozen_segments=frozen_segments,
        )
        sequence = utr_5p_dna + coding_dna + STOP_CODON + utr_3p_dna
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
            name="seq_%03d" % (i + 1),
            antigen_names=names,
            sequence=sequence,
            components=components,
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
