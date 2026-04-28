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
``max_antigens_per_construct``, antigens spill into additional
constructs returned in the same order.
"""

import json
import logging
from dataclasses import dataclass, field
from typing import Optional

from dnachisel import (
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceTranslation,
)
from dnachisel.biotools import reverse_translate

from .mrna_library import (
    LINKERS,
    MITDS,
    SIGNAL_PEPTIDES,
    UTRS_3P,
    UTRS_5P,
)

logger = logging.getLogger(__name__)


STOP_CODON = "TAA"


@dataclass
class MRNAOptions:
    """User-configurable mRNA construct parameters."""
    signal_peptide: Optional[str] = "tPA"
    linker: str = "GS3"
    include_mitd: bool = False
    mitd: str = "HLA_A"
    utr_5p: str = "HBB"
    utr_3p: str = "HBB"
    codon_species: str = "h_sapiens"
    codon_method: str = "use_best_codon"
    max_antigens_per_construct: int = 10
    max_length_nt: int = 4000
    avoid_patterns: tuple = ()


@dataclass
class MRNAConstruct:
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
                   avoid_patterns=()):
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

    Returns
    -------
    str
        Optimized DNA sequence (same protein when translated).
    """
    if not amino_acids:
        return ""
    initial = reverse_translate(amino_acids)
    constraints = [EnforceTranslation()]
    for pattern in avoid_patterns:
        constraints.append(AvoidPattern(pattern))
    problem = DnaOptimizationProblem(
        sequence=initial,
        constraints=constraints,
        objectives=[CodonOptimize(species=species, method=method)],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence


def _antigen_aa_sequences(ranked_vaccine_peptides):
    """Yield ``(name, amino_acid_string)`` for each ranked vaccine peptide."""
    for variant, peptides in ranked_vaccine_peptides:
        if not peptides:
            continue
        peptide = peptides[0]
        fragment = peptide.mutant_protein_fragment
        name = "%s_%s_%s_%s_%s" % (
            fragment.gene_name or 'unknown',
            variant.contig,
            variant.start,
            variant.ref or '.',
            variant.alt or '.',
        )
        yield name, fragment.amino_acids


def _build_protein_string(antigen_aas, signal_peptide_aa, linker_aa,
                          mitd_aa):
    """Concatenate signal peptide + antigens + linker/MITD into one protein.

    Ensures the result begins with M so the back-translated CDS has a
    start codon. When a signal peptide is supplied, its N-terminal M
    serves as the start; otherwise an M is prepended in front of the
    first antigen if it doesn't already start with one.
    """
    parts = []
    if signal_peptide_aa:
        parts.append(signal_peptide_aa)
    body = linker_aa.join(antigen_aas)
    if not signal_peptide_aa and not body.startswith("M"):
        body = "M" + body
    parts.append(body)
    if mitd_aa:
        parts.append(linker_aa)
        parts.append(mitd_aa)
    return "".join(parts)


def _pack_constructs(antigen_pairs, options, signal_peptide_aa, linker_aa,
                     mitd_aa, utr_5p_dna, utr_3p_dna):
    """Greedy bin-packing of antigens into constructs honoring the caps."""
    # When there is no signal peptide, the assembler prepends an ATG to the
    # CDS body if the first antigen doesn't already start with M (see
    # _build_protein_string). Reserve 3 nt up-front to keep packing honest.
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
            or len(current) >= options.max_antigens_per_construct
        )
        if cap_hit and current:
            constructs.append(current)
            current = []
            current_aa_nt = 0
            antigen_nt = len(aa) * 3
            linker_nt = 0
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
    if current:
        constructs.append(current)
    return constructs


def assemble_mrna_constructs(ranked_vaccine_peptides, options=None):
    """Assemble mRNA constructs from ranked vaccine peptides.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(varcode.Variant, list[VaccinePeptide])]
    options : MRNAOptions or None

    Returns
    -------
    list[MRNAConstruct]
    """
    options = options or MRNAOptions()
    signal_peptide_aa = (
        _resolve_named(SIGNAL_PEPTIDES, options.signal_peptide,
                       'signal peptide')
        if options.signal_peptide else ""
    )
    linker_aa = _resolve_named(LINKERS, options.linker, 'linker')
    mitd_aa = _resolve_named(MITDS, options.mitd, 'MITD') if options.include_mitd else ""
    utr_5p_dna = _resolve_named(UTRS_5P, options.utr_5p, "5' UTR")
    utr_3p_dna = _resolve_named(UTRS_3P, options.utr_3p, "3' UTR")

    antigens = list(_antigen_aa_sequences(ranked_vaccine_peptides))
    if not antigens:
        return []

    packed = _pack_constructs(
        antigens, options, signal_peptide_aa, linker_aa, mitd_aa,
        utr_5p_dna, utr_3p_dna)

    constructs = []
    for i, antigen_group in enumerate(packed):
        names = [n for n, _ in antigen_group]
        antigen_aas = [aa for _, aa in antigen_group]
        protein = _build_protein_string(
            antigen_aas, signal_peptide_aa, linker_aa, mitd_aa)
        coding_dna = codon_optimize(
            protein,
            species=options.codon_species,
            method=options.codon_method,
            avoid_patterns=options.avoid_patterns,
        )
        sequence = utr_5p_dna + coding_dna + STOP_CODON + utr_3p_dna
        constructs.append(MRNAConstruct(
            name="construct_%02d" % (i + 1),
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
