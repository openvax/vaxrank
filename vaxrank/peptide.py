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

"""Assemble peptide vaccine constructs from ranked vaccine peptides.

Three assembly modes:

- ``slp`` (default) — one synthetic long peptide per construct, using
  the vaccine peptide's amino-acid window straight from the ranking
  pipeline. This is the canonical neoantigen peptide vaccine format.
- ``minimal_epitope`` — emit only the predicted MHC ligand for each
  vaccine peptide (the top-scoring epitope), one peptide per
  construct. Useful for short-peptide pools.
- ``multi_epitope`` — concatenate antigens with a linker, mirroring
  the mRNA mode. Uncommon for peptide vaccines; included for
  vocabulary symmetry. 2A linkers are accepted but flagged as
  functionally inert (ribosomal skipping is co-translational).

The output API mirrors ``vaxrank/mrna.py`` — see issue #245 for the
symmetric-design rationale.
"""

import csv
import json
import logging
from dataclasses import dataclass, field

from .manufacturability import ManufacturabilityScores
from .vaccine_library import get_linker

logger = logging.getLogger(__name__)


@dataclass
class PeptideOptions:
    """User-configurable peptide construct parameters.

    Defaults match the canonical PGV-001 personalized peptide vaccine
    layout: ~20 synthetic long peptides per pool, one antigen per
    peptide, 15-25 aa per peptide.
    """
    mode: str = "slp"            # 'slp' | 'minimal_epitope' | 'multi_epitope'
    linker: str = "G4S3"         # only used in multi_epitope mode
    min_antigen_length_aa: int = 15
    max_antigen_length_aa: int = 25
    antigens_per_construct: int = 1
    max_constructs: int = 20
    candidates_per_slot: int = 1
    n_terminal_acetylation: bool = False
    c_terminal_amidation: bool = False


@dataclass
class PeptideConstruct:
    """A single assembled peptide construct."""
    name: str
    sequence: str
    antigen_names: list
    components: dict = field(default_factory=dict)
    manufacturability: dict = field(default_factory=dict)


def _antigen_records(ranked_vaccine_peptides, mode, max_antigen_length_aa,
                     candidates_per_slot=1):
    """Yield ``(name, amino_acids)`` per antigen.

    ``candidates_per_slot`` controls how many alternate windows are
    emitted per variant: 1 picks just the top-ranked vaccine peptide,
    higher values walk through additional candidates the ranking
    pipeline produced. Each alternate gets a ``_alt<k>`` suffix on its
    base name for traceability.
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
            if mode == "minimal_epitope":
                top = _top_mutant_epitope(peptide)
                if top is None:
                    logger.info(
                        "Skipping %s in minimal_epitope mode: no mutant "
                        "epitope predictions available.", base_name)
                    continue
                yield base_name + "_epitope", top.peptide_sequence
            else:
                # slp + multi_epitope: pick a mutation-centered window
                yield base_name, _slp_window(
                    fragment, base_name, max_antigen_length_aa)


def _top_mutant_epitope(vaccine_peptide):
    """Pick the highest-scoring mutant epitope; None if none available."""
    predictions = getattr(vaccine_peptide, 'mutant_epitope_predictions', None) or []
    if not predictions:
        return None
    # mutant_epitope_predictions is already sorted by score in the pipeline
    return predictions[0]


def _slp_window(fragment, base_name, max_length_aa):
    """Return a sub-window of the vaccine peptide that preserves the mutation.

    A naive ``amino_acids[:max_length_aa]`` truncation can drop the
    mutated residues entirely when the mutation is past the head of the
    fragment — see issue #245 review. Instead, center a window of size
    ``max_length_aa`` on the mutation. If the mutation itself is longer
    than the cap (rare; e.g. very long inframe insertions), emit the
    full fragment unchanged with a warning.
    """
    aa = fragment.amino_acids
    if len(aa) <= max_length_aa:
        return aa
    mut_start = getattr(fragment, 'mutant_amino_acid_start_offset', 0)
    mut_end = getattr(fragment, 'mutant_amino_acid_end_offset', len(aa))
    mut_len = mut_end - mut_start
    if mut_len > max_length_aa:
        logger.warning(
            "Mutation in %s spans %d aa, longer than --peptide-max-length-aa "
            "(%d); emitting full fragment without truncation.",
            base_name, mut_len, max_length_aa)
        return aa
    midpoint = (mut_start + mut_end) // 2
    half = max_length_aa // 2
    start = max(0, midpoint - half)
    end = min(len(aa), start + max_length_aa)
    start = max(0, end - max_length_aa)
    return aa[start:end]


def _manufacturability_for(sequence):
    """Recompute ManufacturabilityScores from the final emitted sequence."""
    scores = ManufacturabilityScores.from_amino_acids(sequence)
    return {f: getattr(scores, f) for f in ManufacturabilityScores._fields}


def _pack_multi_epitope(records, options, linker):
    """Bin-pack antigens into multi-epitope constructs honoring AA caps.

    The per-construct cap is ``antigens_per_construct * max_antigen_length_aa
    + (antigens_per_construct - 1) * len(linker)`` — generous enough that a
    typical 5-antigen / 25-aa-window construct fits without churn.
    """
    linker_aa = linker.amino_acids
    # Belt-and-suspenders length cap to catch pathologically long antigens.
    max_construct_aa = (
        options.antigens_per_construct * options.max_antigen_length_aa
        + max(0, options.antigens_per_construct - 1) * len(linker_aa)
    )
    constructs = []
    current = []
    current_aa = 0
    for name, aa in records:
        if len(constructs) >= options.max_constructs:
            logger.warning(
                "Reached --peptide-max-constructs (%d); dropping remaining "
                "antigens including %s.", options.max_constructs, name)
            break
        antigen_aa = len(aa)
        linker_extra = len(linker_aa) if current else 0
        projected = current_aa + linker_extra + antigen_aa
        cap_hit = (
            projected > max_construct_aa
            or len(current) >= options.antigens_per_construct
        )
        if cap_hit and current:
            constructs.append(current)
            current = []
            current_aa = 0
            linker_extra = 0
            if len(constructs) >= options.max_constructs:
                logger.warning(
                    "Reached --peptide-max-constructs (%d); dropping "
                    "remaining antigens including %s.",
                    options.max_constructs, name)
                break
        if not current and antigen_aa > max_construct_aa:
            logger.warning(
                "Antigen %s exceeds the multi_epitope construct cap "
                "(%d > %d aa) on its own; emitting as-is.",
                name, antigen_aa, max_construct_aa)
        current.append((name, aa))
        current_aa += linker_extra + antigen_aa
    if current and len(constructs) < options.max_constructs:
        constructs.append(current)
    return constructs


def assemble_peptide_constructs(ranked_vaccine_peptides, options=None):
    """Assemble peptide constructs from ranked vaccine peptides.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(varcode.Variant, list[VaccinePeptide])]
    options : PeptideOptions or None

    Returns
    -------
    list[PeptideConstruct]
    """
    options = options or PeptideOptions()
    if options.mode not in ("slp", "minimal_epitope", "multi_epitope"):
        raise ValueError(
            "Unknown peptide mode '%s'; expected one of slp, "
            "minimal_epitope, multi_epitope." % (options.mode,))

    records = list(_antigen_records(
        ranked_vaccine_peptides, options.mode,
        options.max_antigen_length_aa,
        candidates_per_slot=options.candidates_per_slot))
    if not records:
        return []

    base_components = {
        'mode': options.mode,
        'n_terminal_acetylation': options.n_terminal_acetylation,
        'c_terminal_amidation': options.c_terminal_amidation,
    }

    constructs = []
    if options.mode in ("slp", "minimal_epitope"):
        # One construct per (variant, candidate) record. Cap at max_constructs.
        for i, (name, sequence) in enumerate(records):
            if len(constructs) >= options.max_constructs:
                logger.info(
                    "Reached --peptide-max-constructs (%d); dropping "
                    "remaining %d candidate(s).",
                    options.max_constructs, len(records) - i)
                break
            if len(sequence) < options.min_antigen_length_aa:
                logger.warning(
                    "Antigen %s emitted at %d aa, below "
                    "--peptide-min-antigen-length-aa (%d).",
                    name, len(sequence), options.min_antigen_length_aa)
            constructs.append(PeptideConstruct(
                name="peptide_%03d" % (len(constructs) + 1),
                sequence=sequence,
                antigen_names=[name],
                components=dict(base_components),
                manufacturability=_manufacturability_for(sequence),
            ))
        return constructs

    # multi_epitope
    linker = get_linker(options.linker)
    if linker.inert_in_peptide_mode:
        logger.warning(
            "Linker '%s' is functionally inert in synthesized peptides "
            "(2A skipping is co-translational). Including for vocabulary "
            "symmetry; manifest annotates components.linker_inert=True.",
            options.linker)
    packed = _pack_multi_epitope(records, options, linker)
    for group in packed:
        names = [n for n, _ in group]
        sequence = linker.amino_acids.join(aa for _, aa in group)
        components = dict(base_components)
        components.update({
            'linker': options.linker,
            'linker_inert': linker.inert_in_peptide_mode,
        })
        constructs.append(PeptideConstruct(
            name="peptide_%03d" % (len(constructs) + 1),
            sequence=sequence,
            antigen_names=names,
            components=components,
            manufacturability=_manufacturability_for(sequence),
        ))
    return constructs


def _modification_label(acetyl, amide, sequence):
    """Render the vendor-display string with N-/C-terminal modifications.

    Returns ``None`` when neither modification is requested — the caller
    omits the column entirely so the order form doesn't carry a redundant
    duplicate of the bare sequence.
    """
    if not (acetyl or amide):
        return None
    return ("Ac-" if acetyl else "") + sequence + ("-NH2" if amide else "")


def write_peptide_outputs(constructs, fasta_path, manifest_path=None,
                          order_form_path=None, options=None):
    """Write FASTA, optional JSON manifest, and optional vendor order form.

    Manifest schema matches ``mrna.write_mrna_outputs`` so callers can
    consume both modalities uniformly.
    """
    options = options or PeptideOptions()
    has_mods = options.n_terminal_acetylation or options.c_terminal_amidation

    with open(fasta_path, 'w') as f:
        for c in constructs:
            f.write(">%s antigens=%s length=%d\n" % (
                c.name, ','.join(c.antigen_names), len(c.sequence)))
            for i in range(0, len(c.sequence), 80):
                f.write(c.sequence[i:i + 80] + "\n")

    if manifest_path:
        manifest = [
            {
                'modality': 'peptide',
                'name': c.name,
                'length': len(c.sequence),
                'length_unit': 'aa',
                'antigen_names': c.antigen_names,
                'components': c.components,
                'manufacturability': c.manufacturability,
            }
            for c in constructs
        ]
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)

    if order_form_path:
        with open(order_form_path, 'w', newline='') as f:
            writer = csv.writer(f)
            header = ['name', 'sequence', 'length',
                      'n_terminal_modification', 'c_terminal_modification']
            if has_mods:
                header.append('displayed_sequence')
            header += ['antigen_names', 'notes']
            writer.writerow(header)
            n_term = 'Acetyl' if options.n_terminal_acetylation else 'Free'
            c_term = 'Amide' if options.c_terminal_amidation else 'Free'
            for c in constructs:
                notes = ""
                if c.components.get('linker_inert'):
                    notes = (
                        "Construct contains a 2A linker; ribosomal "
                        "skipping is co-translational and does not "
                        "occur in synthesized peptides.")
                row = [c.name, c.sequence, len(c.sequence), n_term, c_term]
                if has_mods:
                    row.append(_modification_label(
                        options.n_terminal_acetylation,
                        options.c_terminal_amidation,
                        c.sequence))
                row += [';'.join(c.antigen_names), notes]
                writer.writerow(row)
