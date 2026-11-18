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
  the selected antigen window.
- ``minimal_epitope`` — emit only the predicted MHC ligand for each
  vaccine peptide (the top-scoring epitope), one peptide per
  construct. Useful for short-peptide pools.
- ``multi_epitope`` — concatenate antigens with a linker into one
  synthetic peptide molecule. This is distinct from a pool of separate
  peptides. 2A linkers are accepted but flagged as
  functionally inert (ribosomal skipping is co-translational).

The output API mirrors ``vaxrank/mrna.py`` — see issue #245 for the
symmetric-design rationale.
"""

import csv
import json
import logging
from dataclasses import dataclass, field

from .assembly_antigens import antigens_for_assembly, pack_antigens, product_gene_names, record_assembly_decisions
from .manufacturability import ManufacturabilityScores
from .construct_sequence import ConstructPlacement, validate_construct_placements
from .native_serialization import to_native_json
from .vaccine_library import (
    get_linker,
    iter_named_antigens,
)

logger = logging.getLogger(__name__)


@dataclass
class PeptideConstructConfig:
    """User-configurable peptide construct parameters.

    Two orthogonal axes drive the vaccine design (shared with
    ``RNAConstructConfig``):

    - ``antigen_content`` ∈ ``{'mutation_spanning', 'minimal_epitope'}``
      — what each antigen *is*. ``mutation_spanning`` extracts a
      mutation-centered window of up to ``max_antigen_length_aa``
      (the SLP/long-peptide content). ``minimal_epitope`` extracts
      the top ``epitopes_per_antigen`` MHC ligands per VaccinePeptide
      and emits each as its own short antigen.
    - ``antigens_per_construct`` — how many antigens to *concatenate*
      into one construct (with the configured ``linker``).

    The four corners of the matrix:

      content=mutation_spanning, per_construct=1   → SLP (one peptide
                                                     per ranked vaccine
                                                     peptide)
      content=mutation_spanning, per_construct=N   → Multi-SLP /
                                                     "multi-epitope long"
                                                     concatenation
      content=minimal_epitope,  per_construct=1   → Minimal-epitope
                                                     peptide (one short
                                                     ligand per variant)
      content=minimal_epitope,  per_construct=N   → Concatenated
                                                     minimal-ligand
                                                     peptide

    Defaults produce up to 20 separate synthetic long peptides, one antigen
    per molecule, with requested antigen windows of 15-25 aa. A pool contains
    separate peptide molecules; increasing ``antigens_per_construct`` instead
    joins antigen sequences into one molecule, including intervening linkers.

    Legacy ``mode`` field (``slp`` / ``minimal_epitope`` /
    ``multi_epitope``) is preserved as a back-compat shorthand; when
    set non-default it derives ``antigen_content`` and
    ``antigens_per_construct`` in ``__post_init__``.

    ``scale_mg`` / ``purity_percent`` / ``counterion`` drive the
    vendor order form's per-construct columns; per-vaccine constants
    today, kept here so the writer can render them without a separate
    config object.
    """
    mode: str = "slp"            # legacy alias; see __post_init__
    antigen_content: str = "mutation_spanning"  # 'mutation_spanning' | 'minimal_epitope'
    epitopes_per_antigen: int = 1               # for minimal_epitope content
    linker: str = "G4Sx3"        # only used when antigens_per_construct > 1
                                  # G4Sx3 = (G4S)3 = GGGGSGGGGSGGGGS (15 aa)
    min_antigen_length_aa: int = 15
    max_antigen_length_aa: int = 25
    antigens_per_construct: int = 1
    max_constructs: int = 20
    candidates_per_slot: int = 1
    antigen_allocation: dict | None = None
    n_terminal_acetylation: bool = False
    c_terminal_amidation: bool = False
    scale_mg: float = 5.0           # synthesis scale per peptide
    purity_percent: float = 95.0    # HPLC purity target
    counterion: str = "TFA"         # default salt form (TFA / acetate / HCl / free)
    window_selection: dict | None = None

    def __post_init__(self):
        from .antigen_allocation import AntigenAllocation
        AntigenAllocation.resolve(self.antigen_allocation)
        # Normalize ``mode`` to lowercase so YAML configs that say
        # ``SLP`` / ``Slp`` / ``slp`` all resolve identically.
        # Argparse already lowercases via ``type=str.lower`` for the
        # CLI; this catches the YAML-config path and any direct
        # construction.
        if isinstance(self.mode, str):
            self.mode = self.mode.lower()
        # Derive the orthogonal axes from the legacy ``mode`` shorthand
        # when the user passed it (or the default 'slp' was kept). The
        # rule: if the user already set ``antigen_content`` /
        # ``antigens_per_construct`` non-default, those win — mode is
        # only consulted when it would change the dataclass default.
        if self.mode == "slp":
            # default; nothing to derive (mutation_spanning + per_construct=1)
            pass
        elif self.mode == "minimal_epitope":
            self.antigen_content = "minimal_epitope"
            # antigens_per_construct stays at user-set value (default 1)
        elif self.mode == "multi_epitope":
            # mutation_spanning content (legacy semantics);
            # antigens_per_construct must already be > 1 for this to
            # do anything different from slp — left to the user.
            self.antigen_content = "mutation_spanning"
        else:
            raise ValueError(
                "PeptideConstructConfig.mode must be one of "
                "{'slp', 'minimal_epitope', 'multi_epitope'}; got %r" %
                self.mode)
        if self.antigen_content not in (
                'mutation_spanning', 'minimal_epitope'):
            raise ValueError(
                "antigen_content must be 'mutation_spanning' or "
                "'minimal_epitope'; got %r" % self.antigen_content)
        if self.window_selection is not None:
            from .window_selection import WindowSelection
            WindowSelection.resolve(self.window_selection)
            if self.antigen_content != 'mutation_spanning' or self.antigens_per_construct != 1:
                raise ValueError('Peptide window selection requires one SLP per construct')


@dataclass
class PeptideConstruct:
    """A single assembled peptide construct."""
    name: str
    sequence: str
    antigen_names: list
    components: dict = field(default_factory=dict)
    manufacturability: dict = field(default_factory=dict)
    construct_placements: tuple = field(default_factory=tuple)
    antigens: list = field(default_factory=list)

    @property
    def antigen_sources(self):
        """Separate source mappings for each emitted antigen segment."""
        return {a['name']: a['sources'] for a in self.antigens if a.get('sources')}

    def __post_init__(self):
        self.construct_placements = tuple(self.construct_placements)
        validate_construct_placements(self.construct_placements, self.sequence, "peptide")


def peptide_construct_from_sequence(record):
    """Export an explicitly supplied construct without designing or changing it.

    The record is provenance, not a claim that this sequence passed selection
    or a safety policy. Use the construct audit before interpreting predictions.
    """
    return PeptideConstruct(
        name=record.name, sequence=record.sequence,
        antigen_names=[record.native_antigen.display_gene_name
                       if record.native_antigen else record.name],
        manufacturability=_manufacturability_for(record.sequence),
        construct_placements=(ConstructPlacement(record, 0),),
    )


def _manufacturability_for(sequence):
    """Recompute ManufacturabilityScores from the final emitted sequence."""
    scores = ManufacturabilityScores.from_amino_acids(sequence)
    return {f: getattr(scores, f) for f in ManufacturabilityScores._fields}


def assemble_peptide_constructs(
        ranked_vaccine_peptides, options=None, *,
        target_alleles=None, window_audit=None, source_evidence=None, selection_audit=None):
    """Assemble peptide constructs from ranked vaccine peptides.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(source, list[VaccinePeptide])]
    options : PeptideConstructConfig or None
    target_alleles : optional, list[str]
        Patient MHC alleles. When non-empty, the selector reorders
        the ranked list to maximize allele coverage before packing
        — see :func:`vaxrank.coverage.select_antigens_for_coverage`.
        When empty / None, the input order is preserved (today's
        pure-score behavior). Animal-agnostic: HLA, mouse H-2, swine
        SLA all work the same way.
    source_evidence : optional, dict
        Retained source evidence keyed by producer evidence ID, as in
        ``EpitopeDataset.selection['source_evidence']``. Supplies recorded
        source-event identities for product provenance; never changes scores.
    selection_audit : optional, dict
        Filled with exact-sequence deduplication and packing outcomes, including
        eligible antigens excluded by the construct capacity limit.

    Returns
    -------
    list[PeptideConstruct]
    """
    options = options or PeptideConstructConfig()
    selection_audit = selection_audit if selection_audit is not None else {}
    if options.window_selection is not None:
        from .window_selection import optimize_peptide_windows
        ranked_vaccine_peptides = optimize_peptide_windows(ranked_vaccine_peptides, options, window_audit)
    from .assembly_antigens import rank_assembly_windows
    ranked_vaccine_peptides = rank_assembly_windows(ranked_vaccine_peptides)
    if target_alleles:
        from .coverage import select_antigens_for_coverage
        cap = options.antigens_per_construct * options.max_constructs
        # Reorder up to ``cap`` antigens for coverage; leave the
        # rest in their original score order so the report's
        # not-selected list is still meaningful.
        head = select_antigens_for_coverage(
            ranked_vaccine_peptides, target_alleles, cap)
        head_keys = {id(item) for item in head}
        tail = [item for item in ranked_vaccine_peptides
                if id(item) not in head_keys]
        ranked_vaccine_peptides = head + tail

    antigens = antigens_for_assembly(
        ranked_vaccine_peptides, options, 'peptide', source_evidence=source_evidence)
    linker = get_linker(options.linker) if options.antigens_per_construct > 1 else None
    linker_aa = linker.amino_acids if linker else ''
    max_length = (options.antigens_per_construct * options.max_antigen_length_aa
                  + max(0, options.antigens_per_construct - 1) * len(linker_aa))
    def assembled_length(group):
        return sum(len(a.amino_acids) for a in group) + max(0, len(group) - 1) * len(linker_aa)
    packed, exclusions = pack_antigens(antigens, options, assembled_length, max_length, 'aa', selection_audit)
    audits = {}
    for name, _, vp in iter_named_antigens(ranked_vaccine_peptides, options.candidates_per_slot):
        audits.setdefault((name, vp.amino_acids), []).append(vp.window_selection_audit)

    # ``mode`` retained in the manifest for back-compat; ``antigen_content``
    # + ``antigens_per_construct`` are the orthogonal axes that actually
    # drive dispatch.
    base_components = {
        'mode': options.mode,
        'antigen_content': options.antigen_content,
        'antigens_per_construct': options.antigens_per_construct,
        'n_terminal_acetylation': options.n_terminal_acetylation,
        'c_terminal_amidation': options.c_terminal_amidation,
    }

    if linker and linker.inert_in_peptide_mode:
        logger.warning(
            "Linker '%s' is functionally inert in synthesized peptides "
            "(2A skipping is co-translational). Including for vocabulary "
            "symmetry; manifest annotates components.linker_inert=True.", options.linker)
    constructs = []
    for group in packed:
        names = [a.name for a in group]
        sequence = linker_aa.join(a.amino_acids for a in group)
        segments, offset = [], 0
        for antigen in group:
            aa = antigen.amino_acids
            segments.append(dict(name=antigen.name, aa=aa, length_aa=len(aa), start_aa=offset,
                                 end_aa=offset + len(aa), sources=antigen.sources, gene_names=antigen.gene_names,
                                 source_protein_groups=antigen.source_protein_groups))
            offset += len(aa) + len(linker_aa)
        components = dict(base_components)
        if linker:
            components.update(linker=options.linker, linker_inert=linker.inert_in_peptide_mode)
        else:
            audit_queue = audits.get((group[0].sources[0]['candidate_name'], sequence), [])
            audit = audit_queue.pop(0) if audit_queue else {}
            if audit:
                components['window_selection'] = audit
        constructs.append(PeptideConstruct(
            name="peptide_%03d" % (len(constructs) + 1),
            sequence=sequence, antigen_names=names, antigens=segments,
            components=components, manufacturability=_manufacturability_for(sequence)))
    record_assembly_decisions(selection_audit, antigens, constructs, target_alleles, exclusions)
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
    options = options or PeptideConstructConfig()
    constructs = tuple(constructs)
    has_provenance = any(c.construct_placements for c in constructs)
    if has_provenance and not manifest_path:
        raise ValueError("Provenance-bearing peptide output requires a manifest_path sidecar")
    option_mods = options.n_terminal_acetylation or options.c_terminal_amidation
    for c in constructs:
        validate_construct_placements(c.construct_placements, c.sequence, "peptide")
        if c.construct_placements and option_mods:
            raise ValueError(
                "Declare chemical modifications in construct provenance, not writer options")
    has_mods = option_mods or any(
        placement.construct.chemical_modifications
        for c in constructs for placement in c.construct_placements)

    with open(fasta_path, 'w') as f:
        for c in constructs:
            genes = product_gene_names(c)
            f.write(">%s antigens=%s length=%d genes=%s\n" % (
                c.name, ','.join(c.antigen_names), len(c.sequence),
                ','.join(genes)))
            for i in range(0, len(c.sequence), 80):
                f.write(c.sequence[i:i + 80] + "\n")

    if manifest_path:
        manifest = [
            {
                'modality': 'peptide',
                'name': c.name,
                'sequence': c.sequence,
                'length': len(c.sequence),
                'length_unit': 'aa',
                'antigen_names': c.antigen_names,
                **({'antigens': c.antigens} if c.antigens else {}),
                'components': c.components,
                'manufacturability': c.manufacturability,
                **({'construct_provenance': json.loads(to_native_json(c.construct_placements))}
                   if c.construct_placements else {}),
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
            header += ['scale_mg', 'purity_percent', 'counterion',
                       'antigen_names', 'notes']
            if has_provenance:
                header.append('construct_provenance_json')
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
                displayed = _modification_label(
                    options.n_terminal_acetylation, options.c_terminal_amidation,
                    c.sequence) or c.sequence
                row_n_term, row_c_term = n_term, c_term
                if c.construct_placements:
                    n_changes, c_changes = [], []
                    for placement in c.construct_placements:
                        for chemical in placement.construct.chemical_modifications:
                            if chemical.start == chemical.end:
                                final_offset = placement.offset + chemical.start
                                (n_changes if final_offset == 0 else c_changes).append(chemical.name)
                    row_n_term = '; '.join(n_changes) or 'Free'
                    row_c_term = '; '.join(c_changes) or 'Free'
                    displayed = (("[" + row_n_term + "]-") if n_changes else "") + c.sequence
                    if c_changes:
                        displayed += "-[" + row_c_term + "]"
                    notes = (notes + " Explicit construct; see full modification/provenance record. "
                             "Sequence-only predictions do not assess chemical modifications.").strip()
                row = [c.name, c.sequence, len(c.sequence), row_n_term, row_c_term]
                if has_mods:
                    row.append(displayed)
                row += [options.scale_mg, options.purity_percent,
                        options.counterion,
                        ';'.join(c.antigen_names), notes]
                if has_provenance:
                    row.append(to_native_json(c.construct_placements))
                writer.writerow(row)
