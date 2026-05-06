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
from .vaccine_library import (
    get_linker,
    iter_named_antigens,
    select_antigen_window,
    top_mutant_epitopes,
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
                                                     peptide; PGV-001
                                                     canonical)
      content=mutation_spanning, per_construct=N   → Multi-SLP /
                                                     "multi-epitope long"
                                                     concatenation
      content=minimal_epitope,  per_construct=1   → Minimal-epitope
                                                     peptide (one short
                                                     ligand per variant)
      content=minimal_epitope,  per_construct=N   → Concatenated
                                                     minimal-ligand
                                                     pool

    Defaults match the canonical PGV-001 personalized peptide vaccine
    layout: ~20 synthetic long peptides per pool, one antigen per
    peptide, 15-25 aa per peptide.

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
    n_terminal_acetylation: bool = False
    c_terminal_amidation: bool = False
    scale_mg: float = 5.0           # synthesis scale per peptide
    purity_percent: float = 95.0    # HPLC purity target
    counterion: str = "TFA"         # default salt form (TFA / acetate / HCl / free)

    def __post_init__(self):
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


@dataclass
class PeptideConstruct:
    """A single assembled peptide construct."""
    name: str
    sequence: str
    antigen_names: list
    components: dict = field(default_factory=dict)
    manufacturability: dict = field(default_factory=dict)


def _antigen_records(ranked_vaccine_peptides, antigen_content,
                     max_antigen_length_aa, epitopes_per_antigen=1,
                     candidates_per_slot=1):
    """Yield ``(name, amino_acids)`` per antigen.

    Dispatch on ``antigen_content``:
      - ``'mutation_spanning'``: emit one mutation-centered window per
        VaccinePeptide (or per ``candidates_per_slot`` alternates).
      - ``'minimal_epitope'``: emit the top
        ``epitopes_per_antigen`` MHC ligands per VaccinePeptide as
        independent antigens. ``epitopes_per_antigen=1`` (default) is
        the legacy "single top ligand" semantics; >1 packs multiple
        ligands from the same variant.

    Naming + alt-suffix logic comes from ``iter_named_antigens``,
    shared with mRNA assembly so antigen names match across types.
    """
    for base_name, fragment, peptide in iter_named_antigens(
            ranked_vaccine_peptides, candidates_per_slot=candidates_per_slot):
        if antigen_content == "minimal_epitope":
            tops = top_mutant_epitopes(peptide, n=epitopes_per_antigen)
            if not tops:
                logger.info(
                    "Skipping %s in minimal_epitope mode: no mutant "
                    "epitope predictions available.", base_name)
                continue
            for k, ep in enumerate(tops):
                # When epitopes_per_antigen=1 keep the legacy
                # ``<name>_epitope`` suffix; for >1 disambiguate.
                suffix = "_epitope" if len(tops) == 1 else "_epitope%d" % (k + 1)
                yield base_name + suffix, ep.peptide_sequence
        else:
            # mutation_spanning: pick a mutation-centered window
            yield base_name, select_antigen_window(
                fragment, base_name, max_antigen_length_aa)


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
    elif current:
        logger.warning(
            "Dropped final %d antigen(s); --peptide-max-constructs (%d) "
            "reached.", len(current), options.max_constructs)
    return constructs


def assemble_peptide_constructs(
        ranked_vaccine_peptides, options=None, *,
        target_alleles=None):
    """Assemble peptide constructs from ranked vaccine peptides.

    Parameters
    ----------
    ranked_vaccine_peptides : list[(varcode.Variant, list[VaccinePeptide])]
    options : PeptideConstructConfig or None
    target_alleles : optional, list[str]
        Patient MHC alleles. When non-empty, the selector reorders
        the ranked list to maximize allele coverage before packing
        — see :func:`vaxrank.coverage.select_antigens_for_coverage`.
        When empty / None, the input order is preserved (today's
        pure-score behavior). Animal-agnostic: HLA, mouse H-2, swine
        SLA all work the same way.

    Returns
    -------
    list[PeptideConstruct]
    """
    options = options or PeptideConstructConfig()
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

    records = list(_antigen_records(
        ranked_vaccine_peptides, options.antigen_content,
        options.max_antigen_length_aa,
        epitopes_per_antigen=options.epitopes_per_antigen,
        candidates_per_slot=options.candidates_per_slot))
    if not records:
        return []

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

    constructs = []
    if options.antigens_per_construct <= 1:
        # One construct per (variant, candidate) record. Covers both
        # SLP (mutation_spanning, 1 per construct) and minimal-epitope
        # (one short ligand per variant) designs.
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

    # antigens_per_construct > 1: bin-pack into multi-antigen
    # constructs. Works for both content types — concatenated SLPs
    # (the legacy ``multi_epitope`` design) and concatenated minimal
    # ligands (a new combination unlocked by the orthogonal axes).
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
    options = options or PeptideConstructConfig()
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
            header += ['scale_mg', 'purity_percent', 'counterion',
                       'antigen_names', 'notes']
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
                row += [options.scale_mg, options.purity_percent,
                        options.counterion,
                        ';'.join(c.antigen_names), notes]
                writer.writerow(row)
