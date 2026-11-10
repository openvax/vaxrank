# Vaccine designs and construct outputs

**Experimental.** Choose the vaccine type first, then antigen content and
antigens per construct:

| Axis | Values | What it controls |
|---|---|---|
| `--vaccine-type` | `peptide` / `mrna` (multi-valued) | The platform(s); pass both for parallel design |
| `--antigen-content` | `mutation_spanning` / `minimal_epitope` | What each antigen *is* |
| `--peptide-antigens-per-construct` / `--mrna-antigens-per-construct` | `1` / `N` | How many antigens to concatenate per construct |

That matrix yields eight designs, four per vaccine type:

| Type | Content | Per-construct | Design name | Reference |
|---|---|---|---|---|
| peptide | mutation_spanning | 1 | SLP (default) | [PGV-001 (Saxena 2025)](https://pubmed.ncbi.nlm.nih.gov/40094414/) |
| peptide | mutation_spanning | N | Multi-SLP / multi-epitope long peptide | |
| peptide | minimal_epitope | 1 | Minimal-ligand peptide | |
| peptide | minimal_epitope | N | Concatenated minimal-ligand peptide | |
| mrna | mutation_spanning | N | Multi-antigen mRNA (default for mRNA) | [Sahin 2017](https://doi.org/10.1038/nature23003) / [Rojas 2023](https://doi.org/10.1038/s41586-023-06063-y) |
| mrna | mutation_spanning | 1 | Single-antigen mRNA | |
| mrna | minimal_epitope | N | "String of beads" mRNA | [Whitton 1993](https://pubmed.ncbi.nlm.nih.gov/7677954/) |
| mrna | minimal_epitope | 1 | Single-ligand mRNA | |

`--epitopes-per-antigen` sets how many top ligands to take from each ranked
vaccine peptide in `minimal_epitope` mode (default `1`).

## Distinct sequences and source mappings

Peptide and mRNA assembly group identical **emitted antigen sequences** after
window selection or minimal-epitope selection, before packing or applying
construct caps. Ten candidate IDs producing the same sequence occupy one
antigen slot. The first candidate in selection order supplies its display name;
all contributing source mappings remain separate. Scores and RNA support are
not added across those mappings.

`--max-windows-per-antigen-source N` retains up to **3 non-overlapping source
windows** by default. A source is a recorded mutation/event or an identified
translated antigen, not a gene-name spelling or a sample. Boundaries may touch;
windows cannot share a residue on their parent sequence. Selection takes the
best remaining region, removes overlaps, and repeats. The near-best score
comparison restarts each round. A repeated region with no new target
peptide/allele content does not earn another window. Each retained window then
competes for vaccine capacity on its own configured combined score.

Manufacturing shifts belong to one selected region. Their recorded alternatives
must preserve all selected target peptide/allele occurrences and avoid the other
selected regions. They do not consume extra windows or slots. Peptide refinement
enforces this target-preservation contract. Scores and eligibility still come
from the configured DSL and window policy; selection is deterministic and greedy,
not a global optimization or an estimate of immune response.

`--peptide-candidates-per-slot`, `--mrna-candidates-per-slot`, and
`--max-vaccine-peptides-per-variant` apply to legacy overlapping-candidate
selection (`vaccine_peptides.max_windows_per_antigen_source: null`). Saved
pre-3.53 designs and the frozen `openvax-v1` profile retain that behavior. New
biological windows all reach assembly without raising `candidates-per-slot`.
Available sequence/prediction context limits the windows that can be selected;
a previously truncated input cannot recover omitted regions.

Set goals and caps with a separate flag for each antigen type:

| Type | Minimum goal | Maximum cap |
|---|---|---|
| Mutation | `--min-mutation-antigens N` | `--max-mutation-antigens N` |
| Fusion | `--min-fusion-antigens N` | `--max-fusion-antigens N` |
| SV | `--min-sv-antigens N` | `--max-sv-antigens N` |
| CTA | `--min-cta-antigens N` | `--max-cta-antigens N` |
| ERV | `--min-erv-antigens N` | `--max-erv-antigens N` |
| Splice | `--min-splice-antigens N` | `--max-splice-antigens N` |
| Viral | `--min-viral-antigens N` | `--max-viral-antigens N` |
| Isoform | `--min-isoform-antigens N` | `--max-isoform-antigens N` |

Counts are **unique antigen sequences across the whole vaccine**, including
all mRNA constructs. They count selected source windows in long-antigen mode
and individual epitopes in minimal-epitope mode. A shared sequence counts once
under its first ranked source type and retains all other sources. These are
category totals, distinct from `--max-windows-per-antigen-source`, which caps
the number of windows taken from each source protein/event.

Minimums default to zero and are goals subject to eligibility and capacity;
maximums default to no type-specific cap and are hard limits. A maximum of zero
excludes that type from the vaccine. Actual counts and unmet goals appear in
reports and saved decisions. A CLI limit overrides only that type's configured
value; other YAML limits remain in effect.

`--antigen-allocation` chooses **which eligible sequences fill the available
slots within those limits**. Both modes serve minimum goals first and enforce
every maximum. `balanced` (the default for fresh designs) then favors the least
represented eligible type, using the existing score/allele-coverage order for
ties and within each type. `ranked` uses that existing order directly after
serving the goals. Neither mode changes scores or eligibility. Unused capacity
passes to other eligible types that have room under their caps.

For example, suppose ten slots are available, all CTA candidates rank ahead
of mutations, mutations rank ahead of viral candidates, and each type has
enough eligible candidates. With `--min-cta-antigens 2 --max-cta-antigens 4`:

| Mode | CTA | Mutation | Viral |
|---|---:|---:|---:|
| `ranked` | 4 | 6 | 0 |
| `balanced` | 4 | 3 | 3 |

The same limits admit different designs because they do not uniquely specify
the final category counts. The mode supplies the remaining selection rule.

Equivalent Python/YAML settings are `VaccineConfig.max_windows_per_antigen_source`
and `VaccineConfig.antigen_allocation`; both construct option classes also accept
`antigen_allocation`. For example:

```yaml
vaccine_peptides:
  max_windows_per_antigen_source: 3
  antigen_allocation:
    mode: balanced
    min_per_kind: {mutation: 5}
    max_per_kind: {CTA: 8, viral: 4}
```

Supported type names are `mutation`, `fusion`, `SV`, `CTA`, `ERV`, `splice`, and
`viral`, plus explicitly admitted `isoform` inputs. ORF context is retained as
provenance, not an extra category. These allocation settings are design
preferences, not scientifically established optimal proportions.

Both modality manifests contain `antigens`: exact amino-acid sequences, lengths,
zero-based half-open positions in the peptide or mRNA's encoded protein, and a
`sources` list. Source records retain candidate and prediction IDs, genes,
transcripts, proteins, species, input IDs, and recorded event/evidence-group
identities. The Python products expose the same records as `.antigens` and a
derived `.antigen_sources` mapping. Full prediction and reconstruction evidence
remains in the native dataset and its assembly source graph.

Different sequences are not merged across sources. In legacy selection, several
ncORF windows from one event are alternative sequences, not independent target
events. `run_summary.txt` counts actual products, distinct antigen sequences,
recorded source events and candidate mappings without known events separately.
An evidence group identifies shared reconstruction evidence; it is not itself a
claim of an independent biological event.

Native designs saved before 3.47.4 with duplicate slots are rebuilt on replay.
That may require the configured junction predictor again. Save the corrected
native output to replay it offline afterward. Historical designs without
duplicate slots still restore their original products unchanged; their older
manifests may lack the newer source mappings.

The commands below start from `candidates.tsv`, a native export with admitted
antigens and retained sequence/annotation evidence. See [input workflows](input-workflows.md)
to create one.

## Peptide designs

SLP (the default): mutation-spanning long peptide, one antigen per construct,
as in the PGV-001 design.

```sh
vaxrank --input-epitopes candidates.tsv --vaccine-type peptide --output-dir vaccine_out/
# → vaccine_out/{vaccine.fasta, manifest.json, order_form.csv}
```

Concatenate five antigens with an `AAY` spacer:

```sh
vaxrank --input-epitopes candidates.tsv --vaccine-type peptide \
        --output-dir vaccine_out/ \
        --peptide-antigens-per-construct 5 --peptide-linker AAY
```

Select a single short MHC ligand per construct:

```sh
vaxrank --input-epitopes candidates.tsv --vaccine-type peptide \
        --output-dir vaccine_out/ \
        --antigen-content minimal_epitope
```

## mRNA designs

Multi-antigen mRNA is the default for `--vaccine-type mrna`. Antigens are linked with `(G4S)2` and
emitted as CDS, no-polyA and full (with polyA) FASTAs plus a structured
manifest.

```sh
vaxrank --input-epitopes candidates.tsv --vaccine-type mrna --output-dir mrna_out/
# → mrna_out/{cds.fasta, no_polyA.fasta, full.fasta, manifest.json, mrna-sequence-parts.csv}
```

String-of-beads mRNA: short MHC ligands linked together rather than
mutation-spanning windows.

```sh
vaxrank --input-epitopes candidates.tsv --vaccine-type mrna --output-dir out/ \
        --mrna-antigen-content minimal_epitope --mrna-antigens-per-construct 8 \
        --mrna-linker AAY
```

Top-N ligands per variant in a string of beads: pack several top MHC ligands
from each ranked vaccine peptide as separate antigens.

```sh
vaxrank --input-epitopes candidates.tsv --vaccine-type mrna --output-dir out/ \
        --mrna-antigen-content minimal_epitope \
        --mrna-epitopes-per-antigen 2 --mrna-antigens-per-construct 16
```

## Output modes

`--vaccine-type` defaults to `peptide mrna`. Specify one type to emit only that modality. `--output-dir` enables construct files:
a single vaccine type writes directly into the directory, while multiple types
write into `peptide/` and `mrna/` subdirectories. Analysis reports have their
own `--output-*` flags and can be requested without writing construct files.

```sh
# Both at once (per-modality subdirs in mixed_out/)
vaxrank --input-epitopes candidates.tsv --vaccine-type peptide mrna --output-dir mixed_out/

# Reports only (no vaccine constructs)
vaxrank --input-epitopes candidates.tsv --output-pdf-report report.pdf

# Drive vaccine design from a pre-computed LENS report
vaxrank --input-lens patient.lens.tsv --vaccine-type mrna \
        --output-dir mrna_out/ --ensembl-release 102
```

One text/PDF report in the run directory starts with **input coverage and emitted
products**, then compact candidate rankings linked to those products. Product sequences, lengths,
antigen placements and source counts come from the actual assembly outputs.
Candidate windows may be longer than a selected minimal epitope or differ from a
trimmed product. Candidate prediction scores do not establish final-product or
junction coverage. Per-modality directories contain construct files; they do not
repeat the report.

Add `--report-candidate-details` to append grouped source windows and shared
epitope tables. Source IDs link to their complete identities, and tables retain
all source memberships. Differing support, context or predictions remain separate.
This appendix labels target and self epitopes consistently for mutations,
noncanonical ORFs and other admitted sources. Full producer JSON stays in the
native dataset. Compact reports need no annotation queries; the detailed mutation
appendix uses retained annotation or needs `--ensembl-release` when it is missing.
This holds equally for single files, external inputs and manifests.

The Python report API accepts `TemplateDataCreator(..., products={"peptide":
peptides, "mrna": rnas}, candidate_details=False)`, where values are the emitted
construct lists. An empty list means a successful run with no products; omitting
`products` produces an analysis-only report. This replaces the old
`vaccine_constructions`, `mrna_ranking_decisions` and `target_alleles` report
arguments: selection forecasts are no longer displayed as final products.

Pass `selection=dataset.selection` to `TemplateDataCreator` to show the saved
effective configuration and decision audits. The CLI passes this automatically.
Omitting it leaves unavailable settings explicitly unknown; parser arguments
are never treated as authoritative assembly settings. Library assemblers accept
an optional `selection_audit={}` dictionary, filled with deduplication and
packing outcomes. Reports are read-only views of these records and products.


## Peptide construct files

`vaccine.fasta`, `manifest.json` and `order_form.csv`, written into
`--output-dir` (or `--output-dir/peptide/` in multi-mode). The peptide design
comes from `--antigen-content` and `--peptide-antigens-per-construct`.

| Flag | Purpose |
|---|---|
| `--output-dir` | Where to write the construct files |
| `--peptide-linker` | Inter-antigen spacer (e.g. `AAY`, `(G4S)2`); default `G4S3` |
| `--peptide-max-antigen-length-aa` | Truncate antigens longer than this |
| `--peptide-n-terminal-acetyl` | Add an N-terminal acetylation note to the manifest |
| `--peptide-c-terminal-amide` | Add a C-terminal amide note to the manifest |

## mRNA construct files

A directory holding three FASTAs (`cds.fasta`, `no_polyA.fasta`, `full.fasta`)
plus `manifest.json` (per-element view) and `mrna-sequence-parts.csv` (long-format
per-element table with AA and nt). Codon optimization uses
[DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel); 2A
self-cleaving peptides keep their published codon usage automatically.

Construct anatomy:

| Flag | Purpose |
|---|---|
| `--output-dir` | Where to write construct files (or `--output-dir/mrna/` in multi-mode) |
| `--mrna-signal-peptide` | Leader peptide: `HLA-A`, `HLA-B`, `tPA`, `IgK`, `CD8A`, `CD28` |
| `--mrna-5p-utr` | 5' UTR (e.g. `HBB`, `HBB_FI` tandem) |
| `--mrna-3p-utr` | 3' UTR |
| `--mrna-include-mitd` / `--mrna-no-mitd` | Include the BioNTech MITD trafficking domain |
| `--mrna-mitd` | Which MITD variant (`HLA-A` / `HLA-B`) |
| `--mrna-max-length-nt` | Hard cap on construct length (nt) |
| `--mrna-antigens-per-construct` | Antigens per CDS |
| `--mrna-max-constructs` | Stop emitting after this many constructs |

PolyA tail:

| Flag | Purpose |
|---|---|
| `--mrna-poly-a-length` | Length of polyA tail (default `120`) |
| `--mrna-poly-a-segmented` | Use the BNT162b2-style segmented pattern (A30 + linker + A70) |
| `--mrna-poly-a-first-segment` | Length of the first segment when segmented |
| `--mrna-poly-a-segment-linker` | Inter-segment linker sequence |

Codon optimization:

| Flag | Purpose |
|---|---|
| `--mrna-codon-species` | Target organism for codon usage (default `h_sapiens`) |
| `--mrna-codon-method` | DnaChisel optimization strategy |
| `--mrna-csv-no-full-rows` | Skip the `full` polyA'd rows in `mrna-sequence-parts.csv` (saves disk space) |

### Junction linker optimization

A per-junction, MHC-aware linker swap minimizes predicted presentation of
chimeric k-mers spanning antigen junctions.

| Flag | Purpose |
|---|---|
| `--mrna-linker` | Default inter-antigen spacer (e.g. `(G4S)2`) |
| `--mrna-optimize-linkers` / `--mrna-no-optimize-linkers` | Require / disable junction prediction; auto by default |
| `--mrna-junction-predictor` | One MHCtools model for new junction queries; leaves candidate scores unchanged |
| `--mrna-junction-alleles` | Query alleles; defaults to the declared manifest genotype or VCF/BAM `--mhc-alleles` |
| `--mrna-junction-predictor-path` / `--mrna-junction-predictor-models-path` | Independent executable / weights locations |
| `--mrna-junction-candidates` | Candidate linkers considered at each junction |
| `--mrna-junction-rank-strong` | Strong-binder %-rank threshold |
| `--mrna-junction-rank-mild` | Mild-binder %-rank threshold |

VCF + BAM runs reuse a single configured candidate predictor in auto mode. With
multiple candidate predictors they keep the shared linker; select
`--mrna-junction-predictor` to enable junction optimization. External-report
runs use the shared linker without new predictions unless
`--mrna-junction-predictor` is supplied; if those reports declare no genotype,
also supply `--mrna-junction-alleles`, since observed report coverage is not
treated as a genotype. `--mrna-no-optimize-linkers` prevents all junction model
loading.

Cached JSON rerendering (`--input-json-file`) also needs explicit junction model
and allele options for new predictions. YAML equivalents live under `mrna:`;
`optimize_linkers: null` selects auto mode.

The rank policy uses affinity when available, otherwise presentation. It does
not add both endpoint counts or treat allele-independent processing and
haplotype-only results as allele-specific hits. The manifest records the policy,
model names/versions, alleles, lengths and peptide queries under
`elements.junction_swap`; native exports retain it for replay. Disable optional
cleavage annotation with `--no-processing-aware-annotation`.

## Manifest schema

The peptide and mRNA construct JSON manifests share a back-compatible schema
(`modality`, `name`, `length`, `length_unit`, `antigen_names`, `components`,
`manufacturability`). The mRNA manifest additionally exposes `cds`,
`no_polya_nt`, `full_nt`, per-antigen `antigens` (each with AA and nt), and a
structured `elements` dict with one entry per layer: 5' UTR, signal peptide,
antigens, linkers per junction, MITD, stop codon, 3' UTR, polyA. Each layer
carries both AA (where applicable) and nt forms for direct inspection.

## Shared linker library and grammar

Both vaccine types consume the same set of linker names, so one construct
design can be ported between peptide and mRNA backbones.

Static entries:

| Name | Type | Use |
|---|---|---|
| `G2S`, `G3S`, `G4S`, `G5S` | flexible (Gly_n_Ser) | The (Gly4Ser)n family ([Huston *PNAS* 1988](https://doi.org/10.1073/pnas.85.16.5879)); used clinically in BioNTech FixVac / iNeST as `(G4S)2` |
| `EAAAK` | rigid α-helical | When fused antigens need separation rather than flex ([Arai *Protein Eng* 2001](https://doi.org/10.1093/protein/14.8.529)) |
| `RKRR`, `RVKR`, `RKRKR` | furin cleavage | R-X-(K/R)-R motif ([Hosaka *J Biol Chem* 1991](https://pubmed.ncbi.nlm.nih.gov/1905715/)); see the cited cleavage study |
| `AAY` | proteasome-friendly | [Livingston 2001](https://doi.org/10.1016/S0264-410X(01)00233-X); spacer comparisons in [Yang 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514284/) and [Aguilar-Gurrieri 2023](https://doi.org/10.1007/s00262-023-03409-3) |
| `AAA` | alanine spacer | [Aguilar-Gurrieri *Cancer Immunol Immunother* 2023](https://doi.org/10.1007/s00262-023-03409-3): alanine spacers improved SIINFEKL presentation in the tested cell system |
| `GPGPG` | helper-T spacer | Between MHC-II epitopes ([Livingston *J Immunol* 2002](https://doi.org/10.4049/jimmunol.168.11.5499)) |
| `P2A`, `T2A`, `F2A`, `E2A` | self-cleaving 2A | Co-translational ribosomal skipping for mRNA constructs ([Donnelly *J Gen Virol* 2001](https://doi.org/10.1099/0022-1317-82-5-1027); [Kim *PLoS ONE* 2011](https://doi.org/10.1371/journal.pone.0018556)). In peptide mode these are functionally inert and the manifest annotates them as such. |

Compositional forms, parsed at lookup time:

| Form | Meaning | Example |
|---|---|---|
| `(BASE)N` / `(BASE)xN` / `BASExN` | Repeat N times | `(G4S)2` → `GGGGSGGGGS`, `G4Sx2` → same |
| `GnSm` | Literal n glycines + m serines (a single unit, not a repeat) | `G6S` → `GGGGGGS`, `G4S2` → `GGGGSS` |
| `AnY` | n alanines + tyrosine | `A3Y` → `AAAY` |
| `An` | n alanines (no Y) | `A4` → `AAAA` |
| `Gn` | n glycines (no S) | `G4` → `GGGG` |

Repeat counts are capped at 100. The 2A entries are codon-frozen and
positional, so they are rejected in repeat forms; use the base linker once.

`vaccine_library.get_linker(name)` returns a `Linker` with sequence and
primary-source citations. Junction optimization considers `G3S`, `G4S`,
`(G3S)2`, `(G4S)2` and `AAA` by default.

## Legacy peptide-mode flag

`--peptide-mode {slp, minimal_epitope, multi_epitope}` is kept for older
scripts. It is shorthand for the orthogonal axes:

- `slp` ≡ `--antigen-content mutation_spanning --peptide-antigens-per-construct 1`
- `minimal_epitope` ≡ `--antigen-content minimal_epitope --peptide-antigens-per-construct 1`
- `multi_epitope` ≡ `--antigen-content mutation_spanning --peptide-antigens-per-construct N`

Prefer the explicit content and count flags for new designs.
