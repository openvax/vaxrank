# Configuration

Every setting has a compiled-in default. Override it with a YAML file, a
single-key override, or a dedicated CLI flag.

## YAML config file

Store the parameters you repeat on every run in a YAML file:

```sh
vaxrank --config my_config.yaml --input-epitopes candidates.tsv
```

Example `my_config.yaml`:

```yaml
epitopes:
  min_score: 0.00001                        # drop epitopes below this score
  scoring_mode: affinity                    # "affinity" or "percentile_rank"
  logistic_midpoint: 350.0                  # IC50 (nM) at which score = 0.5
  logistic_width: 150.0                     # steepness of logistic curve
  affinity_cutoff: 5000.0                   # IC50 >= this → score 0
  percentile_rank_cutoff: 10.0              # rank >= this → score 0 (percentile mode)

vaccine_peptides:
  preferred_length: 25                      # target amino acids per vaccine peptide
  min_length: 25                            # minimum vaccine peptide length
  max_length: 25                            # maximum vaccine peptide length
  padding_around_mutation: null             # adaptive RNA context; explicit legacy override
  per_mutation: 1                           # peptides to keep per variant
  max_epitopes_per_candidate: 1000          # 0 = keep all
  score_fraction_of_best: 0.99              # drop candidates scoring < 99% of best

peptide:
  manufacturability:                        # GRAVY = mean hydropathy
    max_c_terminal_hydropathy: 1.5          # max GRAVY of C-terminal 7-mer
    min_kmer_hydropathy: 0.0                # min max-7mer GRAVY (floor)
    max_kmer_hydropathy_low_priority: 1.5   # low-priority max-7mer GRAVY cap
    max_kmer_hydropathy_high_priority: 2.5  # high-priority max-7mer GRAVY cap
```

To start from the shipped, commented YAML of supported keys and current
defaults:

```sh
vaxrank --print-default-config > my-config.yaml
```

YAML keys are the public configuration interface; Python dataclass field names
can differ. Epitope filtering and scoring live under `epitopes`, construct
ranking under `vaccine_peptides`, and peptide synthesis thresholds under
`peptide.manufacturability`.

## Overrides and resolution order

CLI arguments override YAML values. `--config-value` overrides individual keys
without editing the file:

```sh
vaxrank --config my_config.yaml \
  --config-value vaccine_peptides.score_fraction_of_best=0.95 \
  --config-value epitopes.percentile_rank_cutoff=5.0
```

Use `--config-text` when the right-hand side should be kept as a raw string
instead of being YAML-parsed.

Values are resolved in this order, later winning:

1. Built-in defaults (see `vaxrank/config/defaults.py`)
2. Saved construct defaults when reloading a native design
3. YAML config files (`--config`, in order)
4. `--config-value` / `--config-text` overrides
5. Dedicated CLI flags (e.g. `--vaccine-peptide-length`)

Native epitope scoring restores its saved policy unless you supply an explicit
epitope configuration. Combining different saved scoring policies requires an
explicit choice. See [saved scoring policies](input-workflows.md#saved-scoring-policies).

## Named policy bundles

`--config builtin:openvax-v1` freezes the existing selection rules; appending
another `--config` derives a policy from it. To inspect what is bundled:

```sh
vaxrank --list-policies
vaxrank --show-policy NAME
```

`--show-policy` prints the resolved formulas and the measurements a policy
requires. [Published comparison policies](published-policies.md) cover LENS,
pVACseq aggregate ordering and experimental TESLA gates.
[Selection policies](selection-policies.md) covers experimental
presentation/processing ensembles, self-content trimming, serum susceptibility
and replayable decision audits.

## RNA reconstruction context

By default the resolved vaccine peptide size K drives a `2*K-1` amino-acid
context request: 29 for 15mers, 49 for 25mers, 59 for 30mers. Actual
RNA-supported context may be shorter. Isovar's `balanced` policy maximizes
mutation-containing peptide windows among candidates retaining at least 85% of
the best candidate's compatible read-name support. Independently, every
retained cDNA base must meet a two-read-object coverage floor.

Compatible reads need not span the whole peptide. Missing RNA context is not
filled with reference sequence. See [Isovar support](isovar-support.md) for
the retained support measurements.

```yaml
isovar:
  protein_sequence_length: null            # derive from peptide size
  protein_context_peptide_length: null     # use vaccine_peptides.preferred_length
  protein_sequence_preference: balanced    # or support / context
  min_protein_sequence_support_fraction: 0.85
  min_variant_sequence_coverage: 2
```

Each setting also has a same-named CLI flag with hyphens in place of
underscores, or can be set with `--config-value isovar.<name>=<value>`.
A coverage floor of zero disables that floor. An explicit
protein length wins over legacy `--padding-around-mutation` or YAML padding;
otherwise an explicit padding request asks for `K + 2*padding`. With neither,
Isovar derives the adaptive target. To reproduce historical Isovar extraction
and support-first selection, pass `--protein-sequence-length 20
--protein-sequence-preference support`.

The target covers all placements of a centered single-residue mutation; wider
mutations may need more context. Deletion windows must span the junction.
Vaxrank rejects candidate windows shorter than `vaccine_peptides.min_length`
(25 by default); lower that minimum only when shorter peptides are part of the
intended design. The RNA support thresholds are unchanged by it. RNA settings
are saved with run arguments and do not rerun reconstruction when rendering
cached reports or ranking external pVACseq/LENS predictions.

For [native Isovar ORF inputs](isovar-inputs.md), use
`isovar.sources.<family>.admission`, `.reconstruction`, `.reads`, and `.protein`.
The same controls work in source JSON, manifest entries, Python
`IsovarSourceOptions`, and `--set` (an alias for `--config-value`). For example,
`--set isovar.sources.variant_orfs.admission.allow_partial=true` permits partial
coding candidates to pass that gate while retaining their uncertainty and
requiring independent tumor-specificity evidence. Reconstruction controls cannot
be applied to a saved producer result. The source guide documents each family's
valid controls, support units, precedence and unresolved outcomes.

The DNA-only fallback remains opt-in, with its own separate padding default. It
examines concrete Varcode outcomes before choosing a usable changed protein,
including alternatives behind unchanged or unresolved fusion candidates. Among
usable outcomes it prefers consequence priority, then protein length and
transcript-pair IDs. Saved fragments retain the chosen candidate's provenance
and protein coordinates. These are sequence predictions; the fallback assigns
no RNA support.

Varcode 11 can attach alternatives to an ordinary annotation: a preferred
intronic comparison can retain a protein-changing exon-skip candidate. Vaxrank
uses Varcode's shared candidate traversal. Reports distinguish the
**producer-preferred effect**, **highest priority effect**, and **candidate
effects**. Producer ordering and consequence severity are not calibrated
probabilities or proof that a product was expressed.

When Isovar supplies an annotation, fallback retains its patient baseline,
phase and RNA hypotheses instead of reannotating the bare variant. An unavailable
baseline or protein remains unavailable. Python callers can supply the same
context with `MutantProteinFragment.from_variant_dna(..., effects=effects)`;
an empty list produces no fragment. Saved fragments and their proteome records
retain the annotation and selected-path evidence through JSON replay.

## Custom filtering and scoring with the Topiary DSL

For anything beyond the scalar logistic and percentile-rank defaults, set
`epitopes.filter_expr` and/or `epitopes.score_expr` to a
[Topiary](https://github.com/openvax/topiary) DSL string. Both accept the full
Topiary 5.0 expression grammar: kind accessors such as `affinity` or
`presentation`, arithmetic, `&` and `|`, `.logistic(...)` and `.clip(...)`
transforms, `column(col_name)` for raw DataFrame columns, and so on.

```yaml
epitopes:
  # Drop rows wholesale before scoring
  filter_expr: "affinity <= 500 & affinity.rank <= 2.0"
  # Compute a per-(peptide, allele) score in [0, 1] (binder-quality score)
  score_expr:  "affinity.logistic_normalized(350, 150)"
```

When `filter_expr` is omitted, no rows are dropped up front. The default
`score_expr` is synthesized from the scalar fields above
and sets `ic50 >= affinity_cutoff` to zero.

Use `affinity.logistic_normalized(m, w)` for a `[0, 1]` binder-quality score.
The plain `affinity.logistic(m, w)` is the raw sigmoid and caps below 1, about
0.912 at the default `m=350, w=150`.

Invalid DSL strings are rejected at config load rather than mid-pipeline, so
typos in the YAML surface before any predictions run.
