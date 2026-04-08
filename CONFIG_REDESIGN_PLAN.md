# Generalized YAML Config and CLI Redesign Plan

## Summary

This document proposes a new configuration architecture for Vaxrank that:

- uses one generalized YAML schema for epitope, mutation, and vaccine peptide selection
- supports Topiary filter/ranking expressions for epitope-level logic
- separates filtering from scoring/ranking at every selection layer
- moves hard-coded thresholds, tie-breakers, and normalization formulas into config
- preserves backward compatibility with the current CLI and current YAML keys during migration

The immediate goal is not to add more one-off flags. The goal is to make the selection logic explicit, inspectable, versioned, and configurable.

## Current Problems

Today the logic is split across:

- `vaxrank/epitope_config.py`
- `vaxrank/vaccine_config.py`
- `vaxrank/cli/epitope_config_args.py`
- `vaxrank/cli/vaccine_config_args.py`
- `vaxrank/epitope_logic.py`
- `vaxrank/core_logic.py`
- `vaxrank/vaccine_peptide.py`
- `vaxrank/epitope_prediction.py`

This creates a few structural problems:

- YAML only covers a small subset of the actual selection logic.
- CLI flags and YAML keys are thin wrappers around hard-coded behaviors.
- Epitope scoring is hard-coded inside `EpitopePrediction.logistic_epitope_score`.
- Vaccine peptide ranking is hard-coded as a lexicographic tuple in `VaccinePeptide.lexicographic_sort_key`.
- Variant ranking is hard-coded as `top_vaccine_peptide.combined_score`.
- The "keep candidates within 99% of best score" rule is hard-coded in `core_logic.py`.
- Manufacturability thresholds and priorities are hard-coded in `VaccinePeptide.peptide_synthesis_difficulty_score_tuple`.
- Expression scoring uses a hard-coded `sqrt(n_alt_reads)` transform.
- Imported epitopes and internally predicted epitopes do not flow through one general ranking/filtering model.

## Design Goals

1. One canonical config schema for all selection logic.
2. Clear separation between:
   - candidate generation
   - hard filtering
   - soft scoring
   - ranking
   - truncation / top-k selection
3. Topiary expressions should be first-class for epitope-level filtering and ranking.
4. Mutation-level and vaccine-peptide-level ranking should use the same mental model as epitope-level ranking, even if the evaluator is Vaxrank-owned.
5. No raw Python `eval`.
6. Default behavior should remain equivalent to the current pipeline until the user opts into new logic.
7. Existing CLI flags and existing YAML keys should continue to work through a compatibility layer.

## Proposed Config Shape

Use a versioned top-level schema:

```yaml
schema_version: 2

inputs:
  source: isovar

self_proteome:
  source: genome
  exclude_gene_ids: []
  exclude_pirlygenes_cta: false
  exclude_fasta: null
  kmer_lengths:
    derive_from_epitopes: true
    min: null
    max: null

epitopes:
  source: predictor
  prediction:
    topiary:
      ranking: "affinity <= 500 | presentation.rank <= 2"
      rank_by: "0.7 * (1 - affinity.norm(mean=500, std=200)) + 0.3 * presentation.score.norm(mean=0.5, std=0.2)"
  filters:
    expression: "overlaps_mutation and not occurs_in_reference"
  scoring:
    derived_fields:
      affinity_score:
        kind: logistic
        input: ic50
        midpoint: 350
        width: 150
        cutoff: 5000
        higher_is_better: true
      percentile_score:
        kind: linear_window
        input: percentile_rank
        best: 0
        worst: 10
        higher_is_better: true
  keep:
    top_n_per_candidate: null

mutations:
  filters:
    expression: "passes_isovar_filters and has_vaccine_peptides"
  ranking:
    score: "best_vaccine_peptide.score"
    tie_break:
      - expr: "sqrt(n_alt_reads)"
        order: desc
  keep:
    top_n: null

vaccine_peptides:
  generation:
    lengths: [25]
    padding_around_mutation: 5
  filters:
    expression: "contains_mutant_epitopes"
  ranking:
    score: "mutant_epitope_score * sqrt(n_alt_reads) - 0.25 * wildtype_epitope_score"
    tie_break:
      - expr: "manufacturability.cysteine_count"
        order: asc
      - expr: "manufacturability.cterm_7mer_gravy"
        order: asc
      - expr: "mutation_distance_from_edge"
        order: desc
  keep:
    per_mutation: 1
    score_fraction_of_best: 0.99

manufacturability:
  hydropathy_window_size: 7
  hard_filters: []
  penalties:
    max_cterm_gravy: 1.5
    max_kmer_gravy_high_priority: 2.5
    max_kmer_gravy_low_priority: 1.5
    min_kmer_gravy: 0.0

compat:
  warn_on_legacy_keys: true
```

## Key Model Change

The new model should be:

- generate candidate rows for epitopes
- annotate those rows with derived features and normalization outputs
- filter and rank epitopes
- aggregate epitope rows into vaccine peptide candidate rows
- filter and rank vaccine peptide candidates
- aggregate vaccine peptide candidates into mutation rows
- filter and rank mutations

That is a better fit than pushing selection behavior into object methods with hidden constants.

## Expression Model

### Epitope level

Use Topiary expressions as the canonical epitope filter/rank language.

Vaxrank should support:

- `epitopes.prediction.topiary.ranking`
- `epitopes.prediction.topiary.rank_by`

These should map directly onto Topiary's existing expression/ranking support, for example:

- `affinity <= 500 | presentation.rank <= 2`
- `0.5 * (1 - affinity.norm(mean=500, std=200)) + 0.5 * presentation.score.norm(mean=0.6, std=0.2)`

This gives Vaxrank a real expression language immediately instead of inventing a new one for epitope scoring.

### Mutation and vaccine peptide levels

Topiary's current DSL is epitope-centric. For mutation and vaccine-peptide ranking, Vaxrank should introduce its own safe expression layer with the same operator/function style:

- arithmetic: `+`, `-`, `*`, `/`
- logical filters: `and`, `or`, `not`
- numeric transforms: `sqrt`, `log`, `log10`, `clip`
- configurable normalizers: `logistic_norm`, `linear_norm`, `gauss_norm`
- null-safe helpers: `coalesce`, `isnan`, `exists`

The syntax does not need to be identical to Topiary internals, but it should feel similar enough that users can move between layers without changing mental models.

## Derived Fields and Normalizers

Normalization should move out of `EpitopePrediction.logistic_epitope_score` and into config-defined derived fields.

Examples:

- `affinity_score = logistic_norm(ic50, midpoint=350, width=150, cutoff=5000, invert=true)`
- `percentile_score = linear_norm(percentile_rank, best=0, worst=10, clamp=true, invert=true)`
- `expression_score = sqrt(n_alt_reads)`
- `combined_score = mutant_epitope_score * expression_score`

The code should still ship with defaults, but those defaults should live in config defaults, not hidden formulas spread across classes.

## Hard-Coded Behavior To Externalize

| Current behavior | Current location | Target config |
|---|---|---|
| IC50 logistic midpoint `350.0` | `epitope_config.py`, `epitope_prediction.py` | `epitopes.scoring.derived_fields.affinity_score.midpoint` |
| IC50 logistic width `150.0` | `epitope_config.py`, `epitope_prediction.py` | `epitopes.scoring.derived_fields.affinity_score.width` |
| IC50 cutoff `5000.0` | `epitope_config.py`, `epitope_prediction.py` | `epitopes.scoring.derived_fields.affinity_score.cutoff` |
| Percentile rank cutoff `10.0` | `epitope_prediction.py` | `epitopes.scoring.derived_fields.percentile_score.worst` |
| Minimum epitope score | `epitope_config.py` | `epitopes.filters.expression` or explicit keep threshold |
| `sqrt(n_alt_reads)` expression score | `vaccine_peptide.py` | `vaccine_peptides.ranking.score` and `mutations.ranking.score` |
| `combined_score = expression_score * mutant_epitope_score` | `vaccine_peptide.py` | `vaccine_peptides.ranking.score` |
| Keep candidates within `0.99` of best | `core_logic.py` | `vaccine_peptides.keep.score_fraction_of_best` |
| Variant ranking by top peptide score | `core_logic.py` | `mutations.ranking.score` |
| Manufacturability hydropathy thresholds `1.5`, `2.5`, `1.5`, `0.0` | `vaccine_peptide.py` | `manufacturability.penalties.*` |
| Manufacturability tuple ordering | `vaccine_peptide.py` | `vaccine_peptides.ranking.tie_break` |
| Top epitopes to keep per peptide | `vaccine_config.py`, `vaccine_peptide.py` | `epitopes.keep.top_n_per_candidate` or `vaccine_peptides.keep.max_epitopes_per_candidate` |
| Vaccine peptide length and mutation padding | `vaccine_config.py` | `vaccine_peptides.generation.*` |

## CLI Strategy

Do not solve this by adding a dedicated flag for every new config field.

The CLI should have three layers:

1. Stable config entrypoint
   - `--config FILE`

2. General-purpose override layer
   - `--set section.path=value`
   - `--expr section.path='expression string'`

3. Backward-compatible convenience flags
   - keep current flags like `--min-epitope-score`, `--vaccine-peptide-length`, `--padding-around-mutation`
   - map them into the new schema
   - emit deprecation warnings when a flag only exists as a legacy shim

This prevents the CLI from turning into another hard-coded schema.

## Internal Architecture

Introduce a real config subsystem:

- `vaxrank/config/schema.py`
  - msgspec structs for the new schema
- `vaxrank/config/defaults.py`
  - central defaults for all tunables
- `vaxrank/config/loader.py`
  - YAML load, validation, override merge
- `vaxrank/config/compat.py`
  - map legacy CLI args and legacy YAML keys into schema v2
- `vaxrank/config/expr.py`
  - safe evaluator for mutation/vaccine expressions
- `vaxrank/selection/epitopes.py`
  - epitope row normalization, Topiary integration, epitope filtering/ranking
- `vaxrank/selection/vaccine_peptides.py`
  - candidate generation plus filter/rank/keep logic
- `vaxrank/selection/mutations.py`
  - mutation-level aggregation and ranking

The existing `EpitopeConfig` and `VaccineConfig` can remain temporarily as adapters that are populated from the new schema during the migration window.

## Imported Epitope Inputs

External epitope inputs should be normalized to the same row schema as internally predicted epitopes before selection runs.

That means:

- imported pVACseq / LENS / native Vaxrank epitopes should flow through one epitope selection path
- missing fields should become explicit `NaN` / `None`
- Topiary ranking/filtering should operate on that normalized table where possible
- Vaxrank-specific fields like `occurs_in_reference` and `overlaps_mutation` should be added after import normalization

## Implementation Phases

### Phase 1: Freeze schema and compatibility rules

- define schema v2 structs
- define legacy-to-v2 key mapping
- define default-equivalent config that reproduces current behavior
- add schema tests and merge tests

### Phase 2: Centralize config loading

- replace direct YAML decoding in `epitope_config_args.py` and `vaccine_config_args.py`
- load one unified config object in `entry_point.py`
- thread one config object through the pipeline

### Phase 3: Epitope layer redesign

- normalize all epitope predictions into one DataFrame-like row model
- integrate Topiary `ranking` / `rank_by`
- replace `logistic_epitope_score` call sites with derived-field evaluation
- preserve current defaults exactly

### Phase 4: Vaccine peptide candidate redesign

- materialize vaccine peptide candidate rows with explicit feature columns
- move filter/ranking/truncation into config-driven logic
- replace `lexicographic_sort_key` hard-coding with ordered configured criteria

### Phase 5: Mutation-level redesign

- materialize mutation-level rows with aggregated features
- make mutation inclusion and ranking configurable
- replace `top_vaccine_peptide.combined_score` as the only ranking rule

### Phase 6: CLI redesign

- add `--set` and `--expr`
- keep legacy flags working
- emit targeted deprecation warnings
- update help text to point advanced users toward YAML plus expressions

### Phase 7: Reports and serialization

- ensure reports read materialized derived scores instead of recomputing hidden defaults
- version saved JSON output with schema version
- ensure backward loading of old JSON still works

### Phase 8: Tests and docs

- config parsing and validation tests
- legacy CLI parity tests
- Topiary expression integration tests
- mutation/vaccine expression tests
- parity tests proving old defaults reproduce current ranking
- new user docs with 3 or 4 realistic config examples

## Acceptance Criteria

The redesign is complete when:

- the full ranking/filtering behavior is represented in YAML rather than hidden in methods
- epitope-level ranking can use Topiary expressions directly
- mutation-level and vaccine-peptide-level filtering/ranking are config-driven
- current default behavior is preserved by a versioned default config
- current CLI flags still work during the migration period
- reports and saved outputs reflect the configured logic, not baked-in formulas

## Non-Goals

- arbitrary Python in config
- supporting every historical internal name forever without warnings
- rewriting all report templates before the selection engine is stabilized

## Recommended First Slice

The best first implementation slice is:

1. Add schema v2 and unified config loading.
2. Materialize epitope derived fields and move current logistic scoring parameters into config defaults.
3. Add Topiary `ranking` / `rank_by` passthrough for epitope selection.
4. Keep current vaccine and mutation ranking logic as the default compatibility behavior.
5. Only after that, generalize vaccine-peptide and mutation ranking expressions.

That keeps the first landing small enough to verify while still establishing the right architecture.
