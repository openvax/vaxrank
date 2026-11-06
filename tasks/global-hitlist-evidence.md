# Global local Hitlist evidence (advances #745)

## Contract

Add `--hitlist-local-evidence` and the corresponding Python keyword to opt
into read-only local MS lookup across candidate inputs. Preserve automatic CTA
lookup and explicit CTA bundle precedence. Do not acquire data, rebuild indexes,
or infer expression, target specificity, or patient-allele restriction from
any-MHC observations. Use the existing exact peptide/species/restriction facets.

Enrich candidate prediction rows before DSL filtering/scoring, after occurrence
identity has been established. Retain raw queries, index identity, returned rows,
and unavailable states in native selection/source evidence. Previously annotated
rows are frozen on replay, including unavailable queries; new rows can be queried
when explicitly opted in. Direct variants use a pre-scoring annotation callback
and retain per-query captures rather than claiming a combined Cartesian scope.
Legacy selected-vaccine JSON cannot be rescored and must reject the new option.

Regression testing reproduced open #804: constant scores restore alleles removed
by the filter. Include the prerequisite fix using Topiary filter decisions and
cover constant scores, allele-free evidence, partial genotype and native replay.

## Plan

- [x] Extract reusable annotation helpers without changing CTA defaults.
- [x] Wire shared scoring and direct prediction, CLI/API, native source evidence.
- [x] Cover exact scope, missing data, opt-in/defaults, allele attribution,
      occurrence identity, pooled sources, and offline replay with tests.
- [x] Document CLI/API behavior; bump the version for the PR.
- [ ] Run lint, full tests, and equivalent CLI smoke; create PR and verify CI.
- [ ] Follow repository merge/release requirements, accounting for queued releases.

## Review

PR #809 includes the existing mhctools 4.9 mouse-allele test fix from
`fix/mhctools-49-mouse-allele-test`, with user-approved release version 3.51.2.
This preserves the queued docs 3.51.3 and packaging 3.51.4 sequence.

Lint and pip check passed. Strict docs build passed using the CI README-to-index
preparation. B16 direct smoke and native replay passed (9 candidates, 2 retained
queries); replay explicitly forbade Hitlist lookup and preserved per-allele scores.
New Hitlist cases cover viral/mutation antigens, all legacy prediction modes,
partial native annotation, default opt-in behavior, and unavailable replay.
Final full suite and GitHub CI are pending.
