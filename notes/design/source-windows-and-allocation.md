# Distinct source windows and mixed-antigen allocation

## Contract

New designs retain at most three non-overlapping target windows per antigen
source, configurable with `--max-windows-per-antigen-source` and
`vaccine_peptides.max_windows_per_antigen_source`. Coordinates are zero-based,
half-open on the actual translated source sequence. Touching boundaries are
allowed; any shared residue is overlap. A source is one recorded mutation/event
or one explicitly identified translated antigen, not a gene-name spelling and
not a specimen. Competing source observations remain separate evidence; use
the existing best-observation selection rather than summing them.

Choose the best available window with the configured score and window policy,
remove overlapping windows, then repeat on the remaining regions. Apply the
near-best comparison anew within each round, so a strong first region cannot
erase a useful second region. Repeated copies with no additional target
peptide/allele content do not consume another source window. This is a
deterministic greedy policy, not a global optimum or an immunogenicity model.

Manufacturing alternatives belong to a selected region and must retain every
target peptide/allele occurrence of that region. They do not create additional
target windows or vaccine slots. Retain coordinates and required content in
the audit so later shift/trim selection can enforce that contract. Refinement
must preserve coordinates on the original parent and non-overlap across final
windows. Existing unversioned legacy windows remain replayable; saved runs
without this setting and the frozen `openvax-v1` profile retain legacy behavior.

## Category diversity

After eligibility checks and exact emitted-sequence deduplication, allocate
unique antigen sequences using a shared `antigen_allocation` policy for peptide
and mRNA. `balanced` (the fresh-design default) rotates through the least
represented available antigen kinds, breaking ties by the existing ranked /
allele-coverage order. `ranked` keeps that order. Optional per-kind minimum
representation goals are served first; optional maxima are hard limits.
Unused capacity passes to other eligible kinds. Never admit a held-out antigen
to meet a goal. Count each unique sequence once, attributed to its first
(highest-priority) source kind; preserve all other source mappings. ORF context
is provenance, not a second allocation category.

Expose the policy in Python, YAML and CLI, record effective settings and
per-product decisions, and report actual category counts and unmet goals after
packing. Diversity is a configurable allocation preference, not a change to
the Topiary DSL scoring expressions or evidence units. Published multivalent
vaccine studies motivate evaluating breadth but do not establish an optimal
CTA/mutation/viral ratio: no such scientific claim is made here.

Primary references checked: Sahin et al., Nature 2017,
https://doi.org/10.1038/nature23003; Saxena et al., Cancer Discovery 2025,
https://doi.org/10.1158/2159-8290.CD-24-0934.

## Verification

Test disjoint versus adjacent/overlapping windows, repeated sequences,
manufacturing target preservation, native/direct/table inputs, two observations
of a source, CLI/YAML equivalence, frozen legacy replay, duplicate mixed-kind
products, caps/goals with missing kinds and limited capacity, and both assembly
modalities. Run lint, full tests, and a mixed-antigen CLI smoke with native replay.
Related: https://github.com/openvax/vaxrank/issues/758.

Existing packing limitation found during implementation: an individual antigen
larger than the configured length cap is still emitted with a warning. Tracked
separately in https://github.com/openvax/vaxrank/issues/826; this change preserves
that existing behavior. Category maxima and construct-count limits are enforced.

Scoring interaction fixed with this work: window utility previously omitted
captured MS increments (#828). It now preserves the existing maximum increment
per peptide sequence alongside the unique peptide/allele prediction terms. An
intact duplicate copy preserves that increment during cleavage assessment.
