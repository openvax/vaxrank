# Generalized epitope inputs: first consumer milestone

Keep ordinary VCF/BAM and single-report commands unchanged. Add
`--input-epitopes FILE` for native reload and `--input-topiary FILE` for normalized
Topiary tables; both also participate in the existing input manifest and
repeatable input options. No new scoring language or source-category branches.

Use one EpitopeDataset containing the Topiary evidence result, native
CandidateEpitope objects, declared input scope, scoring configuration, and
explicit candidate-to-antigen references. Topiary owns combination identities,
representative selection and typed table persistence. Candidate and antigen
payloads use Vaxrank's existing native serializer. A separate scoring view maps
source observation IDs onto Vaxrank prediction IDs without modifying historical
columns. Original predictions are the default; loading never requests inference.

A dataset export keeps all observations and annotations, unknown versus zero,
unknown versus empty flanks, metadata and the selected scoring policy. Native
candidate-only files remain readable. Tables lacking construction evidence remain
rankable and visible with a reason. Explicit admitted antigens use the existing
VaccinePeptide and modality consumers; importing a category does not admit it.
Existing report-specific construction adapters remain until #497's construction
milestone replaces them; this PR does not reconstruct contexts or combine direct
VCF/BAM runs with reports.

Acceptance: simple commands retain their defaults; native and enriched table
save/reload preserve predictions, provenance and DSL scores; typed annotations and
additive features survive; compatible mixed inputs share loading/scoring while
incompatible patient/reference scope fails before scoring; rediscovery cannot
multiply selected candidates; original-only input does not initialize predictors.
Regression tests exercise file/CLI round trips, old native files and unchanged
report loaders. Run lint, the full test script and an equivalent CLI smoke before
opening the version-bumped PR, then CI, merge and deployment.
