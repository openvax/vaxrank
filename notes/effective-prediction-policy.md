# Effective prediction policy (#718)

Use Topiary's public exact-occurrence prediction API for both supplied candidate
and junction queries. Keep the existing mhctools model configuration boundary;
record the actual returned model/version/kind identities rather than inventing a
second predictor abstraction. Add an explicit junction rank endpoint (`auto`,
`pMHC_affinity`, `pMHC_presentation`); `auto` translates the historical
 affinity-first default. Require complete finite percentile-rank coverage for
one selected model/endpoint and the declared genotype. Missing evidence cannot
win an optimization as zero burden.

Persist one effective stage summary in the existing native run configuration,
and render it in the ordinary summary and reports. Separate candidate scoring,
junction optimization and allele-independent processing annotation. Show actual
models/versions, alleles, lengths, endpoint/thresholds, context, cache behavior,
coverage, and enabled/skipped/replayed status. Retain contextual processing
records in native evidence and replay them without inference; no new sidecar.

Use Topiary's public configured-model inventory for both candidate and junction
stages, separately from observed measurement identities. Require the release
fixing Topiary #541 so explicitly unknown model lengths remain reportable.
Keep Topiary selection_policy canonical while documenting the precise legacy
DSL/scalar translation and Vaxrank's separate RNA/window/assembly decisions.

Test mixed endpoints, missing rows/alleles/ranks, explicit model reuse, native
replay without constructing models, and processing occurrence round trips. Run
lint/full tests and real report smoke, review the effective summary, then ship
with a version bump after contextual caching.
