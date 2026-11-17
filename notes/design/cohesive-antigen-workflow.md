# Cohesive antigen workflow follow-up

The common input, prediction, selection, assembly and report paths must preserve
each source's measurements and reference identity. Source switches must give the
same eligible vaccine sequences regardless of duplicate-source input order.

1. Fix #851: choose enabled source occurrences before sequence deduplication,
   retaining all origins and counting each emitted sequence once. Exercise both
   peptide and mRNA products and inspect source limits for the same problem.
2. Fix #852: give independently measured CTA transcript occurrences distinct
   identities; deduplicate their final products only after scoring. Preserve
   expression, admission status and equivalent protein occurrences on replay.
3. Fix #853: retain host-inventory matches independently of CTA exclusion-policy
   matches. Validate both frozen reference contracts without querying live data
   or weakening native host-reference checks.
4. Audit category selection, mixed fresh/saved inputs, reference-policy resolution,
   source evidence, biological reports and native replay together. Add focused
   regressions for concrete gaps; file any newly confirmed defects.
5. Reconcile the documentation and packaging PR feedback (#784, #801, #846),
   preserving newer CLI and antigen docs, contributor attribution and explicit
   build/test dependencies. Verify their current heads before changing them.

Use the existing isolated environment and preserve unrelated local work. Run
focused tests while iterating, then lint, full tests, mixed-source output/report
smokes and visual report checks. Ship through reviewed PRs with green CI, normal
deployment from clean main and a verified PyPI release.

Scientific contract: OncoRef's primary `docs/sequence-bundle.md` distinguishes
raw host inventory from CTA policy roles; all matching occurrences aggregate as
exclude > review > allow. Independent transcript expression measurements remain
independent even when the encoded amino-acid sequence is identical.

The audit also found direct Varcode StructuralVariant fragments being labeled
as ordinary mutations. Preserve the producer's SV class and event endpoints
through source switches, native evidence, both writers and biological reports.
The primary contract is Varcode's `structural_variant.py`; this change does not
infer translations or gene fusions from an SV type.
