# Inspect saved design reports without redundant evidence export (#736)

Add an explicit --no-output-epitopes CLI control that suppresses the default
native prediction-table export under --output-dir. Keep --output-epitopes PATH
mutually exclusive with the disabling control. Default
output bundles retain their existing behavior; explicit report and product
outputs are written normally. The option affects export only, not native input
loading, source evidence, effective policies or replay eligibility.

Exercise saved peptide/mRNA products through HTML/PDF/text reporting with the
large native export disabled and every predictor/reconstructor/assembler
forbidden. Verify exact products, source mappings, policy and original evidence
are preserved and no new native evidence table is written. Keep source input
paths explicit. Document the CLI and existing report API path. No extra sidecar.

This PR follows #740; version and release sequencing must be checked against
current main before merge. Run lint, full tests, real saved-report smoke, review,
all CI, then merge and immediately deploy from clean main.
