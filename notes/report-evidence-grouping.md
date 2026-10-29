# Group report evidence without merging biological observations (#719)

Baseline: the saved COLO report rendered by #746 has 402 pages, 112 candidate
cards and 896 candidate tables. A read-only HTML comparison found 342 distinct
table texts and one displayed target-epitope table repeated 112 times. Equal
formatted text does not prove equal underlying evidence: grouping must also
respect unrounded predictions and their context.

Use short deterministic report source identifiers with full identities in one
linked source directory. Keep every source and its original provenance; these
are display identifiers, never replacements for native keys or product names.

Group source-outcome rows only when their family, gene/coding context, origin,
stage/reasons, admission, retained-record count and actual product outcome agree.
Show all member source IDs and ranking positions. Unknowns remain unknown and
counts remain source-record counts, not independent events.

For candidate details, share repeated evidence components rather than asserting
that every source of the same peptide is equivalent. Full window groups must
respect sequence, target mask, prediction and processing context, unrounded
measurements, comparator evidence, source support and admission. Where sources
share an exact table but differ elsewhere, reference the shared table and keep
the differing fields visible in per-source context/support rows. Include full
linked source memberships. Reference/catalog evidence and shared input/admission
controls can be rendered once, while source-specific RNA support remains
attributed to its original candidate.

Keep this as a report view assembled by small Python helpers and Jinja
templates. Do not change ranking, filters, native evidence, candidate CSV/XLSX,
product manifests or the underlying source catalogs. No new output sidecar.
HTML and PDF links must resolve; PDF readers must be able to inspect the full
identities and differing fields without hidden HTML-only disclosure controls.
Avoid repeating the same full hash in every table.

Tests: identical windows and evidence share a detail block; changing flanks,
target masks, support, raw predictions that round equally, WT comparators or
admission remains visibly distinct. Preserve annotated/noncanonical contexts,
multiple genes/ORFs, held-out sources, historical unknowns, empty reports and
all native/CSV rows. Verify complete source links and exact product parity.
Run Sid and COLO saved-evidence reports with inference/assembly forbidden,
inspect PDF pages and compare table/page counts, then lint/full suite, review,
green CI and clean-main release. Counts alone are not the usability criterion.
