# Expanded osteosarc read corpus

Build a reproducible local analysis dataset covering the 44 unique vaccine
loci represented by the earlier 49 source–variant regression cases, the
previously investigated complex variants, and additional promising coding
SNVs/indels from osteosarc.com. Preserve event identity separately from source
and processing-product identity.

Inventory public RNA libraries by platform, time point, tissue, provider and
processing stage. Retain all aligned records in broad event neighborhoods,
including reference/other alleles and complete stored read sequences, without
read-count subsampling. Include tumor and matched-normal DNA subsets for
compound-variant attribution. Large deletions and rearrangements require both
breakpoint neighborhoods. Record absent assays, unavailable alignments and
failed acquisition separately from measured zero coverage.

Cache source metadata, indexes, subset BAMs and hashes under a timestamped
local run directory. Reuse existing verified downloads where applicable.
Investigate upstream PacBio readsets and alignment provenance before treating
the currently mapped file as all available PacBio evidence. Keep alternate
processing products tied to their original biological library.

Acquisition policy (user decision, 2026-09-18): prioritize complete library and
timepoint coverage, not every reprocessing alternative. Pin reviewed source
groups and a representative per group; retain tagged long-read evidence and
DNA controls. Do not infer library equivalence from timepoint alone. Keep all
alternative URLs and provenance inventoried, preserve already acquired data,
and distinguish deliberately deferred products from failed or zero-coverage
queries. Apply the same selection to both event panels and acquisition retries.

Deliver a manifest, event inventory, per-source/per-event coverage summary,
acquisition commands and a readable limitations report. Keep large data out of
git. Validate subset integrity, coordinate systems, event coverage, deterministic
selection and absence of read-count caps. Run repository lint and tests for
maintainer-code changes. Translation and vaccine ranking are downstream of this
read acquisition task; coverage alone does not establish a coding product.
