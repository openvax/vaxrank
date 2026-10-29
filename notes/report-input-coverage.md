# Input coverage and attrition in reports (#728)

Lead with input coverage, then final products and selection reasons, grouped
supporting evidence, and methods/provenance. Count declared inputs, resolved
antigen records, admitted records, ranked sources and represented source records
separately. Label events, reference models, RNA reconstructions, peptide
occurrences and unique products explicitly; they are not interchangeable.

Preserve upstream VCF/BAM counts and per-variant outcomes at first capture, using
the existing native dataset/source catalog. Keep original observations across
pooling and replay. Historical missing denominators remain unknown. Retain
producer-unresolved ORF models and their reasons, plus held-out candidates and
CTA expression admission decisions with their original units.

Give every unrepresented source a readable identity and stage: admission,
missing recorded prediction, no ranked window, or ranked but absent from actual
products. Explain recorded/verified observations only. Retained predictions do
not establish complete inference coverage. Reasons can overlap; do not add them
as disjoint losses. #730 supplies finer product-cap/optimization rationale.

Show historical input provenance and current replay inputs separately, including
manifest declarations, patient/sample/library/timepoint and reference scope.
Record regional/whole-library scope only when explicit; never infer it from a
filename. Absent CTA/viral inputs mean no supplied evidence, not screening or
exclusion. Introduce no sidecar or admission/ranking/product behavior change.

Regression: mixed direct/native replay, unknown historical counts, absent
families, unresolved ORF models, held-out fusion/SV, missing/filtered predictions
and product-cap rejection. Verify HTML/PDF/text order and product parity, then
lint/full tests, PR review, green CI and clean-main deployment. #717, #719, #730
and #736 remain separate followups within the active goal.

Implementation note: use the existing ExternalRankingOutcome entries when building
the current-stage ledger. They already describe parsed LENS/pVACseq/native
inputs, including sources with no ranked windows. Do not duplicate those
adapters by reconstructing source identity from surviving peptide rows. Source
catalogs supply upstream counts/producer models and admission evidence; current
ranking outcomes and actual product source mappings supply later stages. Keep
unknown legacy denominators separate from measured zeroes.

Direct capture marks deferred window selection as not assessed, never zero.
New optional per-input rna_evidence_extent declarations preserve regional,
whole-library or targeted scope without guessing from paths or changing old
input identities when absent. External adapters also expose their existing
source-policy and parsing exclusions, preserving the same decisions on replay.
