# Admission and source evidence in reports (#729)

Report construct admission as the recorded policy decision. Keep its evidence
type, rationale, overrides and review flag explicit; a cleared admission flag
does not record biological validation or completed scientific review.

Use the existing native antigen/source-evidence records to show sample support,
coding extent, ambiguity and effective admission controls. Name every count's
subject and unit: source RNA, coding intervals, ORF occurrences and peptide
occurrences are different measurements. Preserve unavailable values as unknown.
Catalog/reference evidence must be separate from evidence measured in the sample.
Render every retained evidence record and origin; do not choose one conflicting
origin as authoritative. Do not introduce a new output sidecar or alter selection.

Keep small source-evidence view helpers outside the main formatter. Feed them
the existing linked source catalog and render the same sections in HTML/PDF and
text. New Isovar adaptations retain explicit partial/ambiguous candidate facts;
historical records show only facts recoverable from their recorded evidence.

Verify canonical and generic reports, research-only assumptions, reference versus
sample evidence, missing support, source versus exact-occurrence counts, and
native replay with inference forbidden. Retain the #731 processing-replay check.
Run lint, full tests, real Sid/COLO report smoke and visual inspection, then review
the PR, wait for CI, merge and deploy with a version bump. Coverage/attrition (#728),
ORF labels (#717), decision rationale (#730) and grouping (#719) remain distinct
follow-up requirements in the active goal.
