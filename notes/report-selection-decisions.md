# Explain effective selection and product decisions (#730)

Reports must distinguish resolved design policy from parser/replay defaults.
Read the existing native selection record: run_configuration for epitope and
candidate settings, and construct_assembly for actual modality settings and
products. Expose this through an optional selection argument to the report API;
the CLI supplies dataset.selection. Missing historical fields stay unavailable.
Keep raw invocation arguments labeled as invocation provenance.

Show effective scoring/filter expressions, policy identity, RNA input values and
units, epitope/allele aggregation, window bounds/target retention, ranking/tie
rules and actual construct caps. Show each candidate's score expression with the
values it actually used, distinguishing ranking inputs from displayed context.

Attach concise selection explanations near actual products. Reuse recorded
window-selection policy, metrics, target details and alternatives, and existing
mRNA junction-swap metadata. Display nearest recorded alternatives and whether
self/serum terms were active; zero weight is not a completed assessment.
Recorded target masks and product coordinates remain explicit. Do not infer
unique source coordinates or invent historical tie-break decisions.

Capture missing decisions where they occur, without changing selection behavior:
window ordering/tie keys and selected positions in the existing window audit;
unique-sequence packing/capacity outcomes and complete source mappings in an
optional assembly audit. Embed assembly decisions in the existing source_graph
so its existing payload checksum covers them. Keep historical graph payloads
unchanged when no audit was recorded. Replay preserves the audit without running
selectors, predictors or assembly. No new sidecar or report-specific heuristic.

Reports explain exact-sequence deduplication and distinct overlapping products
without treating them as independent events. Full source links use existing
identity matching. Final capacity rejection and eligibility failure have
different recorded reasons; unavailable reasons are explicit.

Verify: saved cap10 versus parser default20; tied windows and recorded tie keys;
duplicate sources sharing one product; overlapping products from one event;
skipped processing/serum/self terms; capacity exclusion; junction optimization
and unassessed states; historical/empty data. Reports must not mutate native
data or products. Run actual Sid/COLO replay with inference/assembly forbidden,
inspect PDFs, then lint/full suite, review, green CI and clean-main deployment.

Implementation details: ordinary score-fraction/lexicographic selectors also
record their decisions. Per-junction candidate keys are recorded by the linker
optimizer, including the preferred and actually emitted linker after the
whole-construct comparison. Historical assembly records remain unchanged on
replay. Additive explanatory fields are excluded from assembly-input identity
but covered by the complete payload checksum; historical report decisions come
from those saved records rather than newly reranked candidate objects.
