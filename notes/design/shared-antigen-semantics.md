# Shared antigen semantics (3.59.0)

LENS and Isovar must express the same biological facts in the same fields.
Input format is provenance, not a biological event type or an admission policy.
Keep sequence novelty, expression-defined targets, RNA products, coding context,
and supplied evidence distinct; do not infer a splice-site DNA mutation from an
alternative transcript or a fusion from a breakend alone.

1. Reproduce and fix the three 3.58.0 review findings: direct SV prediction and
   reference comparators, selected-variant summaries on initial runs/replay,
   and frozen CTA references in admission-only datasets.
2. Add a small shared, serializable annotation record for source events and
   coding contexts. Populate it at the direct, LENS and Isovar boundaries;
   retain original producer evidence separately. Shared consumers use these
   fields for source classification, selectors, tables, biological reports and
   product provenance. Unknown fields stay unknown.
3. Expose `variant` as the parent of small and structural variants. Add explicit
   SNV/MNV/indel/SV-subtype and variant-coding-context selection, independent of
   producer. Preserve the historical `mutation`, `SV` and `ncORF` selection
   spellings for compatibility and document their meanings. Normalize CLI and
   configuration selectors consistently. RNA `splice`/`isoform` labels describe
   products and never manufacture evidence of a causal DNA mutation.
4. Share target-overlap/source-class logic and selected-source summaries across
   paths. Keep reconstruction, normal/reference comparisons and targetable
   intervals evidence-driven; this PR does not invent missing splice junctions,
   full ORFs or tumor/normal comparisons.
5. Test actual prediction through the SV path, matching LENS/Isovar semantics,
   subtype/context selection in both modalities, direct/replay summary parity,
   empty CTA policy replay, serialization and historical replay. Run lint, the
   full suite, CLI/report smokes and strict documentation checks before shipping.

Use the existing isolated environment and preserve unrelated local work. Open
one new PR with a version bump and the filed defects. Follow normal CI, merge
and deployment gates.

Primary semantic contracts: VCF 4.5 defines symbolic SV types and oriented
breakends; Ensembl/Sequence Ontology distinguishes sequence changes from their
transcript consequences (including splice effects). Local Varcode and Isovar
public models own event reconstruction and coordinate interpretation.
