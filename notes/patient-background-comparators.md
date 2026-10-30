# Exact patient-background comparisons (#716)

Keep `annotation_reference` as the default comparison policy. Add an opt-in
`patient_background` policy with separately declared `focal` or
`supplied_somatic_haplotype` reversion. Use Isovar's released exports and
Topiary's validating reader; do not infer a pair from peptide similarity.

The source export lives once in the existing native source-evidence store.
Candidate evidence refers to its exact event, hypothesis, translation,
transcript and pair IDs and records the peptide's original coordinate interval.
Match against the actual selected Isovar hypothesis and every compatible
nucleotide/transcript alternative, checking sample/read-set provenance. Generic
Isovar protein hypotheses use the same contract. Other source families have no
implied patient-background pair.

Project peptide boundaries through the recorded reverted cDNA edits. Keep
germline and unexplained changes. A boundary inside a reverted edit, a changed
reading frame at a boundary, a stop before the required interval, or a missing
alternative gives an explicit unavailable comparison. Reconcile every matched
alternative before reporting a single sequence; differing sequences or mixed
availability remain ambiguous. Unequal lengths are explicit and must be within
the configured predictor's supported lengths before prediction.

Use the existing `comparators` map for predicted peptides, keeping
`patient_background` distinct from the legacy annotation `wt` role. Exact
bounded equality, absent sequence and absent prediction remain separate states.
Do not drop an identical comparison or replace an unavailable patient comparison
with annotation. Predictions go through Topiary. Default ranking remains the
existing mutant-evidence DSL; comparison evidence is explanatory unless a policy
explicitly uses its supported fields.

The API, YAML and CLI must resolve the same effective policy. Live RNA runs
export evidence from the already reconstructed results, with explicit sample
and read-set identities. Saved comparisons replay without source files or new
inference. Reports show the selected comparison, alternative counts, phase,
construction versus original REF-allele RNA support, and evidence links. A
computed baseline is not observed normal-tissue or protein expression.

Validate real Sid NTF3 provenance, retained germline indels, in-frame and
frameshifting focal edits, multiple transcripts, capped/partial contexts,
missing phase and inconsistent provenance. Verify native replay, CLI/config,
prediction routing and reports; then lint/full tests/CI and ordered release.

Consumer integration found valid unlimited caps and repeated pair records that
Topiary rejected. Topiary #563 fixes #561/#562; this change requires its 5.108.4
release. Report review also reproduced Vaxrank #753: conflicting comparator
affinities were reduced to the first leaf. Human-readable and spreadsheet views
now share exact model/version matching and explicitly preserve conflicts.

Scientific context: Rubinsteyn et al., *Computational Pipeline for the PGV-001
Neoantigen Vaccine Trial*, doi:10.3389/fimmu.2017.01807; Hundal et al.,
*Accounting for proximal variants improves neoantigen prediction*,
doi:10.1038/s41588-018-0283-9. The executable exchange contract is Isovar's
versioned patient-background-pairs documentation and Topiary's validating
`read_isovar_patient_background_pairs` reader.

## Current OncoRef integration

Integrate the released 3.48.0 common CTA admission path and independent 3.48.2
coverage reports before shipping this PR as 3.48.3. Require the reviewed
OncoRef `>=1.23.1,<=1.24.0` range in the ordinary install; the `cta` extra remains
a compatibility alias. PirlyGenes coordinates its supported range so current
Topiary and common CTA capabilities resolve together. Do not lower the Topiary
floor or bypass dependency resolution.

Check the actual common human/dog/mouse snapshot and native replay workflows,
including normal-tissue counterevidence, incomplete or pooled RNA evidence,
full-protein equivalent occurrences and missing host coverage. Preserve the
historical default admission policy and every saved evidence identity. An
owner package upgrade is separate from adopting the reviewed target policy
tracked in #687. Re-run the combined full suite, real report replay, docs,
joint installation and current-head CI before normal clean-main deployment.

Before shipping, test both OncoRef endpoints with PirlyGenes 6.1.0 through
ordinary dependency resolution. That release supports the same bounded range,
so it can publish against the previous Vaxrank before this release enables the
new endpoint. Include current Tsarina 1.36.0: PirlyGenes #664 exposed a mix of
historical and current CTA membership and a null-unsafe comparison in its tests.
Preserve the explicitly historical policy while updating dependency contracts.

The 1.24 source comparison changes the reviewed CTA-selection policy and adds
TCR-T product/cohort/event evidence; it leaves the historical CTA API and common
species evidence/reference assets unchanged. Clinical credit belongs to that
explicit owner policy and must not become implicit Vaxrank admission or ranking.
Retain existing snapshot identities and rerun the public-data, common admission,
native replay and full regression suite on the new endpoint. The scope of
clinical evidence was checked against the primary
[afami-cel](https://doi.org/10.1038/s41591-022-02128-z) and
[IMA203](https://doi.org/10.1038/s41591-025-03650-6) reports: those studies evaluate
particular TCR-T products in selected patients, not arbitrary vaccine windows.
