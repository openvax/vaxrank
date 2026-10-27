# Isovar ORF inputs

## Outcome

Add one documented Isovar input path to the existing native antigen workflow.
It accepts released producer results for offline replay and source requests that
call Isovar directly. Small-variant/ncORF, SV, fusion, viral and targeted-isoform
sequences enter the same prediction, selection, reporting and serialization
pipeline. Ordinary VCF/BAM behavior remains available through its existing API.

## Producer responsibilities

Require Isovar 1.61 and its compatible ORFdata release. Isovar owns RNA assembly,
frame matching, evidence accounting and translation. Vaxrank adapts declared
producer schemas, checks sequence/coordinate/reference consistency, and retains
the complete source evidence and unresolved outcomes. It does not reconstruct
missing flanks, infer tumor specificity from RNA, or select a strongest frame
before the shared ranking stage.

Supported source families are `variant_orfs`, `sv`, `fusion`, `viral`,
`isoforms`, and saved producer results. Small-variant ORF requests use validated
ORFdata transcript models; SV requests use Isovar's nominated-event contract;
viral requests require the alignment reference and explicit source scope;
isoform requests retain annotated competitors and distinguish reference-only
proteins from observed sequence. Existing protein-hypothesis exports should
also enter through this path where their versioned evidence supplies the
necessary sequence and mutation coordinates.

## Shared API and controls

Expose a typed reconstruction/admission configuration and an Isovar source
adapter that returns an `AntigenDataset` with retained evidence. Add `isovar`
to existing input manifests and a matching `--input-isovar` convenience option.
Use the existing YAML and `--set` mechanisms for the same configuration names;
reject unknown or inapplicable controls instead of ignoring them. Record the
fully resolved source options and actual producer version in native outputs.

Separate reconstruction controls (read evidence, assembly, reference frames,
source-specific resource limits) from downstream admission controls (partial
coding regions, ambiguous models, source-specific attestations). Defaults retain
uncertainty and all candidates. Missing attestations and insufficient evidence
remain held out with explicit reasons. Required reference/input mismatches are
errors; absence of reconstructed proteins is a retained scientific outcome.

## Evidence and identity

Store source results once per input, with content identity and candidate links,
rather than repeating large read inventories on every ligand. Preserve sample,
library, reference/build, genetic code, model/ORF IDs, gene/transcript/protein
attributions, source species, sequence extent and support units. Keep tiling,
single-read spans, fragment support, abundance estimates and reference evidence
distinct. Saved antigen and prediction exports must retain unresolved inputs
and replay without running Isovar or consulting the current annotation.

Targetability follows explicit altered intervals or junctions for mutation/SV
inputs, and reviewed source-specific admission for viral/isoform inputs.
Neither an arbitrary longest ORF nor an entire host protein becomes targetable
merely because a coding prediction exists. Handle termination-only events,
synonymous ncORF results, alternate frames and partial ends explicitly.

## Validation and delivery

Use actual released Isovar outputs in bounded offline tests, including supplied
ncORF frames, SV alternatives, viral sample/reference distinctions, targeted
isoforms and no-candidate cases. Test API/CLI/YAML equivalence, invalid options,
reference mismatches, nonstandard genetic codes, all-held-out runs, mixed inputs,
native replay and ordinary VCF/BAM compatibility. Exercise shared ranking and
report generation, not just adapter object shapes. Add a current support matrix,
configuration reference and runnable examples, and correct stale Isovar docs.

Run lint, full tests and a CLI/report smoke test; bump the package version in the
PR. File any missing upstream contract on its owner repository, link it, and
state the supported boundary precisely. Follow repository CI/merge/release gates.
