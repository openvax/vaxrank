# Construct provenance and Sid validation

## Public API

`ConstructSequence` records one native antigen window and its final sequence.
`ConstructSequenceEdit` uses native antigen coordinates; a
`ConstructChemicalModification` uses final amino-acid coordinates. The public
`ConstructEvidence` record independently attributes sequence observations and
modification rationales. The native antigen is never rewritten in place.

Use `save_construct_sequences` / `load_construct_sequences` for lossless,
allowlisted JSON input. `audit_construct_sequence(record, predictor,
genome=genome)` sends the final sequence through Topiary and preserves all
non-excluded exact-self source genes. With no predictor, unresolved native
mapping, or unsupported chemical changes it records **unassessed**, not safe.
`write_construct_audits` emits native JSON and an optional escaped HTML report.

For explicitly supplied peptide products, use
`vaxrank.peptide.peptide_construct_from_sequence(record)` followed by
`write_peptide_outputs`. A provenance-bearing product requires a JSON manifest
sidecar, and its vendor-order CSV also includes the complete provenance record.
Writer options cannot silently add chemical modifications to such a record.

`ConstructPlacement(record, offset)` identifies an exact occurrence in an
assembled peptide or encoded mRNA product. Attach these through the product's
`construct_placements` field. Writers check the actual emitted amino-acid
context, and mRNA output also verifies CDS translation. These records do not
change automatic vaccine selection or manufacture an unrequested tail.

MHC inventory status `predictions_returned` means exactly that; it does not
promise all expected alleles/models/lengths were assessed. The report retains
all returned ligands instead of filtering to strong targets. Processing maps,
serum stability, and historical selection agreement are separate assessments.

The published DYNC1H1 sequence and native context are documented at
https://osteosarc.com/variant/DYNC1H1-chr14-101980529/. Solubility rationale is
attributed separately to the project owner's 2026-09-09 report of JLF practice;
CSBio involvement remains tentative.

## Delivery plan

Work sequence: #421 (native/final sequence representation and reassessment),
#414 (independent historical selection comparisons), #422 (complete cleavage
maps), #423 (documented construct audit). #355 is an independent input-validation
PR. Each release requires review, lint, the full suite, relevant output smoke,
green CI and verified deployment. Partial foundations do not close their parent
issues.

## #421 contract

- An immutable construct sequence records the native antigen, an explicit native
  window, the final amino-acid sequence, modality, and attributed edits. Never
  infer a native occurrence by choosing the first matching substring.
- Sequence edits use zero-based half-open coordinates in the native antigen;
  chemical modifications use final-construct coordinates and do not pretend to
  change the encoded amino-acid sequence. An unknown mapping stays unresolved.
- Insertion/substitution residues have no native/RNA coordinate. Targets map
  only through retained native residues or preserved native deletion junctions.
  Removing a target must remain visible, not create a new target from a tail.
- Sequence evidence and modification rationale have separate provenance. The
  DYNC1H1 page documents JLF V2/V3's final sequence. The user reports that JLF
  uses terminal lysines for solubility; CSBio involvement is tentative. Neither
  statement authorizes automatic tail design or identifies other residue types.
- Round-trip the complete public object graph through the allowlisted native
  serializer. Identity includes chemical form, source occurrence and provenance;
  prediction inputs use the actual final context, not cached native context.
- Preserve these facts in manufacturing exports and readable reports. Reassess
  final contexts through Topiary and the existing antigen-aware self-reference
  machinery. Unsupported chemistry/mapping and absent processing models remain
  explicitly unassessed. No automatic sequence modification or ranking change.

## Following PRs

### #414 pilot specification

Import five original-read cases and the GRCh38 reference from the immutable
Isovar 1.8.1 commit. Keep both DYNC1H1 loci distinct, ONT T1/T2 distinct, and
MAP2's insufficient evidence explicit. Independently record DYNC1H1, EXOC4
and H1-2 provider sequences. Run real NetMHCpan 4.2 predictions once, pin their
inputs and outputs, and exercise offline Topiary cache -> filtering -> ranking
for multiple competing windows. The clinical typing table, null HLA allele,
class-I-only scope, incomplete reference proteome, and unknown historical
selection settings must be visible. Compare exact sequence and target offsets;
do not call a present-day output snapshot historical agreement. Native peptide
lengths are compared separately from documented manufacturing additions.

Draft-review correction: record modality explicitly on each independent source
record; provider/display labels must not decide whether a construct is mRNA or
a manufactured peptide. Regress the mRNA minimal-epitope case. Regenerate real
prediction artifacts into a new directory and compare their observations before
adopting updated metadata; do not change RNA gates or historical expectations.
This local preparation does not authorize publication of genomic fixtures.

Topiary #296, the cached-coordinate blocker, is fixed upstream by Topiary #299
and released as Topiary 5.55.1 on PyPI (2026-09-11); the topiary floor is
raised to `>=5.55.1` accordingly. Topiary #300, a separate gap this pilot
surfaced (`serum_half_life`/`blood_half_life` declared but unreachable from
the ranking DSL), is fixed in the same release. Reaching that fix also
required raising vaxrank's own mhctools floor to `>=3.39.0` to match
topiary's real requirement; below that floor mhctools' `Kind` class predates
both kinds, so topiary's fix was reachable but not exercised. All 18 pilot
tests and the full suite pass against the real published release, verified
directly against PyPI's index rather than the topiary source tree.

- #414: pin Isovar >=1.8.1 fixtures with their original sample/timepoint metadata;
  independently snapshot documented provider-specific sequences; compare actual
  reconstruction, candidate generation and final selection. Missing historical
  settings and disagreements are explicit, not fabricated matching predictions.
- #422: immutable full-context site scores with exact bond convention, model and
  compartment provenance; retain repeated occurrences; overlay target ligands;
  test endpoints, invalid arrays, changed contexts and missing coverage.
- #423: real, versioned prediction caches and source-aware non-CTA matches for a
  documented Sid pilot, plus an explicit inventory of remaining unassessed data.
- #355: reject malformed allele-scoped evidence at input with actionable source
  location, while report rendering can represent unavailable display evidence.

## Verification

Independent DYNC1H1 native/final expectations; terminal and internal insertions,
deletions, substitutions, chemical modifications, ambiguous/unresolved mappings,
CTA/viral targets, retained deletion junctions, exact JSON round-trips, cache
identity, new boundary ligands, and non-CTA shared-source preservation. Check
actual emitted files, not only helper return values. Each parent issue remains
open until its full acceptance criteria are verified.
