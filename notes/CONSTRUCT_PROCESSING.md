# #422 follow-up: context and patient-MHC overlays

Build on the full per-bond profiles and report in PR #433. Do not change ranking.
Do not close #422 until native, actual final and complete translated products
are covered; #423 remains responsible for independent Sid data and real caches.
Include the reproduced legacy endpoint-sentinel bug in #434: retain internal
evidence, report endpoint/unassessed C-boundary status, and never turn a forced
endpoint zero into a measured probability or a composite processing deficit.

## Invariants

- Preserve the exact `ConstructSequence` and every `ConstructPlacement`, with
  native source objects, gene/transcript/species, evidence, edits and chemistry.
  Complete mRNA inputs must be checked against actual CDS translation and emitted
  nucleotide views. Linkers and encoded elements are included, not cropped away.
- Native RNA/protein contexts may be windows, not complete proteins or exposed
  peptide molecules. Report that scope. Pepsickle can provide explicitly padded
  sequence-context scores, but terminal peptidase rules must not treat artificial
  window ends as exposed molecular termini. Native terminal trimming stays
  unassessed without explicit molecular/terminal context.
- Final chemistry and model-input assumptions must be separate from documented
  manufacture. Unsupported or unresolved chemistry is never silently replaced
  with an unmodified bare sequence. Conditional assumptions, if supplied, must
  remain visibly conditional in JSON and the readable report.
- Score distinct contexts once; retain repeated occurrences and all source
  placements. Do not reuse native-context scores for an edited final sequence.
- Predict unfiltered patient-MHC output through Topiary and retain every returned
  peptide/offset/model/kind/allele. Reuse existing complete-window/self assessment
  where a source antigen exists; do not invent one gene or antigen kind for a
  multi-antigen product merely to satisfy an API shape.
- Overlay all returned ligand occurrences and target masks on exact bonds.
  Show internal and boundary observations separately, including zero-width target
  junctions. Sequence preservation, target overlap, cleavage-model evidence and
  historical selection agreement are different observations.
- Missing/unreturned HLA/model/length combinations stay explicit. Binding or a
  cleavage score does not establish presentation, T-cell recognition, intact
  delivery, or clinical safety. All report text is source-agnostic (target antigen,
  not always mutant), and escaped.
- Keep serum/extracellular, proteasomal, other cytosolic, ER and endolysosomal
  evidence separate. An enzyme having one compartment annotation is not matrix
  calibration or a claim that the molecule reaches that compartment.

## Suggested delivery

Add a typed complete-context audit around the existing immutable records, with
explicit native/final/product identity, processing profiles, prediction coverage
and ligand/target overlays. Extend the report using small glue and templates.
Avoid a second hand-maintained native serializer or an ad-hoc ranking filter.

Tests: repeated motifs with different flanks, both ends, truncation versus real
termini, additions/deletions and target junctions, mutation/CTA/viral masks,
multi-antigen product provenance, full translated CDS validation, empty/error/
partial outputs, chemistry, native JSON and offline report regressions. Run
lint, full tests with coverage, real-model/report smoke, CI and release gates.

## Public API

`native_sequence_context(antigen)` retains the entire supplied native context;
it does not assert that reconstruction produced an entire protein. Use
`final_sequence_context(construct)` for an exact manufactured antigen, and
`translated_sequence_context(rna_product, evidence)` for the complete translated
CDS, including linkers, signal peptides and other encoded elements. The latter
retains the entire `RNAConstruct` graph through the native serializer and checks
its CDS translation, all placements, and emitted nucleotide views. It snapshots
mutable assembly data and revalidates it before inference/report export.

```python
from mhctools import NetMHCpan42
from mhctools.peptidases import get_cleavage_model
from vaxrank import (
    MHCRequest, native_sequence_context, final_sequence_context,
    audit_sequence_contexts, write_sequence_context_audits,
)

# `construct` is an attributed ConstructSequence, not a guessed manufactured
# sequence. Supply actual patient typing and the installed model's output
# identity/version in `requests`. Versions must match output, not assumptions.
contexts = [native_sequence_context(construct.native_antigen),
            final_sequence_context(construct)]
# Defaults preserve unknown chemistry. An explicitly conditional analysis may
# use n_term='free', c_term='free', chemistry_basis='Conditional ...' instead.
# Documented arbitrary chemistry is still unsupported; this cannot override it.
audits = audit_sequence_contexts(
    contexts, mhc_predictor=NetMHCpan42(alleles=patient_alleles,
                                     default_peptide_lengths=[8, 9, 10, 11]),
    mhc_requests=requests, pepsickle=True,
    peptidase_predictors=[get_cleavage_model('cpn-basic')])
write_sequence_context_audits(audits, json_path='contexts.json',
                             html_path='contexts.html')
```

Each `MHCRequest(kind, predictor_name, predictor_version, allele,
peptide_length)` describes one requested output combination. Coverage lists all
missing occurrence offsets, including completely unreturned combinations.
Unrequested outputs are retained and flagged, not silently discarded. The
request set is explicit input provenance, not auto-inferred from successful
rows; incomplete requests must not be called complete patient coverage.

Native windows and unknown-terminus translated products can receive Pepsickle
sequence-context scores, labeled `sequence_context_only`; their model's padding
does not establish exposed termini. Terminal peptidase rules do not run on native
windows or embedded encoded antigen segments. Unknown terminal chemistry stays
unknown for enzyme inference. Enzyme compartment metadata is not serum-matrix
calibration. No ranking, manufacture or target-selection decisions are changed.

The JSON contains native `audits` plus derived `coverage` and `overlays` entries.
Every target mask and returned ligand has exact internal/boundary observations;
endpoints remain `sequence_endpoint`. The HTML lists all model/site observations,
all ligand occurrences, requested/unreturned combinations, source graphs,
chemistry assumptions and limitations. MHC values use mhctools' native units.

This API does not identify independently intended target ligands, adjudicate
historical sequence choices, or apply self/CTA/tissue policy. These require the
independent #423 audit and #414 selection validation. Existing
`audit_construct_sequence` remains the attributed single-antigen self-reference
entry point; no fake gene or antigen kind is created for a multi-antigen product.

## Legacy endpoint correction (#434)

The ordinary HTML/ASCII processing columns now print `sequence endpoint` at the
end of supplied context, retain maximum internal-cut evidence, and leave the
composite unavailable. Native JSON carries the explicit boundary status and
null missing scores. A real zero at an internal bond remains zero. Exact-length,
finite probability validation prevents malformed arrays from shifting the
apparent endpoint. This does not change vaccine ranking.

Primary scope: [Pepsickle documentation](https://github.com/pdxgx/pepsickle)
describes eight-residue contexts and padding; it does not turn absent context
into molecular-terminal evidence. The released mhctools enzyme models retain
their primary references, assay descriptions and limitations in every profile.
