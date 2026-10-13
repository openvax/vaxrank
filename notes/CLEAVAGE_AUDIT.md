# #422: complete sequence processing evidence

This is a non-mutating audit, not a ranking policy or a clinical verdict.

1. Preserve the exact native antigen context, final manufactured sequence, or
   complete translated product. Reuse `ConstructSequence`/`ConstructPlacement`
   for edits, target masks, repeated occurrences, chemistry and source evidence;
   never infer an occurrence by choosing the first sequence match.
2. Preserve immutable per-bond results using mhctools' published cleavage
   contract. Keep native scores and motif recognition distinct. The Pepsickle
   adapter records its complete raw vector, including the final residue output
   separately from actual internal bonds. Record backend/version/model/settings,
   model-asset identity where available, compartment and missing-context reasons.
3. Run each exact context/chemistry/model once per audit. Native context scores
   cannot substitute for a changed final sequence. Unsupported chemistry,
   missing models, malformed arrays, empty output and backend failure remain
   explicit; unreturned bonds are unassessed, never assigned zero.
4. Inventory unfiltered patient-MHC output through Topiary and overlay every
   returned ligand occurrence and target interval. Separate internal from
   boundary cuts and distinguish sequence preservation from model evidence.
   Preserve partial HLA/model/length coverage without claiming completeness.
5. Write full native JSON and an escaped Jinja report with all positions,
   sources, settings and limitations. Keep proteasomal, extracellular/serum,
   cytosolic trimming, ER and endolysosomal evidence distinguishable. Whole-
   peptide half-life numbers do not stand in for site evidence.

Implementation may be delivered as a foundational site-contract/report PR and
a follow-up whole-product orchestration PR. #422 stays open until both native
and complete manufactured/translated contexts and ligand overlays are covered.
#423 then supplies independently documented Sid inputs and versioned real-model
regressions; this issue does not substitute synthetic scores for that evidence.

Validation: malformed/nonfinite/out-of-range arrays, both sequence ends,
repeated motifs, changed tails/junctions, deletion junction targets, non-mutation
targets, allele-free versus MHC output, missing models, unsupported chemistry,
native JSON round trips, and packaged report rendering. Run lint, full tests,
CLI/report smoke, CI, review and release gates for each PR.

Primary scope: Pepsickle's epitope model is proteasome-type agnostic, uses
sequence context and pads missing terminal context; its scores do not prove
destruction or presentation. Record the input-domain limits rather than hiding
them: https://github.com/pdxgx/pepsickle and
https://doi.org/10.1093/bioinformatics/btab628. mhctools 3.41.0's peptidase
recognition rules are not serum half-life or calibrated cleavage probabilities.

## Foundational site API (first PR; #422 remains open)

```python
from mhctools import CleavageInput
from mhctools.peptidases import get_cleavage_model
from vaxrank import (
    audit_pepsickle_inputs, audit_peptidase_inputs, write_cleavage_profiles,
)

# Explicit model input assumption: canonical linear L-peptide, free termini.
# Do not assert free termini if the manufactured chemistry is unknown.
inputs = [CleavageInput("GKRFHATISFDTDTGLKQALETKK", source_id="supplied-final-sequence")]
profiles = audit_pepsickle_inputs(inputs)
profiles += audit_peptidase_inputs(inputs, get_cleavage_model("cpn-basic"))
write_cleavage_profiles(profiles, json_path="sites.json", html_path="sites.html")
```

The Pepsickle dependency remains optional. Missing model assets are reported as
unassessed, not replaced by a predictor. Model weights are fingerprinted on the
active interpreter's import path; they are not redistributed. Terminal chemistry
outside a backend's domain remains unassessed. Arbitrary chemical modifications
are outside `CleavageInput` and must not be reduced to a bare sequence.
The remaining upstream Pepsickle per-bond/provenance API gap is tracked in
[mhctools #329](https://github.com/openvax/mhctools/issues/329). Once it ships,
Vaxrank should consume that result directly and remove its local model-metadata
adapter, without keeping a second long-term provenance implementation.

This first API does not infer patient typing, intended target ligands, historical
selection, or gene provenance from sequence. `profile.interval_evidence(start,
end)` provides exact internal and boundary observations for a caller-specified
half-open interval; the follow-up construct orchestration supplies validated
native/final/product contexts and patient-MHC overlays. The report lists every
internal bond, including those no model assessed. A terminal-only enzyme rule
does not become a full-sequence prediction or a serum-stability estimate.

## Complete-context API (3.17.0)

`audit_sequence_contexts` and `write_sequence_context_audits` now provide native,
final and validated complete-CDS contexts, full target/ligand overlays, and
explicit missing model/HLA/length/occurrence coverage. See
[CONSTRUCT_PROCESSING.md](CONSTRUCT_PROCESSING.md) for source-preserving inputs,
conditional chemistry, requested prediction provenance and JSON/report usage.
Historical Sid selection validation and the independent actual-construct audit
remain separate work under #414 and #423.
