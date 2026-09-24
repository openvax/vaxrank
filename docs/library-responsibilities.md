# Varcode, Isovar, and Vaxrank

Vaxrank evaluates protein/peptide candidates. It should preserve the evidence
behind a sequence, not reinterpret an RNA junction or select a transcript merely
because its effect class ranks first.

| Library | Responsibility |
|---|---|
| **Varcode** | Generate structural/transcript hypotheses and predict their coding consequences. |
| **Isovar** | Reconstruct RNA-supported sequences, compare them with those hypotheses, and preserve unresolved alternatives. |
| **Vaxrank** | Evaluate the resulting protein/peptide candidates, retaining their evidence. |

This is the shared responsibility split; it does not imply that all candidate
paths are already connected to the default pipeline.

## What Vaxrank receives

A sequence candidate needs its nucleotide/transcript source, compatible
transcript IDs, reference/annotation identity, completeness, frame assumptions,
and sample/library/read evidence where available. Keep RNA-supported and
DNA-only predicted sequences distinguishable.

Vaxrank uses these candidates for peptide generation, prediction, filtering,
ranking, and vaccine design. Binding or manufacturability scores do not resolve
which RNA isoform is expressed. RNA support does not itself prove translation,
antigen presentation, tumor specificity, or immunogenicity.

If candidates encode the same protein, grouping them must retain all source
structures and provenance. Different proteins remain distinct candidates.
Missing RNA is not proof of no expression; insufficient evidence should remain
visible rather than becoming an unsupported confidence claim.

## Available today and remaining work

- **Ordinary RNA path:** `from_isovar_result` uses
  `IsovarResult.top_protein_sequence`. Isovar can retain more alternatives in
  `sorted_protein_sequences`; the ordinary Vaxrank path does not evaluate all
  of them automatically.
- **Supplied-fusion path:** `fusion_antigens_from_isovar` retains coding
  hypotheses and reports unresolved/ambiguous status. It is separate from the
  small-variant pipeline; see
  [fusion inputs](https://github.com/openvax/vaxrank/blob/main/README.md#upstream-inputs).
- **Opt-in DNA fallback:** consequences come from Varcode, not RNA
  reconstruction. Fusion candidates are checked for a usable changed protein
  before effect selection ([#482](https://github.com/openvax/vaxrank/issues/482),
  fixed in 3.20.2). A set-level protein-change flag does not guarantee its first
  protein is changed.
- **RNA discovery/reconciliation:** Isovar's supplied-fusion translator is not
  an automatic soft-clip assembly workflow. That upstream connection is
  [Isovar #305](https://github.com/openvax/isovar/issues/305).

The intended contract is to retain supported alternatives and their uncertainty
through evaluation. A consumer that requires one result must make the selection
policy explicit; effect severity, RNA support, and peptide ranking answer
different questions.

## Other library guides

- [Varcode: hypotheses and coding consequences](https://openvax.github.io/varcode/library_roles/).
- [Isovar: RNA reconstruction and reconciliation](https://github.com/openvax/isovar/blob/master/docs/library-responsibilities.md).

