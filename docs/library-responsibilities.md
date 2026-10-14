# Varcode, Isovar, and Vaxrank

| Library | Responsibility |
| --- | --- |
| **Varcode** | Predict coding consequences and transcript hypotheses from variants. |
| **Isovar** | Reconstruct RNA-supported sequences and compare them with those hypotheses. |
| **Vaxrank** | Score peptide candidates, select antigen windows and assemble vaccines. |

Sequences carry their reference, transcript IDs, frame/completeness and available
RNA evidence. Identical proteins retain all contributing sources; different
proteins remain separate candidates.

## Available today and remaining work

- The ordinary RNA path uses `IsovarResult.top_protein_sequence`, rather than all
  entries in `sorted_protein_sequences`.
- The [supplied-fusion adapter](how-it-works.md#fusion-antigens-from-supplied-assemblies)
  retains coding hypotheses and held-out reasons. It runs separately from the
  small-variant pipeline.
- Opt-in DNA fallback selects a usable changed protein from Varcode outcomes.
  These candidates have no RNA support.
- Automatic fusion discovery from RNA remains upstream work in
  [Isovar #305](https://github.com/openvax/isovar/issues/305).

## Other library guides

- [Varcode: hypotheses and coding consequences](https://openvax.github.io/varcode/library_roles/).
- [Isovar: RNA reconstruction and reconciliation](https://github.com/openvax/isovar/blob/master/docs/library-responsibilities.md).
