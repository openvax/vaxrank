# OncoRef 1.23 compatibility

Permit the verified OncoRef 1.23.0 release in Vaxrank's tested dependency
window so PirlyGenes 6.0.19 and Isovar's direct CTA catalogue can coexist.
Retain 1.19.2 as the lower endpoint for the currently published PirlyGenes
6.0.18. Historical annotation-bound CTA admission and exclusion semantics
remain unchanged; adopting the new species-scoped selection policy is #687.

Verify both endpoints with CTA identity, admission, species, expression and
exclusion regressions. Run the full suite on the current endpoint and check
the jointly installed OpenVax dependencies. Release 3.47.25 only after the
coordinated report sequence has published 3.47.24, then publish PirlyGenes
with its exact owner pin. Preserve data and source snapshot provenance.
