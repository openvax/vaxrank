# Experimental host references for shared antigen inputs

Draft implementation for #625. The host reference and each antigen's source
species are separate: a human design can contain viral antigen sequences.

## First slice

Add optional `host_taxon_id` to shared input scope. Validate its type, agreement
between inputs, the configured/retained annotation, and parsed MHC allele species
before prediction. Use ORFdata's released taxonomy API and mhcgnomes allele
parsing, without maintaining a second species lookup table. Preserve the explicit
host declaration in native evidence and offline replay. Historical inputs with
no host declaration keep their existing identity and behavior.

`format` currently selects an input adapter (`cta_expression`, `antigens`,
`epitopes`, or an external producer format), rather than a file extension. Keep
the existing key compatible and explain it in the input guide. A future schema
can separate source kind from encoding if that becomes necessary.

## Remaining work before this PR is ready

- [ ] Consume released ORFdata composition and snapshot APIs from
  [#60](https://github.com/openvax/orfdata/pull/60) and
  [#68](https://github.com/openvax/orfdata/issues/68). Pin the tested release.
- [ ] Declare a composed host reference in the manifest with verified source
  hashes and explicit ORF evidence/coverage selection. Keep viral target layers
  separate from host self layers.
- [ ] Use the same resolved reference in candidate self matches and actual
  product-window audits. Preserve all mapped occurrences and explicit source
  exclusions; a gene exclusion must not silently exclude its additional ORFs.
- [ ] Retain the reference snapshot and its identity for deterministic native
  replay, corruption checks and reference-cache invalidation.
- [ ] Consume OncoRef's versioned species-specific CTA evidence from
  [#593](https://github.com/pirl-unc/oncoref/issues/593). Retain current human
  behavior until the new contract is available. Canine admission policy remains
  [Vaxrank #596](https://github.com/openvax/vaxrank/issues/596).
- [ ] Exercise released producer artifacts, uORF-only self matches, retained
  non-CTA sharing, missing nonhuman catalogue coverage and mixed product replay.

The first slice adds host validation only. Nonhuman CTA admission and expanded
self-reference checks remain pending above.
