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

## Review fixes and chimeric hosts (#630, #631)

- Normalize integer-valued numeric host columns and blanks at table ingestion;
  reject booleans, fractional IDs and infinities before candidate adaptation.
  Manifest/API taxon IDs remain strict positive integers.
- Add optional `mhc_species` (canonical concrete MHCgnomes species names). It
  declares the intended MHC repertoire independently of the animal/reference.
  With no override, check alleles against the known host/reference species.
  With an override, check every declared, observed and requested allele against
  that list; still check host taxon against annotation. Do not infer the
  transgene or knockouts from the animal's species.
- Support the same declarations in manifests, API scope and CLI. Validate
  input agreement, fresh predictions, native replay, direct VCF references and
  junction queries before model acquisition. Preserve both declarations and
  historical identities when they are absent.
- Regression coverage: plain TSV/CSV blanks, invalid numeric IDs, annotation-free
  dog/macaque conflicts, matching dog groups, HLA-only and mixed HLA/H2 mice,
  incompatible annotations/third-species alleles, export and offline replay.

HHD mice provide a concrete reason to separate these roles: their transgene
contains human peptide-binding domains and mouse H2-D domains. See the
[MGI transgene record](https://www.informatics.jax.org/allele/MGI%3A3711164).
This scope declaration does not compose a multi-species self proteome; the
reference snapshot work above still needs explicit layers and source evidence.

MHCgnomes casing follow-up: [#205](https://github.com/pirl-unc/mhcgnomes/issues/205).
The boundary uses canonical scientific-name casing for current releases.

Verification follow-up: exercise reference-only offline junction checks and
native table/comparator alleles too. Resolve any full-suite compatibility
failures before final validation; do not weaken the explicit species checks.

Reference-selection follow-up: honor the species of `--genome` when applying
`--ensembl-release`, reject incompatible assemblies, and test fresh HLA-mouse
antigen CLI loading with the actual resolver. Keep the human default only when
no genome is supplied. Document the mouse annotation flags explicitly.
