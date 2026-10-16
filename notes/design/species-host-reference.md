# Experimental host references for shared antigen inputs

Experimental implementation for #625. The host reference and each antigen's source
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

## Scope and follow-ups

- [x] Consume released ORFdata composition and snapshot APIs from
  [#60](https://github.com/openvax/orfdata/pull/60) and
  [#68](https://github.com/openvax/orfdata/issues/68). Pin the tested release.
- [x] Declare a composed host reference in the manifest with verified source
  hashes and explicit ORF evidence/coverage selection. Keep viral target layers
  separate from host self layers.
- [x] Use the same resolved reference in candidate self matches and actual
  product-window audits. Preserve all mapped occurrences and explicit source
  exclusions; a gene exclusion must not silently exclude its additional ORFs.
- [x] Retain the reference snapshot and its identity for deterministic native
  replay, corruption checks and reference-cache invalidation.
- [ ] Consume OncoRef's versioned species-specific CTA evidence from
  [#593](https://github.com/pirl-unc/oncoref/issues/593). Retain current human
  behavior until the new contract is available. Canine admission policy remains
  [Vaxrank #596](https://github.com/openvax/vaxrank/issues/596).
- [ ] Exercise released producer artifacts, uORF-only self matches, retained
  non-CTA sharing, missing nonhuman catalogue coverage and mixed product replay.

Host validation and ORFdata snapshot consumption ship in this slice. Nonhuman CTA
admission and donor/transgene self layers remain pending.

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

Input simplification: use `host_species: mouse` and `mhc_species: human` in
ordinary examples. Accept common/scientific names and taxonomy IDs at the
boundary; require a list only for multiple MHC species. `host_species` is an
input alias for the existing canonical host taxon, not a second stored claim.
Keep canonical IDs and sorted MHC tuples internally so aliases, existing native
archives and input identities remain equivalent. Reject conflicting aliases.
Exercise manifest, API, CLI and sparse table inputs through export/replay.

Review fixes (#633–#635): normalize MHC common names with ORFdata's aliases;
resolve taxonomy IDs through an upstream public reverse lookup, independently
of Ensembl's supported annotations ([ORFdata #84](https://github.com/openvax/orfdata/pull/84),
released in 0.27.1). Read legacy native table declarations before
candidate adaptation, reconcile them with saved scope and CLI/manifest assertions,
and retain them through export/replay. Cover common-name/ID equivalence, cattle
without annotation, and sparse native CSV/TSV species declarations.

Pooled replay fix (#637): reconcile all saved input declarations before checking
each native source's host and alleles. Retain the existing per-source genotype,
reference and pooling checks. Carry the validated shared species scope into
every report and candidate so each exported source can reload independently.
Cover raw/native mixtures in either order, conflicts and CLI export/replay.

## Snapshot consumer

Accept an optional `host_reference` ORFdata snapshot path in the shared manifest
and `--host-reference` for standalone inputs. A `{path, identity}` mapping can
pin its content identity. Use the released `Proteome.load/save` contract; keep
annotated-only defaults unchanged. Validate the snapshot against the host,
assembly and annotation before models run. Reject conflicting snapshots.

Match host entries only, retaining every occurrence and the snapshot identity.
CTA gene exclusions apply to annotated proteins, never automatically to extra
ORFs. Keep complete ORF evidence and mapping history in the bundled snapshot.
Use this same reference for native candidates and final construct audits;
include its identity in replay decisions. Existing predictions must not silently
acquire new self-reference semantics on reload. Bundle the snapshot beside
native outputs with verified relative paths for offline relocation and replay.

Test noncanonical-only matches, multiple placements, retained non-CTA sharing,
viral layer separation, wrong references, corrupt/missing snapshots, and mixed
native replay. Synthetic cases test the adapter, not species catalogue coverage.
ORF definitions follow [Mudge et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9757701/).

## Release review

Ship host/MHC declarations and frozen host snapshots as experimental. Keep
species-specific CTA admission and donor/transgene layers as separate follow-ups.
Before release:

- Reconcile shared species before validating any new sequence source, including
  mixed saved/new inputs and either input order. Keep per-source genotype checks.
- Validate snapshot references using their own metadata, without requiring an
  Ensembl assembly alias, during prediction, native replay and junction queries.
- Apply the same final-product species checks to live genomes and snapshots;
  explicit MHC species overrides must work for both reference forms.
- Preserve pre-species assembly identities when the new fields are absent.
  Verify peptide/mRNA restoration from an actual 3.44.3 native export.
- Add regression cases, run lint/full tests and CLI replay, require green CI,
  then merge and deploy 3.45.0 from clean main.
