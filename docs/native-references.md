# Moving native mutation datasets

Native epitope exports retain their evidence, scores, mutation fragments and
supporting transcript objects. Restoring transcript identity no longer opens an
annotation database; Vaxrank requires PyEnsembl 2.22.2 for that contract.

When exporting a custom Genome, Vaxrank copies available annotation resources
into a sibling directory, for example:

```
result.tsv
result.tsv.references/
```

Move both together. Each bundled resource has a relative path, byte size and
SHA-256 digest in the versioned `vaxrank.native_reference.v1` manifest. Identical
files share a copy. This can include complete GTF/transcript/protein/reference
DNA files, so the bundle size depends on the supplied custom reference.
Export does not download remote sources. Standard Ensembl references retain
their exact release/species identity; explicitly attached local reference DNA
is bundled too.

Loading verifies available bundle files before rebinding paths. A checksum
conflict fails explicitly, even if an original annotation still exists. A
missing bundle file remains missing; Vaxrank does not silently substitute the
old source path or another annotation release. Gene/transcript/reference
identity and saved protein IDs survive relocation and re-save.

Scoring and construct replay use the retained evidence without predictor calls
or annotation downloads. Even a table moved without its resource directory can
replay those measurements; annotation-dependent work needs the resources.
Availability is inspectable through `dataset.native_references.availability`.

To prepare bundled local annotation indexes for full reports:

```bash
vaxrank --input-epitopes result.tsv --index-native-references --output-dir reports
```

The equivalent Python call is `dataset.index_native_references()`. Indexing
verifies the resources again and acquires no data. Custom-reference indexes are
local to the relocated bundle; standard Ensembl releases use PyEnsembl's normal
release cache. Missing annotation resources produce an explicit error.

Earlier native files retain their original reference identities and now restore
metadata lazily. Where the saved antigen already records protein IDs for the
same supporting transcripts, replay reuses them. Legacy files have no verified
resource manifest; Vaxrank cannot infer a replacement annotation or recover an
unavailable source. Re-export with the recorded reference available to bundle
it.
