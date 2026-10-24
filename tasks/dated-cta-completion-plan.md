# Complete Vaxrank #665 after OncoRef #609

## Scope and compatibility

Finish the existing PR for #662; do not create a competing implementation. Rebase
its current remote head onto latest main, resolving version conflicts with the
next unused patch version. Preserve the reviewed resolved-assembly behavior and
exact int/string annotation provenance. Old numbered archives must retain their
saved fingerprints and field layout. Do not modify patient input, cached genome
files, retained admissions, source identities or the cache hashing algorithm.

Consume the published OncoRef fix at ==1.19.2 instead of the old exact 1.9.15
pin. Retain the exact pin so future upstream curation changes are deliberate. Verify current frozen human CTA APIs retain expected membership and alias
resolution. Keep PyEnsembl >=2.29 to supply resolved assembly metadata. Update
version/release documentation and remove stale comments/mocks for OncoRef #609.

## Evidence

- Existing PR identity baseline: 32 tests pass against its old OncoRef pin.
- Remove the monkeypatch substituting OncoRef's annotation validator in the dated
  admission/replay test, so that path uses the real public upstream implementation.
- Add a tiny synthetic GTF/cDNA/peptide fixture downloaded/indexed by a real dated
  EnsemblRelease with an explicit mirror server. Exercise public expression
  admission for dated GRCh37 and GRCh38 with independent accession/custom labels.
- Assert biological assembly, date string, original transcript/protein IDs,
  source exclusions and exact serialized replay. Preserve file bytes and stable
  file metadata after reopening/query/admission; no upstream identity lookups on
  persisted result replay. No additional scientific alias or protein rules.
- Keep the saved numbered fingerprint regression and fail early for unknown
  species/build/date metadata. Compare current reference membership hashes to the
  numbered baseline using the same published OncoRef tables.

## Validation and release

Use isolated temporary environments and APFS copies of the human75/human93/
mouse93/mouse102 and shared test-read/proteome caches. Never write to originals.
Run ./lint.sh, ./test.sh, focused identity and input-scope tests through the script,
CLI smoke, docs, pip check and package validation. Wait for every GitHub check on
the final head (Python3.11–3.14 and supported msgspec matrix). Update the existing
PR body to the final result, mark ready, merge and run ./deploy.sh on clean main.
Respect release_deploy.py's owned remote claim; verify published artifacts/tag.
Then proceed to PyEnsembl #456 and the related #497 validation regression.

## Confirmed upstream data delta to disclose

Comparing the actual published 1.9.15 data with the 1.19.2 code/data shows identical
624 admitted IDs and identical alias/canonical reference-table hashes. The full
candidate universe already gained SPAG9 (ENSG00000008294) in OncoRef's published
#613 recalibration: 2532 -> 2533. This is upstream data provenance, not a new
Vaxrank biological rule. Update the current-data expectation and document the
addition; preserve saved numbered record fingerprints and existing cached files.
Fresh candidate-universe hashes appropriately describe the selected OncoRef
version. Do not hide this delta with an ad-hoc excluded-gene list or rewrite old
snapshots to match new data.
