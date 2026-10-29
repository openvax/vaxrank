# Isovar consumer fixtures

These are unchanged JSON producer results, compressed with deterministic gzip.
`manifest.json` pins the SHA256 of each uncompressed JSON payload, its size and
its generating source. Upstream is Isovar 1.61.0, commit
`7bc88001704e5e8ff7c7bf2a0cd49cf1ea376b9e` (Apache-2.0):
https://github.com/openvax/isovar/tree/7bc88001704e5e8ff7c7bf2a0cd49cf1ea376b9e

The later `patient_background` test bundle is pinned by the repository-wide
`tests/data/manifest.json`. It contains public Isovar **1.76.0** exports, using `test_edit_attribution.py`
helpers from commit `511ae460a383485beff6f47c9c1f8e22dfdd8346`. It contains complete
protein-hypothesis and paired-background exports for synthetic focal/haplotype,
germline deletion, in-frame deletion, frameshift, tiled-read and two-transcript
cases. The Sid NTF3 case uses the original upstream read corpus and its local
reference reconstruction, case `28-NTF3-chr12-5494381-53f498a544883d51`.
The adjacent change is **supplied as germline for this regression**, exercising
retention of caller-declared edits; that classification is not independently
established by the RNA fixture. Sample/read-set labels `patient` / `RNA` are
test identities. No MHC measurements are included in this fixture.

- `ncorf`, `synonymous`, `tiled`, `unresolved`: the tiny noncoding transcript
  from `tests/test_variant_rna_orfs.py`, plus-strand fixture and its `assess`
  helper. They cover a changed residue, synonymous context, tiled RNA without
  an end-to-end witness, and no alternate reads. All sequence is synthetic.
- `sv`: `tests/test_sv_rna.py`'s synthetic `Scenario`/`reads`/`run` result.
  `sv_orfs` is the public `export_sv_rna_orfs` output from that same result.
- `sv_occurrences`: the synthetic two-path export from
  `tests/test_sv_rna_orf_export.py`'s union-witness regression. It records one
  fragment across two placements without the observation identity ledger.
- `isoforms`, `isoforms_partial`: `tests/test_antigen_isoforms.py`'s synthetic
  target models with two reads of the first complete transcript or its first
  120 bases. The unobserved competing reference transcript remains present.
- `viral`, `viral_partial`: two synthetic reads from the unchanged public HPV16
  reference fixture used by upstream `tests/test_viral_coding.py`. These exercise
  full/partial reference-anchored coding and source provenance. They are **not**
  real SiHa RNA, tumor-specificity evidence or expression measurements.
- `colo829_dynll2`: one complete producer result from the released
  `reports/hcc1395_colo829/ncorf_rna_audit/results.json.gz` audit. The selected
  DYNLL2 result has the smallest candidate catalog in that audit (300 records).
  No protein candidate, competing model, RNA alternative, witness set, support
  value or producer setting was removed. The upstream report documents the
  original COLO829 public RNA and its reference/call provenance.

Temporary source paths in synthetic results are retained as the producer wrote
them. Tests replay the frozen outputs without opening those paths, querying an
annotation, downloading source reads, or claiming the synthetic fixtures are
biological observations. Separate tests execute the live public APIs using tiny
local BAMs/reference models, and verify CLI design/report/native replay.
