# Sid test data

Vaxrank's Sid regression tests use **1,148 selected alignment records in 58
cohorts** from the public
[CC0 Sid dataset](https://registry.opendata.aws/sid-osteosarc/). Vaxrank stores
no reads. The records come from **openvax-v1**, the OpenVax libraries' shared
Sid test data, published by osteosarc 0.11
([iskandr/osteosarc#56](https://github.com/iskandr/osteosarc/issues/56)). Each
cohort is the openvax-v1 member `vaxrank/<path>`, e.g.
`vaxrank/osteosarc/shared-v1/00-ABCF2-chr7-151218156-f30f618fb76a0e49.bam`.

The reviewed recipe ships in the package, in `vaxrank/data/sid-recipe`:

* `selection.json.gz` lists every cohort's records by exact digest, in order
  and with multiplicity;
* `headers/` holds each cohort's reviewed SAM header (no records);
* `fusion/` holds the two fusion inputs' RNA windows and reference hypotheses;
* `catalog.json` pins the native `osteosarc.File` and `Variant` identities from
  snapshot `e4224eeedc18b9a0e66d9e57afcb5f9d613ed76a775e58cd56d320937cd6cf6f`;
* `support/` holds the non-read inputs: small Ensembl reference subsets,
  documented biological expectations, cached NetMHCpan outputs and report
  expectations.

The selection is explicit:

* 49 native-coordinate retrieval cases cover the 44 original vaccine loci.
  Each retains one template, including its already-selected mate/alternative
  alignments. NTF3 retains a template carrying the compound AG>GT evidence.
* Five reconstruction/ranking cohorts and two context cohorts retain their
  reviewed records. Their counts, competing alignments, low-quality reads and
  ambiguity are part of the regressions, so these are not resampled.
* Two fusion cases retain seven original SAM records and their reviewed RNA
  windows/reference hypotheses.

These are selected test cohorts, **not coverage, VAF or expression estimates**.
Mates and supplementary records are explicitly listed; this does not claim all
mates or all alignments of every template in the original BAM were recovered.

Osteosarc supplies the Sid reads and variant identities; it does not supply
the independent Ensembl annotations or MHC prediction expectations.

## Fixtures Osteosarc does not supply

Osteosarc publishes the human Sid dataset, so two groups of fixtures stay in
this repository permanently rather than pending migration:

- `tests/data/b16.f10/` — mouse B16-F10 reads and VCFs. They drive the CLI and
  report smoke (`run-vaxrank-b16-test-data.sh`) and the mutant-protein-sequence
  tests. Osteosarc has no mouse content, so there is nothing to migrate to.
- `tests/data/epitope_fixtures/` — LENS and pVACseq epitope tables for the
  external-input, rescoring and report paths. Osteosarc publishes reads, not
  epitope tables.

Both groups are pinned by sha256 in `tests/data/manifest.json` and verified by
`tests/test_data_manifest.py`, which also fails on a fixture nobody pinned. That
check is the substitute for osteosarc's provenance, not an equivalent of it:
these files carry digests, not a recorded upstream source. Regenerate the
digests with `python tests/test_data_manifest.py` when changing a fixture is the
point of the commit.

## Use

```python
from vaxrank.sid_test_data import sid_test_data, sid_reads, sid_variants

root = sid_test_data()  # builds the test files once per process
cases = root / "osteosarc/shared-v1/manifest.json"
variants = sid_variants(["MAP2-chr2-209694768"])
# sid_reads("osteosarc/shared-v1/<case>.bam") returns osteosarc.ReadSubset.
```

`sid_test_data()` builds the files into a temporary process directory. It
reads the 58 members from openvax-v1 with `osteosarc.bundle_file`, which
exports each into the osteosarc cache once, selects each cohort's records by
digest (`osteosarc.cohort_bundle.select_records`), and writes them with the
recipe's reviewed header and order. Missing, changed or extra records fail the
build. The files are the same bytes as the zip Vaxrank shipped through 3.25.

The first build downloads and verifies openvax-v1 (28 MB) into the osteosarc
cache (`OSTEOSARC_CACHE`, else the shared OpenVax cache) and exports the
members there; later builds reuse them offline. `provenance.json` in the built directory records the recipe hashes,
the openvax-v1 manifest checksum, each cohort's source identity and the
osteosarc version.

`vaxrank-test-data` exports the 49 retrieval cases (`python -m
vaxrank.download_test_data` is equivalent):

```bash
vaxrank-test-data --output /tmp/sid-retrieval-tests --offline
vaxrank-test-data --output /tmp/sid-retrieval-tests --verify-only
```

Existing output is verified, never replaced. Explicit custom `--manifest`
downloads retain the generic shared OpenVax cache API for existing callers.
They do not participate in the Sid tests.

Every run prints one JSON object and exits non-zero only on failure, so a
caller reads a status rather than parsing prose:

```json
{
  "action": "export",
  "assets": 98,
  "cached_paths": null,
  "data_version": "minimal-vaccine-rna-v2",
  "dataset": "osteosarc",
  "error": null,
  "files": {"00-ABCF2-chr7-151218156-f30f618fb76a0e49.bam": {"status": "available", "verified": true}},
  "output": "/tmp/sid-retrieval-tests",
  "status": "available"
}
```

`status` is `available` only when every asset's digest and size were checked
and matched. Otherwise it names the problem — `corrupt`, `missing`,
`inaccessible`, `unexpected_contents`, `manifest_mismatch`, `symlink`,
`not_a_directory`, or `unavailable` when the manifest itself could not be
read — and `error` carries the prose. Without `--output`, `status` is `cached`
and `cached_paths` gives each asset's location in the shared cache.
`--progress`, `--timeout` and `--max-retries` apply to custom manifests, whose
assets are the only ones this command downloads.

The cache root is `--cache-root`, else `OSTEOSARC_CACHE`, else
`OPENVAX_DATA_CACHE`, else the platform cache directory — the order osteosarc
itself resolves, so the bundled reads and any custom assets in one invocation
always come from the same root.

From Python:

```python
from vaxrank.download_test_data import (
    download_test_data, inspect_dataset, verify_dataset)

report = inspect_dataset("/tmp/sid-retrieval-tests")  # never raises
if report["status"] != "available":
    print(report["detail"], report["files"])

verify_dataset("/tmp/sid-retrieval-tests")  # raises on the first problem
download_test_data("/tmp/sid-retrieval-tests", offline=True)
```

`inspect_dataset` answers "what is wrong with this dataset", `verify_dataset`
answers "is it usable". Both check structure that datacache cannot express:
`inspect_files` follows symlinks and inspects only the inventory it is handed,
so it would accept a symlinked asset and never notice an extra file.

## Changing the recipe

The allowlist stores alignment digests, order and multiplicity, not read bases.
Native BAM digests include float-tag bits lost by SAM text formatting.
Historical SAM/fusion inputs use their original text representation. A cohort
can only use records that openvax-v1 holds; adding records means adding them
to the shared test data in osteosarc first.

To update the source catalogue deliberately, first review/update the recipe and
create the named osteosarc metadata snapshot, then run:

```bash
python examples/osteosarc_test_data/pin_catalog.py --cache /path/to/snapshot-cache
```

The published-allele baseline is explicit (`corrections=False`). In particular,
MAP2's old deletion remains the historical selection baseline; the separate
corrected-complex-allele regression runs against the same selected original
reads.

After a change, run `./lint.sh` and `./test.sh`. The tests check every
selected record, native mitochondrial/GRCh37 retrieval, reconstruction and
ranking, and the actual wheel/sdist payload.
