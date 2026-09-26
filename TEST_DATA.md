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

## Use

```python
from vaxrank.sid_test_data import sid_test_data, sid_reads, sid_variants

root = sid_test_data()  # builds the test files once per process
cases = root / "osteosarc/shared-v1/manifest.json"
variants = sid_variants(["MAP2-chr2-209694768"])
# sid_reads("osteosarc/shared-v1/<case>.bam") returns osteosarc.ReadSubset.
```

`sid_test_data()` builds the files into a temporary process directory. It
exports the 58 members from openvax-v1, selects each cohort's records by
digest (`osteosarc.cohort_bundle.select_records`), and writes them with the
recipe's reviewed header and order. Missing, changed or extra records fail the
build. The files are the same bytes as the zip Vaxrank shipped through 3.25.

The first build downloads and verifies openvax-v1 (28 MB) into the osteosarc
cache (`OSTEOSARC_CACHE`, else the shared OpenVax cache); later builds work
offline. `provenance.json` in the built directory records the recipe hashes,
the openvax-v1 manifest checksum, each cohort's source identity and the
osteosarc version.

The existing CLI exports the 49 retrieval cases:

```bash
python -m vaxrank.download_test_data --output /tmp/sid-retrieval-tests --offline
python -m vaxrank.download_test_data --output /tmp/sid-retrieval-tests --verify-only
```

Existing output is verified, never replaced. Explicit custom `--manifest`
downloads retain the generic shared OpenVax cache API for existing callers.
They do not participate in the Sid tests.

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
