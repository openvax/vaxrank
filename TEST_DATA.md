# Bundled Sid test data

Vaxrank packages **1,148 selected alignment records in 58 cohorts** in
`vaxrank/data/sid-test-data.zip`. All Sid read fixtures are acquired through
**osteosarc 0.1.2** from the public
[CC0 Sid dataset](https://registry.opendata.aws/sid-osteosarc/).
The archive contains no whole-source BAM, BAM index, dataset snapshot, or
unselected regional reads. Its BAM indexes describe only the tiny selected BAMs.

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

Small Ensembl reference subsets, documented biological expectations, cached
NetMHCpan outputs, and report expectations are included as separate supporting
inputs. Osteosarc supplies the Sid acquisition/variant identities; it does not
supply the independent Ensembl annotations or MHC prediction expectations.

## Offline use

```python
from vaxrank.sid_test_data import sid_test_data, sid_reads, sid_variants

root = sid_test_data()  # verifies and opens installed package resources
cases = root / "osteosarc/shared-v1/manifest.json"
variants = sid_variants(["MAP2-chr2-209694768"])
# sid_reads("osteosarc/shared-v1/<case>.bam") returns osteosarc.ReadSubset.
```

The tests use this same installed-package loader. It writes only a temporary
process directory, needs no persistent cache, and performs no downloads.
The existing CLI exports the 49 retrieval cases offline:

```bash
python -m vaxrank.download_test_data --output /tmp/sid-retrieval-tests --offline
python -m vaxrank.download_test_data --output /tmp/sid-retrieval-tests --verify-only
```

Existing output is verified, never replaced. Explicit custom `--manifest`
downloads retain the generic shared OpenVax cache API for existing callers.
They do not participate in the bundled Sid tests.

## Regenerate the package subset

Install the declared requirements and SAMtools 1.21 or later, then run from the
repository root (the recipe and generator are also included in the sdist):

```bash
python examples/osteosarc_test_data/build.py --cache /tmp/sid-acquisition
```

This works with an empty acquisition cache. The checked-in small catalogue
contains the selected native `osteosarc.Asset` and `Variant` identities from
snapshot `e4224eeedc18b9a0e66d9e57afcb5f9d613ed76a775e58cd56d320937cd6cf6f`,
plus pinned index receipts. It does not require the full historical metadata
snapshot. `osteosarc.extract_reads` retrieves indexed regions and verifies
source/index identity; the generator then retains **only** the records listed
in `recipe/selection.json.gz`. A minimal set of one-base retrieval anchors covers the selected alignment
spans, including explicitly required off-locus mates or supplementary records.
The `retrieval_regions` helper derives this point cover when preparing a recipe,
avoiding a separate remote seek for every selected read. The temporary regional results remain in the acquisition cache and
never enter the package.

The allowlist stores alignment digests, order and multiplicity, not read bases.
Native BAM digests include float-tag bits lost by SAM text formatting. Historical
SAM/fusion inputs use their original text representation. Missing, changed or
extra identical selected records fail generation. A freshly acquired source
cannot silently redefine a regression baseline.

Repeat against the same verified cache without network access:

```bash
python examples/osteosarc_test_data/build.py --cache /tmp/sid-acquisition --offline
```

The read selection is deterministic. Acquisition receipts preserve real source
identity, software versions and retrieval evidence; archive bytes can differ
between cache locations/tool versions even when all selected records agree.
`bundle.json` verifies every packaged file, and `provenance.json` records the
recipe hashes and full osteosarc extraction lineage.

To update the source catalogue deliberately, first review/update the selector
recipe and create the named osteosarc metadata snapshot, then run:

```bash
python examples/osteosarc_test_data/pin_catalog.py --cache /path/to/snapshot-cache
```

The published-allele baseline is explicit (`corrections=False`). In particular,
MAP2's old deletion remains the historical selection baseline; the separate
corrected-complex-allele regression runs against the same selected original
reads. This migration does not change either biological expectation.

After regeneration, run `./lint.sh` and `./test.sh`. The tests check every
selected record, native mitochondrial/GRCh37 retrieval, reconstruction and
ranking, bundle integrity, and the actual wheel/sdist payload.
