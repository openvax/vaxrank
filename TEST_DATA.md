# Reproducible osteosarc test data

## What travels with the tests

`tests/data/osteosarc/shared-v1/` contains **49 cases covering all 44 original
vaccine loci**, including both DYNC1H1 loci, plus a manifest. Its 98 original
BAM/index files total 4,126,730 bytes. Additional cases represent mitochondrial
platform/reference combinations and native-GRCh37 NR2F2, not additional vaccine
loci. These are Isovar's deliberately selected regression reads, not complete
platform coverage, unbiased expression measurements or validated vaccine targets.

The downloader generates this directory by exporting verified upstream subset
assets **without rewriting or further downsampling their reads**. Native alleles,
original GRCh38 identities where applicable, source BAM URLs, processing-product
identifiers, selection rules and source-region digests remain in the manifest.
Its upstream commit and manifest digest identify the fuller selection records
and upstream fixture generator. The reads originate in the public osteosarc
dataset ([CC0 source registry](https://registry.opendata.aws/sid-osteosarc/)).

This small bundle is separate from the ongoing multi-gigabyte, uncapped regional
analysis acquisition. It is not a replacement for that acquisition or an SV
discovery benchmark. Large analysis BAMs do not enter Git.

## Download, export and verify

Install Vaxrank/dependencies first. The packaged manifest lives at
`vaxrank/data/osteosarc-test-data-v1.json`; ordinary use needs neither an Isovar
checkout, a reference-genome download, samtools, nor the original large BAMs.

```sh
# Populate the shared cache only (explicit network operation).
python -m vaxrank.download_test_data

# Generate a new offline subset from those cached objects.
python -m vaxrank.download_test_data --offline --output /path/to/new/subset

# Verify the checked-in subset without network or cache mutations.
python -m vaxrank.download_test_data --verify-only \
  --output tests/data/osteosarc/shared-v1

# Regenerate into a NEW directory for byte-for-byte comparison/review.
python -m vaxrank.download_test_data --output /path/to/regenerated/subset
```

Existing identical exports are verified and reused; modified or unrelated output
directories are never repaired/overwritten. An export is published only after
every asset verifies. Interrupted acquisition leaves completed cache objects
reusable, but no finished-looking partial export. Individual interrupted file
transfers restart on retry; this does not claim byte-range download resumption.
Publication uses Linux/macOS atomic no-replace rename, so even an empty directory
created concurrently is preserved and validated, never replaced. If the OS or
filesystem lacks that operation, export fails safely without publishing output.
Corrupt cache hits fail visibly. `--repair-cache` explicitly redownloads invalid
objects; it cannot be combined with `--offline` and never repairs an export.

The source is pinned to Isovar commit
`0cad5b275c852263a1c77722aa463fe5568b2c76`, not a mutable branch or release lookup.
Every download and cache hit is checked for both size and SHA-256. The complete
manifest is shipped in both the package and exported dataset; tests require no
network or user cache to read this fixture. The fixture is repository-only;
sdists/wheels contain the downloader and its manifest, not the test BAMs.

## Shared OpenVax cache contract

Dataset identity: `osteosarc / vaccine-rna-v1` (independent of package version).
Default root: datacache's platform cache directory for **`openvax`**.
Override with `--cache-root /path/to/cache` or `OPENVAX_DATA_CACHE`.

Each object is stored at:

```text
<cache-root>/objects/sha256/<SHA-256><original suffixes>
```

For example, a BAM ends in `.bam` and its index in `.bam.bai`. Consumers can use
datacache 1.9.1+ directly, without importing Vaxrank:

```python
from pathlib import Path
from datacache import Cache

# asset is one entry from the pinned manifest's assets list.
cache = Cache("openvax", cache_root=Path(shared_root) / "objects" / "sha256")
filename = asset["sha256"] + "".join(Path(asset["filename"]).suffixes)
path = cache.fetch(asset["url"], filename=filename, timeout=60,
                   expected_sha256=asset["sha256"],
                   expected_size=asset["size_bytes"])
```

Content-identical assets reuse the same key across consumers and dataset
revisions. This is a cache convention, not a claim that Isovar, Varcode or
Topiary have already adopted it. Keep generated outputs outside the object
cache. A generic versioned bundle registry belongs in
[datacache #59](https://github.com/openvax/datacache/issues/59), not Vaxrank.
Adoption is tracked in [Isovar #308](https://github.com/openvax/isovar/issues/308),
[Varcode #464](https://github.com/openvax/varcode/issues/464), and
[Topiary #349](https://github.com/openvax/topiary/issues/349).

## Updating the pin

Maintainers can reproduce the packaged manifest from immutable Git objects in
an Isovar checkout containing the pinned commit:

```sh
python examples/osteosarc_test_data/pin_manifest.py --isovar-repo /path/to/isovar
```

The script verifies each upstream blob against the upstream manifest. Changing
the upstream revision or case membership requires an explicit reviewed dataset
revision, not silently overwriting `vaccine-rna-v1`. Generate a fresh export,
review provenance/size/biological scope, and run:

```sh
./lint.sh
./test.sh tests/test_download_test_data.py
./test.sh
```
