# Development

```bash
git clone git@github.com:openvax/vaxrank.git
cd vaxrank
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -e ".[hitlist]" pytest pytest-cov ruff mkdocs-material
# Reference releases exercised by the test suite
pyensembl install --release 75 93 --species human
pyensembl install --release 93 102 --species mouse
```

Install the [PDF dependencies](installation.md#pdf-reports) and `samtools`
for the original-read integration tests. Run linting and tests:

```bash
./lint.sh && ./test.sh
```

See [CONTRIBUTING.md](https://github.com/openvax/vaxrank/blob/main/CONTRIBUTING.md)
for pull-request expectations and
[RELEASING.md](https://github.com/openvax/vaxrank/blob/main/RELEASING.md) for
the release procedure.

## Test suite and test data

The test suite is intentionally repository-only and is not included in source
distributions: its integration tests depend on helpers and genomic fixtures
that are not package payloads. Run tests from a Git checkout rather than an
unpacked PyPI archive.

The [Sid test data](https://github.com/openvax/vaxrank/blob/main/TEST_DATA.md)
contains only the reads the regression suite explicitly selects. They come from
openvax-v1, the OpenVax libraries' shared Sid test data published by osteosarc;
the reviewed recipe ships in the package and verifies every selected record.
The first run downloads required reads and reference data. Cached runs can
reuse them; the full suite also contains explicit download and integration checks.

The first run may take a while to build the reference proteome kmer index.
Later runs use the cached index.

## Reference membership caches

Caches are keyed by installed annotation/sequence content, or by actual
protein content for custom reference providers, as well as species and release. This supports explicit offline PyEnsembl `Genome` references and keeps
a subset from sharing a full-reference cache. Upgrading from
species/release-only caches triggers a one-time rebuild; old cache files are
not trusted and not automatically deleted.

## Scripts

- `develop.sh`: installs the package in editable mode.
- `lint.sh`: runs Ruff on the package, tests and release helpers.
- `test.sh`: runs pytest; CI additionally collects coverage.
- `deploy.sh`: runs lint, tests, build and upload through one Python
  environment, verifies PyPI artifacts by SHA-256, and safely resumes an
  interrupted upload from its original artifacts.
- `check_release_version.py --base origin/main`: compares the committed PR
  version with the fetched target branch and release tags, without importing
  either version file.

### Deploying a release

Deploy is restricted to a synchronized `main`/`master`. Bump the version in a
PR, then run `./deploy.sh` after merging. An optional version argument (with or
without a `v` prefix) must match `vaxrank/version.py`; it never edits or
commits a version bump.

After the quality gates, `release_deploy.py` rechecks live main and claims
`refs/tags/release-locks/vVERSION` on `origin` before building. Independent
clones wait for that owner, then verify the completed release instead of
rebuilding it. Annotated release tags record the source commit and the exact
wheel/sdist SHA-256 hashes, which must match PyPI. `--dry-run` validates and
builds an archive of committed sources in temporary directories without
uploading, tagging, pushing, or replacing a running deploy's build files.

An interrupted upload is resumed with its original `dist/` files and local
annotated release tag; missing or changed artifacts are rejected. Normal
failures release the remote claim, but a killed process or lost connection may
leave one behind. Before removing a stale claim, confirm its owner has stopped
and recover its original artifacts and tag. The error message prints a deletion
command guarded by that claim's exact object ID, so it cannot delete a newer
owner's claim. Published versions with older lightweight tags have no recorded
hashes and need manual artifact/tag verification; they are never rebuilt
automatically.

### Release version check in CI

The Release version PR check fetches the current target branch and tags on
every run, including when a PR is retargeted. A PR must increase the version
and choose a version with no existing release tag. Fetch the target branch and
tags before running `check_release_version.py` locally
(`git fetch --tags origin`).

Closed PR events are skipped. If validation fails while an open PR's job is
running, the workflow rechecks the live PR state before reporting a failure;
open PRs and failed state lookups still fail validation. A green check only
covers the base and tags available during that run. Require the check with
GitHub's "Require branches to be up to date before merging" option to prevent
an old green result surviving another PR's release; otherwise rerun it after
the target branch advances. The workflow itself does not change branch
protection.

### CI installs

CI uses `scripts/pip_install.sh` for installs that resolve sibling releases. It
makes at most three attempts, waiting 20 then 40 seconds, and bypasses pip's
cache on retries. All requirements and constraints stay unchanged, and
persistent failures retain pip's exit status. `PYTHON`,
`PIP_INSTALL_ATTEMPTS` and `PIP_INSTALL_RETRY_DELAY` override the interpreter
and retry settings.
