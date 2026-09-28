# Releasing Vaxrank

Use Python 3.10 or newer. Set `PYTHON` to the environment used for lint and
tests so the release uses the same dependencies. The package reads README.md
directly as Markdown; Pandoc and pypandoc are not required.

1. Bump the version in `vaxrank/version.py` as part of the PR, including for documentation-only changes.
2. After lint, tests and GitHub CI pass, merge the PR, check out main, and pull the merge. Confirm the working tree is clean.
3. Run `./deploy.sh` without a version argument. It uses one Python environment for lint, tests, build,
   and upload; verifies every published artifact by SHA-256; and pushes the
   release tag only after PyPI contains the complete matching release.

If an upload is interrupted after the local release tag is created, leave the
tag and `dist/` intact and rerun `./deploy.sh` from the same clean checkout. The
script reuses those exact artifacts, verifies any file already on PyPI, and
uploads only a missing file. It refuses a retry when an original artifact is
missing or its bytes differ from PyPI.
