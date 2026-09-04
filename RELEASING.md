# Releasing Vaxrank

This document explains what to do once your [Pull Request](https://www.atlassian.com/git/tutorials/making-a-pull-request/) has been reviewed and all final changes applied. Now you're ready to merge your branch into main and release it to the world:

0. Make sure that you have `pandoc` and `pypandoc` installed: this is needed for readme markdown on PyPI. (See [here](http://pandoc.org/installing.html) and [here](https://pypi.python.org/pypi/pypandoc), respectively, for instructions.)
1. Bump the [version](http://semver.org/) in `vaxrank/version.py`, as part of the PR you want to release.
2. Merge your branch into main, check out main, and pull the merge.
3. Run `./deploy.sh`. It uses one Python environment for lint, tests, build,
   and upload; verifies every published artifact by SHA-256; and pushes the
   release tag only after PyPI contains the complete matching release.

If an upload is interrupted after the local release tag is created, leave the
tag and `dist/` intact and rerun `./deploy.sh` from the same clean checkout. The
script reuses those exact artifacts, verifies any file already on PyPI, and
uploads only a missing file. It refuses a retry when an original artifact is
missing or its bytes differ from PyPI.
