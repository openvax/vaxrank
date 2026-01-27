# Releasing Vaxrank

This document explains what to do once your [Pull Request](https://www.atlassian.com/git/tutorials/making-a-pull-request/) has been reviewed and all final changes applied. Now you're ready merge your branch into master and release it to the world:

0. Make sure that you have `pandoc` and `pypandoc` installed: this is needed for readme markdown on PyPI. (See [here](http://pandoc.org/installing.html) and [here](https://pypi.python.org/pypi/pypandoc), respectively, for instructions.)
1. Bump the [version](http://semver.org/) in `vaxrank/version.py`, as part of the PR you want to release.
2. Merge your branch into master.
3. Run `./deploy.sh` (builds with `build`, uploads with `twine`, and tags the release), or run the steps manually:
   - `python -m pip install --upgrade build twine`
   - `python -m build`
   - `python -m twine check dist/*`
   - `python -m twine upload dist/*`
   - `git tag "v${VERSION}"`
   - `git push --tags`
