# Releasing Vaxrank

This document explains what do once your [Pull Request](https://www.atlassian.com/git/tutorials/making-a-pull-request/) has been reviewed and all final changes applied. Now you're ready merge your branch into master and release it to the world:

1. Bump the [version](http://semver.org/) on __init__.py, as part of the PR you want to release.
2. Merge your branch to master.
2. Run `python setup.py sdist upload`, which pushes the newest release to PyPI.
