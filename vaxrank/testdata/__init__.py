"""Vaxrank's regression test data: how it is described, built and downloaded.

* ``sid`` builds the small Sid regression files from a reviewed recipe
  (``sid-recipe/``, shipped as package data) and the shared ``openvax-v1`` reads
  published by osteosarc. The recipe stores digests, headers and expectations,
  never reads.
* ``download`` is the ``vaxrank-test-data`` command, which exports verified test
  data into the shared OpenVax cache.

This lives in the package, not in ``tests/``, because the command ships with it
and ``tests/`` is excluded from the sdist.
"""
