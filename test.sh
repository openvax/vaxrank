#!/bin/bash
# Coverage is disabled on Python 3.13 due to a pysam/Cython 3.1 bug:
# pysam's .pyx files enable Cython profiling, which uses sys.monitoring
# on 3.13+, causing "ValueError: Firing event 10 with no exception set"
# when coverage.py is active. Fixed on pysam master (commit 54f3d41) but
# not yet released. Remove this workaround once pysam >=0.24 is out.
PYTHON_VERSION=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
if [[ "$PYTHON_VERSION" == "3.13" ]]; then
  echo "Python 3.13 detected — running tests without coverage (pysam/Cython bug)"
  pytest tests
else
  pytest --cov=vaxrank/ --cov-report=term-missing tests
fi
