#!/usr/bin/env bash
set -euo pipefail

# See deploy.sh — Apple Silicon Homebrew + macOS SIP combine to strip dyld's
# search path when bash runs this script, so WeasyPrint can't find pango.
# Re-export here so direct `./test.sh` runs work the same as deploy-driven ones.
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/lib \
      && ":${DYLD_FALLBACK_LIBRARY_PATH:-}:" != *":/opt/homebrew/lib:"* ]]; then
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
fi

# -n 2 caps pytest-xdist at 2 workers. Each worker is its own Python
# process and loads its own pyensembl / mhctools / topiary state
# (~3 GB resident), so naive `-n auto` on an 8-core box can balloon
# to ~24 GB. Two workers is a sweet spot: meaningful parallelism
# without the memory blowup.
pytest -n 2 --cov=vaxrank/ --cov-report=term-missing tests
