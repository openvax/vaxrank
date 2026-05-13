#!/usr/bin/env bash
set -euo pipefail

# See deploy.sh — Apple Silicon Homebrew + macOS SIP combine to strip dyld's
# search path when bash runs this script, so WeasyPrint can't find pango.
# Re-export here so direct `./test.sh` runs work the same as deploy-driven ones.
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/lib \
      && ":${DYLD_FALLBACK_LIBRARY_PATH:-}:" != *":/opt/homebrew/lib:"* ]]; then
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
fi

# Use pytest-xdist with a memory-aware worker count. Each worker loads
# its own pyensembl / mhctools / topiary state (~3 GB resident), so
# naive `-n auto` on an 8-core box can balloon to ~24 GB — and gets
# even worse when sibling pytest suites (hitlist, trufflepig) are
# running concurrently. We pick the smaller of (cores, available_RAM
# / per_worker_gb), and drop to 1 worker when the machine is already
# under strain. xdist is not a required dependency, so fall back to
# serial pytest when it isn't installed (e.g., CI).

readonly PER_WORKER_GB=3.0

# macOS available-RAM heuristic: free + inactive + speculative pages
# are reclaimable on demand, so they count as headroom.
macos_available_bytes() {
  local page_size pages
  page_size=$(sysctl -n hw.pagesize)
  pages=$(vm_stat | awk '
    /Pages free/        { gsub(/\./, "", $3); free     = $3 }
    /Pages inactive/    { gsub(/\./, "", $3); inactive = $3 }
    /Pages speculative/ { gsub(/\./, "", $3); spec     = $3 }
    END                 { print free + inactive + spec }
  ')
  echo $(( pages * page_size ))
}

# Pick a pytest -n value that respects both CPU and available RAM.
pytest_workers() {
  local cpus avail_bytes
  cpus=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.logicalcpu)
  if [[ "$(uname)" == "Darwin" ]]; then
    avail_bytes=$(macos_available_bytes)
  else
    avail_bytes=$(awk '/MemAvailable/ { print $2 * 1024 }' /proc/meminfo)
  fi
  awk -v cpus="$cpus" -v bytes="$avail_bytes" -v budget="$PER_WORKER_GB" '
    BEGIN {
      by_memory = int(bytes / 1024^3 / budget)
      n = (cpus < by_memory ? cpus : by_memory)
      print (n < 1 ? 1 : n)
    }
  '
}

XDIST_FLAG=""
if python -c "import xdist" 2>/dev/null; then
  workers=$(pytest_workers)
  XDIST_FLAG="-n $workers"
  echo "Running pytest with -n ${workers} (per-worker budget ≈ ${PER_WORKER_GB} GB)"
fi
pytest ${XDIST_FLAG} --cov=vaxrank/ --cov-report=term-missing tests
