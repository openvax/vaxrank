#!/usr/bin/env bash
# Run the vaxrank test suite.
#
# Defaults to serial pytest. This is counterintuitive — the previous
# revision of this script went to lengths to pick a memory-aware
# pytest-xdist worker count — but measurement on the current suite
# shows xdist is consistently slower than serial regardless of worker
# count:
#
#   serial, no cov     49s   ← default
#   serial + --cov     78s
#   xdist -n 1         60s   (worker fork + IPC overhead)
#   xdist -n 2         95s   (slower than -n 1; disk contention)
#   xdist -n 4 + cov   often 20+ min (memory pressure → swap)
#
# Two structural reasons:
#
# 1. ``test_reference_peptide_logic`` runs ~20s on its own — first-time
#    mouse kmer-set pickle load. xdist's wall-time floor IS that test.
#    Parallelism over the other 831 fast tests can't beat the single
#    long pole.
# 2. Every xdist worker is a fresh Python process. Each independently
#    loads pyensembl + mhctools + topiary (~3 GB resident) AND
#    decompresses its own copy of the mouse kmer-set pickle.
#    Concurrent loads contend on disk I/O and (with --cov on) thrash
#    memory into swap.
#
# When/if those preconditions change (e.g. the slow test is split, or
# many new genuinely-parallel tests land), opt back into xdist with
# ``TEST_SH_XDIST=1``. The memory-aware worker math below still applies.
#
# Tunables (env vars):
#   TEST_SH_XDIST    set to 1 to use pytest-xdist (default: serial)
#   PER_WORKER_GB    per-worker memory budget in GB (default: 3.0)
#   TEST_SH_MIN      floor on workers when xdist is enabled (default: 1)
#   TEST_SH_MAX      hard ceiling on workers when xdist is enabled (default: unset)

set -eo pipefail

# See deploy.sh — Apple Silicon Homebrew + macOS SIP combine to strip dyld's
# search path when bash runs this script, so WeasyPrint can't find pango.
# Re-export here so direct `./test.sh` runs work the same as deploy-driven ones.
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/lib \
      && ":${DYLD_FALLBACK_LIBRARY_PATH:-}:" != *":/opt/homebrew/lib:"* ]]; then
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
fi

TEST_SH_XDIST="${TEST_SH_XDIST:-0}"
PER_WORKER_GB="${PER_WORKER_GB:-3.0}"
TEST_SH_MIN="${TEST_SH_MIN:-1}"
TEST_SH_MAX="${TEST_SH_MAX:-0}"

log() { printf '[test.sh] %s\n' "$*" >&2; }

case "$(uname -s)" in
    Darwin) OS=macos ;;
    Linux)  OS=linux ;;
    *)      OS=unknown ;;
esac

cpu_count() {
    local n=""
    if command -v getconf >/dev/null 2>&1; then
        n=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
    fi
    if [[ -z "$n" && "$OS" == "macos" ]]; then
        n=$(sysctl -n hw.logicalcpu 2>/dev/null || true)
    fi
    if [[ -z "$n" && -r /proc/cpuinfo ]]; then
        n=$(grep -c '^processor' /proc/cpuinfo 2>/dev/null || true)
    fi
    if [[ -z "$n" || "$n" -lt 1 ]]; then n=1; fi
    echo "$n"
}

cpu_cap() {
    local c="$1"
    if   (( c <= 1 )); then echo 1
    elif (( c <= 3 )); then echo $(( c - 1 ))
    else                    echo $(( c - 2 ))
    fi
}

mac_available_bytes() {
    local page_size
    page_size=$(sysctl -n hw.pagesize 2>/dev/null) || return 1
    vm_stat 2>/dev/null | awk -v ps="$page_size" '
        /Pages free/        { gsub(/\./, "", $3); free     = $3 }
        /Pages inactive/    { gsub(/\./, "", $3); inactive = $3 }
        /Pages speculative/ { gsub(/\./, "", $3); spec     = $3 }
        END { print (free + inactive + spec) * ps }
    '
}

linux_available_bytes() {
    [[ -r /proc/meminfo ]] || return 1
    awk '
        /^MemAvailable:/ { print $2 * 1024; found=1; exit }
        END              { if (!found) exit 1 }
    ' /proc/meminfo
}

available_bytes() {
    case "$OS" in
        macos) mac_available_bytes ;;
        linux) linux_available_bytes ;;
        *)     return 1 ;;
    esac
}

XDIST_FLAGS=()
if (( TEST_SH_XDIST == 1 )) && python -c "import xdist" 2>/dev/null; then
    CPUS=$(cpu_count)
    CPU_CAP=$(cpu_cap "$CPUS")

    avail=""
    if avail=$(available_bytes 2>/dev/null) && [[ -n "$avail" ]]; then
        MEM_CAP=$(awk -v b="$avail" -v g="$PER_WORKER_GB" 'BEGIN {
            n = int(b / 1024^3 / g)
            if (n < 1) n = 1
            print n
        }')
        AVAIL_GB=$(awk -v b="$avail" 'BEGIN { printf "%.1f", b / 1024^3 }')
        mem_note="ram_free=${AVAIL_GB}GB mem_cap=${MEM_CAP}"
    else
        MEM_CAP=$CPU_CAP
        mem_note="ram_free=? (probe unavailable) mem_cap=cpu_cap"
    fi

    if (( CPU_CAP < MEM_CAP )); then WORKERS=$CPU_CAP; else WORKERS=$MEM_CAP; fi
    if (( TEST_SH_MAX > 0 && WORKERS > TEST_SH_MAX )); then WORKERS=$TEST_SH_MAX; fi
    if (( WORKERS < TEST_SH_MIN )); then WORKERS=$TEST_SH_MIN; fi

    XDIST_FLAGS=(-n "$WORKERS")
    log "platform=${OS} cpus=${CPUS} cpu_cap=${CPU_CAP} ${mem_note} per_worker=${PER_WORKER_GB}GB"
    log "TEST_SH_XDIST=1 → xdist with -n ${WORKERS}"
    log "→ exec pytest -n ${WORKERS} --cov=vaxrank/ --cov-report=term-missing tests $*"
else
    if (( TEST_SH_XDIST == 1 )); then
        log "TEST_SH_XDIST=1 set but pytest-xdist not installed; running serial"
    else
        log "running serial (set TEST_SH_XDIST=1 to use pytest-xdist)"
    fi
    log "→ exec pytest --cov=vaxrank/ --cov-report=term-missing tests $*"
fi

exec pytest "${XDIST_FLAGS[@]}" --cov=vaxrank/ --cov-report=term-missing tests "$@"
