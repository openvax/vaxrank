#!/usr/bin/env bash
# Bounded wkhtmltopdf installer for Ubuntu CI runners (openvax/vaxrank#806).
# Preserves the pinned 0.12.6-1 release and failure checks while ensuring a
# stalled apt or upstream GitHub download fails quickly with a clear message.
set -euo pipefail

WKHTMLTOPDF_DEB_URL="${WKHTMLTOPDF_DEB_URL:-https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6-1/wkhtmltox_0.12.6-1.bionic_amd64.deb}"
WKHTMLTOPDF_DEB_PATH="${WKHTMLTOPDF_DEB_PATH:-wkhtmltox_0.12.6-1.bionic_amd64.deb}"
WKHTMLTOPDF_APT_TIMEOUT="${WKHTMLTOPDF_APT_TIMEOUT:-120}"
WKHTMLTOPDF_DOWNLOAD_TIMEOUT="${WKHTMLTOPDF_DOWNLOAD_TIMEOUT:-120}"
WKHTMLTOPDF_DPKG_TIMEOUT="${WKHTMLTOPDF_DPKG_TIMEOUT:-60}"
WKHTMLTOPDF_CONNECT_TIMEOUT="${WKHTMLTOPDF_CONNECT_TIMEOUT:-20}"
WKHTMLTOPDF_READ_TIMEOUT="${WKHTMLTOPDF_READ_TIMEOUT:-30}"
WKHTMLTOPDF_RETRIES="${WKHTMLTOPDF_RETRIES:-3}"
WKHTMLTOPDF_WAITRETRY="${WKHTMLTOPDF_WAITRETRY:-5}"
PYTHON_BIN="${PYTHON:-python3}"

for var_name in \
    WKHTMLTOPDF_APT_TIMEOUT \
    WKHTMLTOPDF_DOWNLOAD_TIMEOUT \
    WKHTMLTOPDF_DPKG_TIMEOUT \
    WKHTMLTOPDF_CONNECT_TIMEOUT \
    WKHTMLTOPDF_READ_TIMEOUT \
    WKHTMLTOPDF_RETRIES \
    WKHTMLTOPDF_WAITRETRY
do
    val="${!var_name}"
    if [[ ! "$val" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: ${var_name} must be a positive integer (got: ${val})" >&2
        exit 2
    fi
done

run_bounded() {
    local seconds="$1"
    local label="$2"
    shift 2
    local status=0
    if command -v timeout >/dev/null 2>&1; then
        timeout --foreground "${seconds}" "$@" || status=$?
    else
        "${PYTHON_BIN}" - "$seconds" "$@" <<'PY' || status=$?
import subprocess, sys
timeout_s = int(sys.argv[1])
cmd = sys.argv[2:]
try:
    proc = subprocess.run(cmd, timeout=timeout_s, check=False)
    sys.exit(proc.returncode)
except subprocess.TimeoutExpired:
    sys.exit(124)
PY
    fi
    if (( status == 124 )); then
        echo "ERROR: wkhtmltopdf install step '${label}' timed out after ${seconds}s" >&2
        exit 124
    fi
    if (( status != 0 )); then
        echo "ERROR: wkhtmltopdf install step '${label}' failed with exit code ${status}" >&2
        exit "$status"
    fi
}

APT_GET_CMD=()
if [[ -n "${WKHTMLTOPDF_APT_GET_BIN:-}" ]]; then
    APT_GET_CMD=("${WKHTMLTOPDF_APT_GET_BIN}")
else
    APT_GET_CMD=(sudo apt-get)
fi

WGET_CMD="${WKHTMLTOPDF_WGET_BIN:-wget}"

DPKG_CMD=()
if [[ -n "${WKHTMLTOPDF_DPKG_BIN:-}" ]]; then
    DPKG_CMD=("${WKHTMLTOPDF_DPKG_BIN}")
else
    DPKG_CMD=(sudo dpkg)
fi

APT_NET_OPTS=(
    -o "Acquire::Retries=${WKHTMLTOPDF_RETRIES}"
    -o "Acquire::http::Timeout=${WKHTMLTOPDF_READ_TIMEOUT}"
    -o "Acquire::https::Timeout=${WKHTMLTOPDF_READ_TIMEOUT}"
)

# Refresh the runner's package index before installing font dependencies.
run_bounded "${WKHTMLTOPDF_APT_TIMEOUT}" "apt-get update" \
    "${APT_GET_CMD[@]}" "${APT_NET_OPTS[@]}" update
run_bounded "${WKHTMLTOPDF_APT_TIMEOUT}" "apt-get install fonts" \
    "${APT_GET_CMD[@]}" "${APT_NET_OPTS[@]}" install -y xfonts-base xfonts-75dpi
run_bounded "${WKHTMLTOPDF_DOWNLOAD_TIMEOUT}" "wget wkhtmltox package" \
    "${WGET_CMD}" \
    --tries="${WKHTMLTOPDF_RETRIES}" \
    --connect-timeout="${WKHTMLTOPDF_CONNECT_TIMEOUT}" \
    --read-timeout="${WKHTMLTOPDF_READ_TIMEOUT}" \
    --timeout="${WKHTMLTOPDF_READ_TIMEOUT}" \
    --waitretry="${WKHTMLTOPDF_WAITRETRY}" \
    --retry-connrefused \
    -O "${WKHTMLTOPDF_DEB_PATH}" \
    "${WKHTMLTOPDF_DEB_URL}"
run_bounded "${WKHTMLTOPDF_DPKG_TIMEOUT}" "dpkg install wkhtmltox" \
    "${DPKG_CMD[@]}" -i "${WKHTMLTOPDF_DEB_PATH}"
