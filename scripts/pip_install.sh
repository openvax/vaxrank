#!/usr/bin/env bash
# Adapted from openvax/topiary scripts/pip_install.sh (#427).
# Retry CI dependency resolution without relaxing any supplied requirement.
set -uo pipefail

attempts="${PIP_INSTALL_ATTEMPTS:-3}"
delay="${PIP_INSTALL_RETRY_DELAY:-20}"
python="${PYTHON:-python}"
if [[ ! "$attempts" =~ ^[1-9][0-9]*$ || ! "$delay" =~ ^(0|[1-9][0-9]*)$ ]]; then
    echo "PIP_INSTALL_ATTEMPTS must be positive and PIP_INSTALL_RETRY_DELAY nonnegative integers." >&2
    exit 2
fi

for (( attempt = 1; attempt <= attempts; attempt++ )); do
    install_args=(-m pip install)
    if (( attempt > 1 )); then
        install_args+=(--no-cache-dir)
    fi
    if "$python" "${install_args[@]}" "$@"; then
        exit 0
    else
        status=$?
    fi
    # Do not retry an interrupted command or a missing/nonexecutable interpreter.
    if (( status == 126 || status == 127 || status >= 128 )); then
        exit "$status"
    fi
    if (( attempt < attempts )); then
        pause=$(( attempt * delay ))
        echo "pip install failed (attempt ${attempt}/${attempts}); retrying in ${pause}s without pip's cache" >&2
        sleep "$pause" || exit $?
    fi
done
echo "pip install failed after ${attempts} attempts" >&2
exit "$status"
