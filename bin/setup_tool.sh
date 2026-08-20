#!/usr/bin/env bash
set -euo pipefail

TOOL_NAME="${1}"
REPO_URL="${2}"

# Resolve projectDir/tools relative to this script (bin/ -> ../tools)
TOOLS_ROOT="$(cd "$(dirname "$0")/../tools" && pwd)"
TOOL_DIR="${TOOLS_ROOT}/${TOOL_NAME}"
LOCK_DIR="${TOOLS_ROOT}/.${TOOL_NAME}.setup.lock"

mkdir -p "${TOOLS_ROOT}"

acquire_lock() {
    while true; do
        # Try to create the lock directory atomically
        if mkdir "${LOCK_DIR}" 2>/dev/null; then
            return 0
        fi
        # If another process finished the setup while we waited, no need to clone
        if [ -f "${TOOL_DIR}/DESCRIPTION" ]; then
            return 1
        fi
        echo "Waiting for ${TOOL_NAME} setup..." >&2
        sleep 5
    done
}

if [ ! -f "${TOOL_DIR}/DESCRIPTION" ]; then
    if acquire_lock; then
        # Second check after acquiring the lock
        if [ ! -f "${TOOL_DIR}/DESCRIPTION" ]; then
            # Make writable in case a previous clone left read-only .git objects
            chmod -R +w "${TOOL_DIR}" 2>/dev/null || true
            rm -rf "${TOOL_DIR}"

            TMP_DIR="$(mktemp -d)"
            git clone --depth 1 "${REPO_URL}" "${TMP_DIR}/${TOOL_NAME}"
            cp -a "${TMP_DIR}/${TOOL_NAME}" "${TOOL_DIR}"
            rm -rf "${TMP_DIR}"
        fi
        # Release lock
        rmdir "${LOCK_DIR}" 2>/dev/null || true
    fi
fi

if [ ! -f "${TOOL_DIR}/DESCRIPTION" ]; then
    echo "ERROR: ${TOOL_NAME} was not set up in ${TOOL_DIR}" >&2
    exit 1
fi

echo "${TOOL_DIR}"
