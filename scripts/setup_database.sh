#!/usr/bin/env bash
# One-time database setup for the PGIRL pipeline.
#
# Creates the PostgreSQL database, loads the schema, fetches the canonical
# Ebola reference proteomes, and syncs curated reference data.
# It is safe to run multiple times; existing tables are not dropped.
#
# Usage:
#   PGIRL_DB_URL="postgresql://localhost:5432/pgirl" ./scripts/setup_database.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

DB_URL="${PGIRL_DB_URL:-postgresql://localhost:5432/pgirl}"
DB_NAME="${DB_URL##*/}"

export PYTHONPATH="${PROJECT_ROOT}"

# Prefer the active Conda environment's Python, then common installations,
# then the system python3. This makes the script work across machines.
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/python3" ]; then
    PYTHON="${CONDA_PREFIX}/bin/python3"
elif [ -x "${HOME}/anaconda3/bin/python3" ]; then
    PYTHON="${HOME}/anaconda3/bin/python3"
elif [ -x "${HOME}/miniconda3/bin/python3" ]; then
    PYTHON="${HOME}/miniconda3/bin/python3"
elif [ -x "/opt/miniconda3/bin/python3" ]; then
    PYTHON="/opt/miniconda3/bin/python3"
else
    PYTHON="python3"
fi

# Find psql binary (common Homebrew path + PATH)
if command -v /opt/homebrew/opt/postgresql@16/bin/psql >/dev/null 2>&1; then
    PSQL="/opt/homebrew/opt/postgresql@16/bin/psql"
elif command -v psql >/dev/null 2>&1; then
    PSQL="psql"
else
    echo "ERROR: psql not found. Install PostgreSQL first." >&2
    exit 1
fi

# Find createdb binary (common Homebrew path + PATH)
if command -v /opt/homebrew/opt/postgresql@16/bin/createdb >/dev/null 2>&1; then
    CREATEDB="/opt/homebrew/opt/postgresql@16/bin/createdb"
elif command -v createdb >/dev/null 2>&1; then
    CREATEDB="createdb"
else
    echo "ERROR: createdb not found. Install PostgreSQL first." >&2
    exit 1
fi

echo "=== PGIRL database setup ==="
echo "Database URL: ${DB_URL}"
echo "Python:       ${PYTHON}"
echo "psql:         ${PSQL}"

# Create database if it does not exist
if ! "${PSQL}" -d "${DB_NAME}" -c "SELECT 1" >/dev/null 2>&1; then
    echo "Creating database ${DB_NAME} ..."
    "${CREATEDB}" "${DB_NAME}"
else
    echo "Database ${DB_NAME} already exists."
fi

# Load schema
schema_file="${PROJECT_ROOT}/database/schema.sql"
if [ ! -f "${schema_file}" ]; then
    echo "ERROR: schema file not found: ${schema_file}" >&2
    exit 1
fi
echo "Loading schema from ${schema_file} ..."
"${PSQL}" -d "${DB_NAME}" -f "${schema_file}"

# Fetch reference genomes / proteomes from NCBI
echo "Fetching Ebola reference proteomes ..."
"${PYTHON}" "${PROJECT_ROOT}/database/ebola/protein_variants/fetch_reference_proteomes.py"

# Optionally sync curated reference data from external sources.
# This step can take 10-30 minutes because it fetches data from NCBI, UniProt,
# and PubTator. Set PGIRL_SYNC_SOURCES=none to skip it and only load the schema
# and canonical reference proteomes (enough to run the deterministic pipeline).
SYNC_SOURCES="${PGIRL_SYNC_SOURCES:-all}"
if [ "${SYNC_SOURCES}" = "none" ]; then
    echo "Skipping curated data sync (PGIRL_SYNC_SOURCES=none)."
else
    echo "Syncing curated reference data (sources: ${SYNC_SOURCES}) ..."
    "${PYTHON}" "${PROJECT_ROOT}/scripts/db_sync.py" \
        --db-url "${DB_URL}" \
        --pathogens ebola \
        --sources ${SYNC_SOURCES}
fi

echo "=== PGIRL database setup complete ==="
