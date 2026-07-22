"""
config.py — Central configuration for the Genomic Epidemic Intelligence System.

This is the ONLY file you need to edit when moving to a new machine or location.
All other scripts import from here — no hardcoded paths anywhere else.

To move to a new computer:
  1. Copy the whole project folder
  2. Install PostgreSQL and restore from backup:
         psql pgirl < database/pgirl_backup.sql
  3. Update DB_URL below if needed (usually unchanged)
  4. Done — everything else works automatically
"""

import os
from pathlib import Path

# ─── Base directory ───────────────────────────────────────────────────────────
# Automatically resolves to wherever this project folder lives.
# No manual editing needed when moving machines — it's always relative to this file.
BASE_DIR = Path(__file__).parent.resolve()

# ─── Sub-directories ──────────────────────────────────────────────────────────
DATABASE_DIR    = BASE_DIR / "database"
API_SOURCES_DIR = DATABASE_DIR / "api_sources"
EBOLA_DATA_DIR  = DATABASE_DIR / "ebola"
SCRIPTS_DIR     = BASE_DIR / "scripts"
INTELLIGENCE_DIR = BASE_DIR / "intelligence_engine"

# ─── Bioinformatics tool databases ────────────────────────────────────────────
# Kraken2 viral database for taxonomic classification.
# Override by setting environment variable PGIRL_KRAKEN_DB before running any script.
KRAKEN_DB = os.environ.get("PGIRL_KRAKEN_DB", str(Path.home() / "dabase" / "kraken_database"))

# ─── Database connection ──────────────────────────────────────────────────────
# Override by setting environment variable PGIRL_DB_URL before running any script.
# Example: export PGIRL_DB_URL=postgresql://user:password@remotehost:5432/pgirl
DB_URL = os.environ.get("PGIRL_DB_URL", "postgresql://localhost:5432/pgirl")

# ─── Ebola curated data files ─────────────────────────────────────────────────
EBOLA_OUTBREAKS_CSV         = EBOLA_DATA_DIR / "outbreaks" / "outbreaks.csv"
EBOLA_MUTATIONS_CSV         = EBOLA_DATA_DIR / "mutations" / "mutation_catalogue.csv"
EBOLA_MOLECULAR_EPI_CSV     = EBOLA_DATA_DIR / "molecular_epidemiology" / "molecular_epidemiology.csv"
EBOLA_GENOTYPE_PHENOTYPE_CSV = EBOLA_DATA_DIR / "genotype_phenotype" / "genotype_phenotype_associations.csv"

# ─── LLM configuration ────────────────────────────────────────────────────────
# Only Ollama is supported for privacy — no data leaves the machine.
# To use Ollama:
#   curl -fsSL https://ollama.com/install.sh | sh
#   ollama pull qwen2.5:7b
#   (or pull qwen2.5:14b for higher-quality but slower narrative)
#
# Override any of these via environment variables.
LLM_PROVIDER   = os.environ.get("PGIRL_LLM_PROVIDER", "ollama")  # ollama only
OLLAMA_URL     = os.environ.get("OLLAMA_URL", "http://localhost:11434")
OLLAMA_MODEL   = os.environ.get("OLLAMA_MODEL", "qwen2.5:7b")
