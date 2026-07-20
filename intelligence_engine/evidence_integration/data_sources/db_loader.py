"""Targeted database loader for the Genomic Intelligence Engine.

Instead of bulk-loading whole reference tables, this module requests only the
rows the engine actually needs for a given sample: one lineage record, contextual
genomes for that lineage, protein variants matching detected mutations, and
phenotype associations for detected variants/motifs plus intervention-relevant
categories.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
import sys

sys.path.insert(0, str(PROJECT_ROOT))

try:
    from config import DB_URL as CONFIG_DB_URL
except ImportError:  # pragma: no cover - defensive fallback
    CONFIG_DB_URL = "postgresql://localhost:5432/pgirl"

log = logging.getLogger(__name__)

# Categories needed for negative-finding reports (e.g., "no known vaccine-escape").
INTERVENTION_CATEGORIES = [
    "vaccine_escape",
    "vaccine_effectiveness",
    "drug_resistance",
    "drug_susceptibility",
    "diagnostic_sensitivity",
    "increased_transmission",
]


def _pg_array_literal(values: list[str]) -> str:
    """Render a Python list as a PostgreSQL array literal string."""
    if not values:
        return "{}"
    parts = []
    for v in values:
        s = str(v).strip()
        if not s:
            continue
        if any(c in s for c in ", \"'"):
            escaped = s.replace('"', '""')
            parts.append(f'"{escaped}"')
        else:
            parts.append(s)
    return "{" + ",".join(parts) + "}"


def _arrays_to_strings(df: pd.DataFrame) -> pd.DataFrame:
    """Convert list columns to PostgreSQL array literals for CSV compatibility."""
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, list)).any():
            df[col] = df[col].apply(
                lambda x: _pg_array_literal(x) if isinstance(x, list) else (x if pd.notna(x) else "")
            )
    return df


def _query(db_url: str, query: str, params: tuple) -> pd.DataFrame:
    """Run a SQL query and return a DataFrame."""
    return pd.read_sql_query(query, db_url, params=params)


def _try_query(query_fn, fallback_value: Optional[pd.DataFrame] = None):
    """Wrap a query and return fallback_value on any error."""
    try:
        return query_fn()
    except Exception as exc:
        log.warning("DB query failed: %s", exc)
        return fallback_value


def load_lineage(
    pathogen_id: str,
    species_id: str,
    lineage_id: Optional[str] = None,
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Load a single lineage row, or all lineages for the pathogen/species."""
    url = db_url or CONFIG_DB_URL

    def _query_fn():
        if lineage_id:
            q = "SELECT * FROM lineages WHERE pathogen_id = %s AND species_id = %s AND lineage_id = %s"
            p = (pathogen_id, species_id, lineage_id)
        else:
            q = "SELECT * FROM lineages WHERE pathogen_id = %s AND species_id = %s"
            p = (pathogen_id, species_id)
        df = _query(url, q, p)
        return df

    return _try_query(_query_fn)


def load_genome_metadata_for_lineage(
    pathogen_id: str,
    species_id: str,
    lineage_id: Optional[str] = None,
    limit: Optional[int] = None,
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Load contextual genome metadata for a lineage (or all for the species)."""
    url = db_url or CONFIG_DB_URL

    def _query_fn():
        if lineage_id:
            q = (
                "SELECT * FROM genome_metadata "
                "WHERE pathogen_id = %s AND species_id = %s AND lineage_id = %s"
            )
            p = (pathogen_id, species_id, lineage_id)
        else:
            q = "SELECT * FROM genome_metadata WHERE pathogen_id = %s AND species_id = %s"
            p = (pathogen_id, species_id)
        if limit:
            q += f" LIMIT {int(limit)}"
        df = _query(url, q, p)
        return df

    return _try_query(_query_fn)


def load_protein_variants_for_variants(
    pathogen_id: str,
    species_id: str,
    variants: list[dict],
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Load curated variant records matching the detected variants."""
    url = db_url or CONFIG_DB_URL

    if not variants:
        return None

    genes = sorted({str(v.get("gene", "")).strip().upper() for v in variants if v.get("gene")})
    positions = [v.get("position") for v in variants if v.get("position") is not None]
    alts = sorted({str(v.get("alt_aa", "")).strip().upper() for v in variants if v.get("alt_aa")})

    if not genes:
        return None

    def _query_fn():
        # Query all variants in the detected genes and positions; the analyzer
        # will match the exact alt_aa. This keeps the query small and fast.
        gene_placeholders = ", ".join(["%s"] * len(genes))
        params: list[Any] = [pathogen_id, species_id]
        q = (
            "SELECT * FROM v_variant_summary "
            "WHERE pathogen_id = %s AND species_id = %s "
            f"AND upper(gene) IN ({gene_placeholders})"
        )
        params.extend(genes)

        if positions:
            pos_placeholders = ", ".join(["%s"] * len(positions))
            q += f" AND position IN ({pos_placeholders})"
            params.extend(positions)
        if alts:
            alt_placeholders = ", ".join(["%s"] * len(alts))
            q += f" AND upper(alt_aa) IN ({alt_placeholders})"
            params.extend(alts)

        df = _query(url, q, tuple(params))
        return df

    return _try_query(_query_fn)


def load_genotype_phenotype_for_variants(
    pathogen_id: str,
    species_id: str,
    variants: list[dict],
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Load phenotype associations for detected variants, motifs, and key categories."""
    url = db_url or CONFIG_DB_URL

    genes = sorted({str(v.get("gene", "")).strip().upper() for v in variants if v.get("gene")})
    positions = [v.get("position") for v in variants if v.get("position") is not None]
    alts = sorted({str(v.get("alt_aa", "")).strip().upper() for v in variants if v.get("alt_aa")})

    if not genes:
        return None

    def _query_fn():
        params: list[Any] = [pathogen_id, species_id]
        q = (
            "SELECT * FROM genotype_phenotype "
            "WHERE pathogen_id = %s AND species_id = %s AND ("
        )

        parts = []

        # 1. Exact amino-acid change associations in detected genes.
        if positions and alts:
            gene_placeholders = ", ".join(["%s"] * len(genes))
            pos_placeholders = ", ".join(["%s"] * len(positions))
            alt_placeholders = ", ".join(["%s"] * len(alts))
            parts.append(
                f"(upper(protein) IN ({gene_placeholders}) "
                f"AND position IN ({pos_placeholders}) "
                f"AND upper(alt_aa) IN ({alt_placeholders}))"
            )
            params.extend(genes)
            params.extend(positions)
            params.extend(alts)

        # 2. Motif-level associations for detected genes (position is NULL).
        gene_placeholders = ", ".join(["%s"] * len(genes))
        parts.append(
            f"(upper(protein) IN ({gene_placeholders}) AND position IS NULL)"
        )
        params.extend(genes)

        # 3. Intervention categories needed for negative findings.
        if INTERVENTION_CATEGORIES:
            cat_placeholders = ", ".join(["%s"] * len(INTERVENTION_CATEGORIES))
            parts.append(f"(phenotype_category IN ({cat_placeholders}))")
            params.extend(INTERVENTION_CATEGORIES)

        q += " OR ".join(parts) + ")"

        df = _query(url, q, tuple(params))
        return df

    return _try_query(_query_fn)


# Keep broad loaders as a fallback when no sample-specific filters are available.

def load_table(
    table_or_view: str,
    pathogen_id: str,
    species_id: str,
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Load a full reference table/view filtered by pathogen and species."""
    url = db_url or CONFIG_DB_URL

    def _query_fn():
        q = f"SELECT * FROM {table_or_view} WHERE pathogen_id = %s AND species_id = %s"
        df = _query(url, q, (pathogen_id, species_id))
        return df

    return _try_query(_query_fn)


def load_genotype_phenotype(
    pathogen_id: str,
    species_id: str,
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    return load_table("genotype_phenotype", pathogen_id, species_id, db_url)


def load_protein_variants(
    pathogen_id: str,
    species_id: str,
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    return load_table("v_variant_summary", pathogen_id, species_id, db_url)


def load_lineages(
    pathogen_id: str,
    species_id: str,
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    return load_table("lineages", pathogen_id, species_id, db_url)


def load_genome_metadata(
    pathogen_id: str,
    species_id: str,
    db_url: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    return load_table("genome_metadata", pathogen_id, species_id, db_url)
