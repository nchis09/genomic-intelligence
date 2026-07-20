"""storage.py — Persist the normalized epi master object as Parquet + DuckDB.

Writes each entity list (outbreaks, molecular_epidemiology, demographics,
clinical, interventions, diagnostics, therapeutics, vaccines, surveillance,
genomic_links, knowledge_assertions, references) to its own Parquet file
under ``<output_dir>/tables/``. Singleton sections (metadata, pathogen_profile,
transmission) are written as single-row tables. Provides a helper to open a
DuckDB connection with all tables registered as views for ad-hoc SQL analysis.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pyarrow as pa
import pyarrow.parquet as pq

from intelligence_engine.data_engine.analytics.schemas import (
    ClinicalFeature,
    DemographicSummary,
    DiagnosticMethod,
    EpiMetadata,
    GenomicLink,
    InterventionRecord,
    KnowledgeAssertion,
    MolecularEpidemiology,
    NormalizedEpiObject,
    OutbreakRecord,
    PathogenProfile,
    Reference,
    SurveillanceSummary,
    TherapeuticProduct,
    TransmissionParams,
    VaccineProduct,
)

log = logging.getLogger(__name__)

# List-valued entity sections
_LIST_TABLE_MODELS = {
    "molecular_epidemiology": MolecularEpidemiology,
    "outbreaks": OutbreakRecord,
    "demographics": DemographicSummary,
    "clinical": ClinicalFeature,
    "interventions": InterventionRecord,
    "diagnostics": DiagnosticMethod,
    "therapeutics": TherapeuticProduct,
    "vaccines": VaccineProduct,
    "surveillance": SurveillanceSummary,
    "genomic_links": GenomicLink,
    "knowledge_assertions": KnowledgeAssertion,
    "references": Reference,
}

# Singleton entity sections (0 or 1 row)
_SINGLETON_TABLE_MODELS = {
    "metadata": EpiMetadata,
    "pathogen_profile": PathogenProfile,
    "transmission": TransmissionParams,
}

_TABLE_NAMES = tuple(_LIST_TABLE_MODELS) + tuple(_SINGLETON_TABLE_MODELS)

# Python type -> Arrow type used to build an empty schema when a table has
# zero rows (Arrow/DuckDB cannot read a Parquet file with no columns).
_PY_TO_ARROW = {
    str: pa.string(),
    int: pa.int64(),
    float: pa.float64(),
    bool: pa.bool_(),
}


def _empty_table_for_model(model_cls) -> pa.Table:
    """Build a zero-row Arrow table with one column per model field."""
    fields = []
    for name, info in model_cls.model_fields.items():
        annotation = info.annotation
        # Unwrap Optional[X] -> X
        py_type = getattr(annotation, "__args__", (annotation,))[0] if hasattr(annotation, "__args__") else annotation
        arrow_type = _PY_TO_ARROW.get(py_type, pa.string())
        fields.append(pa.field(name, arrow_type))
    schema = pa.schema(fields)
    return pa.Table.from_pylist([], schema=schema)


def write_parquet_tables(dataset: NormalizedEpiObject, output_dir: str) -> dict[str, str]:
    """Write each entity section of *dataset* to a Parquet file under
    output_dir/tables/.

    List-valued sections (outbreaks, demographics, etc.) become one row per
    list item. Singleton sections (metadata, pathogen_profile, transmission)
    become a 0-or-1-row table.

    Returns:
        Dict mapping table name to the written file path.
    """
    tables_dir = Path(output_dir) / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    written = {}

    for table_name, model_cls in _LIST_TABLE_MODELS.items():
        records = [r.model_dump() for r in getattr(dataset, table_name)]
        path = tables_dir / f"{table_name}.parquet"
        if not records:
            arrow_table = _empty_table_for_model(model_cls)
        else:
            arrow_table = pa.Table.from_pylist(records)
        pq.write_table(arrow_table, path)
        written[table_name] = str(path)
        log.info(f"Wrote {len(records)} rows to {path}")

    for table_name, model_cls in _SINGLETON_TABLE_MODELS.items():
        value = getattr(dataset, table_name)
        path = tables_dir / f"{table_name}.parquet"
        if value is None:
            arrow_table = _empty_table_for_model(model_cls)
            row_count = 0
        else:
            arrow_table = pa.Table.from_pylist([value.model_dump()])
            row_count = 1
        pq.write_table(arrow_table, path)
        written[table_name] = str(path)
        log.info(f"Wrote {row_count} row(s) to {path}")

    return written


def open_duckdb(output_dir: str, db_path: Optional[str] = None):
    """Open a DuckDB connection with all normalized tables registered as views.

    Args:
        output_dir: Directory containing tables/<name>.parquet files
                    (as written by write_parquet_tables).
        db_path: Optional path for a persistent DuckDB file. If None, uses
                 an in-memory database.

    Returns:
        A duckdb.DuckDBPyConnection with views: outbreaks, transmission_params,
        interventions, demographics, genomic_links, references.
    """
    import duckdb

    con = duckdb.connect(db_path or ":memory:")
    tables_dir = Path(output_dir) / "tables"
    for table_name in _TABLE_NAMES:
        parquet_path = tables_dir / f"{table_name}.parquet"
        if parquet_path.exists():
            con.execute(
                f'CREATE OR REPLACE VIEW "{table_name}" AS SELECT * FROM read_parquet(\'{parquet_path}\')'
            )
    return con
