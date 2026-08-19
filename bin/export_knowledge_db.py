#!/usr/bin/env python3
"""
export_knowledge_db.py

Export every base table of the pipeline-lifetime shared PostgreSQL knowledge
warehouse into a single, portable DuckDB file. Run once, gated on every
species' BUILD_KNOWLEDGE_DB call having finished, so PATHOGEN_IDENTIFICATION
can read a self-contained local file instead of holding a live connection to
the shared Postgres server for the rest of the run -- the server can then be
stopped immediately afterwards instead of staying alive until the pipeline's
very end.

No table list is hardcoded: every base table in the `public` schema is
copied generically, so this keeps working if database/knowledge_schema.sql
gains new tables in the future.
"""

import argparse
import json
import sys
from pathlib import Path

import duckdb
import pandas as pd
import psycopg2


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--dbname", default="genomic_intelligence")
    parser.add_argument("--user", default="postgres")
    parser.add_argument("--output", required=True, help="Output .duckdb file path")
    return parser.parse_args()


def list_base_tables(pg_conn) -> list[str]:
    with pg_conn.cursor() as cur:
        cur.execute(
            """
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'public' AND table_type = 'BASE TABLE'
            ORDER BY table_name
            """
        )
        return [row[0] for row in cur.fetchall()]


def main():
    args = parse_args()
    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()

    print(f"Connecting to shared PostgreSQL at {args.host}:{args.port}...", file=sys.stderr)
    pg_conn = psycopg2.connect(
        host=args.host, port=args.port, dbname=args.dbname, user=args.user
    )

    tables = list_base_tables(pg_conn)
    if not tables:
        print("WARNING: no base tables found to export.", file=sys.stderr)

    duck_con = duckdb.connect(str(out_path))
    try:
        for table in tables:
            print(f"Exporting table '{table}'...", file=sys.stderr)
            df = pd.read_sql_query(f'SELECT * FROM "{table}"', pg_conn)
            # PostgreSQL jsonb columns arrive as Python dicts via psycopg2;
            # convert them to JSON text so DuckDB stores a plain VARCHAR
            # that R-side jsonlite::fromJSON can parse reliably.
            for col in df.columns:
                if df[col].dtype == object and len(df) > 0:
                    sample_val = df[col].dropna().iloc[0] if not df[col].dropna().empty else None
                    if isinstance(sample_val, dict):
                        df[col] = df[col].apply(lambda v: json.dumps(v) if isinstance(v, dict) else v)
            # DuckDB's Python API picks up local variables via its
            # replacement scan, so `df` here is directly queryable in SQL.
            duck_con.execute(f'CREATE TABLE "{table}" AS SELECT * FROM df')
        duck_con.commit()
    finally:
        duck_con.close()
        pg_conn.close()

    print(f"Wrote {len(tables)} table(s) to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
