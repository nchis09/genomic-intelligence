#!/usr/bin/env python3
"""
run_schemaspy.py

Post-pipeline CLI that generates an interactive HTML schema report for the
PostgreSQL knowledge warehouse using SchemaSpy.

Requires:
  - Java (openjdk)
  - SchemaSpy JAR (e.g. schemaspy-*.jar)
  - PostgreSQL JDBC driver JAR (e.g. postgresql-*.jar)

These JARs can be provided via:
  1. --schemaspy-jar and --jdbc-jar command-line options
  2. SCHEMASPY_JAR and PGJDBC_JAR environment variables
  3. assets/schemaspy/ directory in the repository

Example:
    conda run -n pgirl_schemaspy python bin/run_schemaspy.py --outdir results
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import psycopg2

# Import the shared TemporaryPostgres helper from build_knowledge_db.py
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))
from build_knowledge_db import TemporaryPostgres


def find_jar(name_pattern, env_var, default_dir, description):
    """Locate a JAR file using env var, default directory, or CLI."""
    if env_var and os.environ.get(env_var):
        return Path(os.environ[env_var]).resolve()

    default_path = REPO_ROOT / default_dir
    if default_path.is_dir():
        matches = sorted(default_path.glob(name_pattern))
        if matches:
            return matches[0]

    return None


def run_schemaspy(args):
    outdir = Path(args.outdir).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else outdir / "pipeline_info" / "schemaspy"
    output_dir.mkdir(parents=True, exist_ok=True)

    data_dir = outdir / "knowledge_warehouse" / "_shared_pg_data"
    if not data_dir.exists():
        print(
            f"ERROR: knowledge warehouse data directory not found: {data_dir}\n"
            "Run the pipeline first, or provide --db-host/--db-port to connect to an existing DB.",
            file=sys.stderr,
        )
        return 1

    pg = None
    started_here = False
    try:
        if args.db_host and args.db_port:
            host = args.db_host
            port = args.db_port
        else:
            log_file = outdir / "knowledge_warehouse" / "_schemaspy_postgres.log"
            pg = TemporaryPostgres(data_dir, log_file=log_file, host="127.0.0.1", port=None)
            pg.start()
            started_here = True
            host = pg.host
            port = pg.port

        # Verify connection
        conn = psycopg2.connect(host=host, port=port, dbname=pg.db_name if pg else "genomic_intelligence", user="postgres")
        conn.close()

        schemaspy_jar = Path(args.schemaspy_jar).resolve() if args.schemaspy_jar else find_jar(
            "schemaspy*.jar", "SCHEMASPY_JAR", "assets/schemaspy", "SchemaSpy JAR"
        )
        jdbc_jar = Path(args.jdbc_jar).resolve() if args.jdbc_jar else find_jar(
            "postgresql*.jar", "PGJDBC_JAR", "assets/schemaspy", "PostgreSQL JDBC driver JAR"
        )

        missing = []
        if not schemaspy_jar or not schemaspy_jar.exists():
            missing.append("SchemaSpy JAR (--schemaspy-jar or SCHEMASPY_JAR)")
        if not jdbc_jar or not jdbc_jar.exists():
            missing.append("PostgreSQL JDBC driver JAR (--jdbc-jar or PGJDBC_JAR)")
        if missing:
            print(f"ERROR: missing required JARs: {', '.join(missing)}\n", file=sys.stderr)
            print(
                "Place them in assets/schemaspy/ or set the environment variables:\n"
                "  SCHEMASPY_JAR=/path/to/schemaspy-x.y.z.jar\n"
                "  PGJDBC_JAR=/path/to/postgresql-x.y.z.jar\n",
                file=sys.stderr,
            )
            return 1

        cmd = [
            "java", "-jar", str(schemaspy_jar),
            "-t", "pgsql",
            "-dp", str(jdbc_jar),
            "-db", "genomic_intelligence",
            "-s", "public",
            "-host", host,
            "-port", str(port),
            "-u", "postgres",
            "-o", str(output_dir),
        ]
        print(f"Running SchemaSpy...\n  output: {output_dir}", file=sys.stderr)
        subprocess.run(cmd, check=True)
        print(f"SchemaSpy report generated: {output_dir / 'index.html'}", file=sys.stderr)
        return 0

    except subprocess.CalledProcessError as e:
        print(f"ERROR: SchemaSpy failed (exit {e.returncode}): {e}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    finally:
        if started_here and pg is not None:
            pg.stop()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate an interactive SchemaSpy HTML report for the knowledge warehouse.",
    )
    parser.add_argument(
        "--outdir",
        default="results",
        help="Pipeline output directory (default: results)",
    )
    parser.add_argument(
        "--output-dir",
        help="Destination for the SchemaSpy report (default: {outdir}/pipeline_info/schemaspy)",
    )
    parser.add_argument(
        "--db-host",
        help="Existing PostgreSQL host (if not set, a temporary server is started from {outdir}/knowledge_warehouse/_shared_pg_data)",
    )
    parser.add_argument(
        "--db-port",
        type=int,
        help="Existing PostgreSQL port (required with --db-host)",
    )
    parser.add_argument(
        "--schemaspy-jar",
        help="Path to the SchemaSpy JAR",
    )
    parser.add_argument(
        "--jdbc-jar",
        help="Path to the PostgreSQL JDBC driver JAR",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run_schemaspy(parse_args()))
