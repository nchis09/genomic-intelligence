#!/usr/bin/env python3
"""
start_shared_db.py

Start (or reuse, if already running) the pipeline-lifetime shared PostgreSQL
instance used by BUILD_KNOWLEDGE_DB (all species) and QUERY_KNOWLEDGE_DB
(FIGURES stage). Runs once, early in the pipeline, from a Nextflow process
whose task can complete immediately afterwards -- `pg_ctl start` forks the
postmaster as an independent daemon bound to a fixed, non-work-dir data
directory and host:port, so it keeps running after this task finishes.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_knowledge_db import TemporaryPostgres  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, help="Fixed, persistent Postgres data directory (outside any task work dir)")
    parser.add_argument("--log-file", required=True, help="Postgres server log file")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    data_dir = Path(args.data_dir).resolve()
    data_dir.parent.mkdir(parents=True, exist_ok=True)

    pg = TemporaryPostgres(data_dir, log_file=args.log_file, host=args.host, port=args.port)

    print(f"Initializing shared PostgreSQL data directory at {data_dir}...", file=sys.stderr)
    pg.init(reset=False)
    print(f"Starting shared PostgreSQL on {args.host}:{args.port}...", file=sys.stderr)
    pg.start()
    print("Creating warehouse database (if needed)...", file=sys.stderr)
    pg.create_db()
    print("Shared PostgreSQL is up.", file=sys.stderr)


if __name__ == "__main__":
    main()
