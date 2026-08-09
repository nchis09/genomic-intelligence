#!/usr/bin/env python3
"""
stop_shared_db.py

Dump the pipeline-lifetime shared PostgreSQL database to a single, portable
.sql file, then stop the server. Run once, gated on the completion of every
species' BUILD_KNOWLEDGE_DB and the FIGURES stage's QUERY_KNOWLEDGE_DB calls.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_knowledge_db import TemporaryPostgres, dump_database  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, help="Same data directory passed to start_shared_db.py")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--outdir", required=True, help="Directory to write the final .sql dump into")
    parser.add_argument("--prefix", default="genomic_intelligence", help="Dump file prefix")
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    pg = TemporaryPostgres(Path(args.data_dir).resolve(), host=args.host, port=args.port)

    print("Dumping shared database...", file=sys.stderr)
    dump_database(pg, outdir, args.prefix)

    print("Stopping shared PostgreSQL...", file=sys.stderr)
    pg.stop()
    print("Shared PostgreSQL stopped.", file=sys.stderr)


if __name__ == "__main__":
    main()
