#!/usr/bin/env python3
"""
build_knowledge_db.py

Build a temporary PostgreSQL knowledge warehouse for the Genomic Intelligence pipeline.

The script:
  1. Initializes a local PostgreSQL data directory.
  2. Starts a temporary server on a free local port.
  3. Creates the database and loads the schema.
  4. Registers every pipeline output file in pipeline_outputs.
  5. Ingests reference genomes, samples, mutations, trees, phenotype annotations,
     and epidemiological data.
  6. Dumps the populated database to a portable .sql file.
  7. Stops the server.

All state is written under --outdir, so the whole directory can be deleted.
"""

import argparse
import csv
import glob
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

import psycopg2


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_free_port(host="127.0.0.1"):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((host, 0))
        return s.getsockname()[1]


def run_cmd(cmd, check=True, env=None, **kwargs):
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    return subprocess.run(cmd, check=check, env=merged_env, **kwargs)


def find_executable(name):
    # Prefer the active conda environment so we don't accidentally pick a
    # system or Homebrew psql/pg_ctl that is incompatible with the server
    # we just started from the pipeline's own conda packages.
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate = os.path.join(conda_prefix, "bin", name)
        if os.path.isfile(candidate):
            return candidate
    exe = shutil.which(name)
    if exe:
        return exe
    raise FileNotFoundError(f"Required PostgreSQL binary not found: {name}")


class TemporaryPostgres:
    def __init__(self, data_dir, log_file=None, host="127.0.0.1", port=None):
        self.data_dir = Path(data_dir).resolve()
        self.log_file = Path(log_file or self.data_dir.parent / "postgres.log").resolve()
        self.host = host
        self.port = port or get_free_port(host)
        self.db_name = "genomic_intelligence"
        self.superuser = "postgres"

    def is_running(self):
        try:
            with socket.create_connection((self.host, self.port), timeout=2):
                return True
        except OSError:
            return False

    def init(self, reset=True):
        """Initialize the data directory. `reset=False` makes this idempotent
        for a shared, pipeline-lifetime instance: if the directory already
        looks initialized (has PG_VERSION), leave it alone."""
        if self.data_dir.exists():
            if not reset:
                if (self.data_dir / "PG_VERSION").exists():
                    return
                shutil.rmtree(self.data_dir)
            else:
                shutil.rmtree(self.data_dir)
        initdb = find_executable("initdb")
        run_cmd([
            initdb,
            "-D", str(self.data_dir),
            "-U", self.superuser,
            "--no-locale",
            "--encoding=UTF8",
            "--auth=trust",
        ], check=True)

    def start(self):
        """Start the server. If the port is occupied by a different
        (stale) postgres, kill it. If it's our own postgres, just reuse it."""
        pg_ctl = find_executable("pg_ctl")
        if self.is_running():
            status = run_cmd([
                pg_ctl, "status", "-D", str(self.data_dir)
            ], check=False, capture_output=True, text=True)
            if status.returncode == 0 and "is running" in status.stdout:
                return
            # Something else is holding the port; kill it.
            self._kill_stale_listener()
            # Make sure our own data dir is not registered as running.
            run_cmd([pg_ctl, "stop", "-D", str(self.data_dir), "-m", "fast"], check=False)
        options = f"-p {self.port} -h {self.host} -k /tmp"
        run_cmd([
            pg_ctl,
            "-w",  # wait for the server to start (or fail)
            "start",
            "-D", str(self.data_dir),
            "-l", str(self.log_file),
            "-o", options,
        ], check=True)
        self._wait_for_server()

    def _kill_stale_listener(self):
        """Kill any process currently listening on our TCP port."""
        try:
            result = run_cmd(["lsof", "-ti", f"TCP:{self.port}"],
                             check=False, capture_output=True, text=True)
        except FileNotFoundError:
            return
        pids = [p for p in result.stdout.strip().splitlines() if p.isdigit()]
        for pid in pids:
            run_cmd(["kill", "-9", pid], check=False)

    def _wait_for_server(self, timeout=60):
        """Wait for the server to accept SQL connections, not just a TCP socket."""
        psql = find_executable("psql")
        env = {"PGUSER": self.superuser}
        deadline = time.time() + timeout
        while time.time() < deadline:
            if self.is_running():
                try:
                    run_cmd([
                        psql, "-h", self.host, "-p", str(self.port), "-U", self.superuser,
                        "-d", "template1", "-tAc", "SELECT 1",
                    ], check=True, env=env, capture_output=True, text=True, timeout=5)
                    return
                except subprocess.SubprocessError:
                    pass
            time.sleep(0.5)
        raise RuntimeError("PostgreSQL server did not become ready within timeout")

    def stop(self):
        pg_ctl = find_executable("pg_ctl")
        run_cmd([pg_ctl, "stop", "-D", str(self.data_dir), "-m", "fast"], check=False)

    def create_db(self):
        """Create the warehouse database if it doesn't already exist."""
        env = {"PGUSER": self.superuser}
        psql = find_executable("psql")
        for _ in range(20):
            try:
                check = run_cmd([
                    psql, "-h", self.host, "-p", str(self.port), "-U", self.superuser,
                    "-d", "template1",
                    "-tAc", f"SELECT 1 FROM pg_database WHERE datname = '{self.db_name}'",
                ], check=True, env=env, capture_output=True, text=True)
                if check.stdout.strip() == "1":
                    return
                break
            except subprocess.CalledProcessError:
                time.sleep(0.5)
        else:
            raise RuntimeError("PostgreSQL server did not accept SQL before create_db")
        createdb = find_executable("createdb")
        run_cmd([
            createdb,
            "-h", self.host,
            "-p", str(self.port),
            "-U", self.superuser,
            self.db_name,
        ], check=True, env=env)

    def connect(self):
        return psycopg2.connect(
            host=self.host,
            port=self.port,
            dbname=self.db_name,
            user=self.superuser,
        )


def with_advisory_lock(conn, key, fn):
    """Run `fn()` while holding a Postgres session-level advisory lock,
    serializing concurrent callers (e.g. multiple species' BUILD_KNOWLEDGE_DB
    processes racing to CREATE TABLE against the same shared database)."""
    with conn.cursor() as cur:
        cur.execute("SELECT pg_advisory_lock(%s)", (key,))
    try:
        return fn()
    finally:
        with conn.cursor() as cur:
            cur.execute("SELECT pg_advisory_unlock(%s)", (key,))


# Fixed advisory-lock key used to serialize schema creation across species
# processes connecting to the shared knowledge-warehouse database.
SCHEMA_LOCK_KEY = 927364501


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def load_schema(conn, schema_path):
    with open(schema_path) as f:
        sql = f.read()
    with conn.cursor() as cur:
        cur.execute(sql)
    conn.commit()


def load_views(conn, views_path):
    if not views_path.exists():
        print(f"  SKIP: views file not found: {views_path}", file=sys.stderr)
        return
    with open(views_path) as f:
        sql = f.read()
    if not sql.strip():
        return
    with conn.cursor() as cur:
        cur.execute(sql)
    conn.commit()


def create_run(conn, run_id):
    with conn.cursor() as cur:
        cur.execute(
            "INSERT INTO analysis_runs (run_id) VALUES (%s) ON CONFLICT DO NOTHING RETURNING run_id",
            (run_id,),
        )
    conn.commit()


def normalize_date(value):
    if not value or str(value).strip() in ("", "NA", "N/A", "None", "null"):
        return None
    value = str(value).strip()
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return value
    formats = ["%Y-%m-%d", "%Y-%m", "%Y", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d"]
    for fmt in formats:
        try:
            return time.strftime("%Y-%m-%d", time.strptime(value, fmt))
        except ValueError:
            continue
    return None


def normalize_text(value):
    if value is None:
        return None
    value = str(value).strip()
    if value.lower() in ("", "na", "n/a", "none", "null"):
        return None
    return value


def parse_int(value):
    if value is None:
        return None
    try:
        return int(float(str(value).strip()))
    except (ValueError, TypeError):
        return None


def parse_float(value):
    if value is None:
        return None
    try:
        return float(str(value).strip().replace("%", ""))
    except (ValueError, TypeError):
        return None


def normalize_country(value):
    value = normalize_text(value)
    if value and value.lower() in ("unknown", "na", "n/a"):
        return None
    return value


def detect_delimiter(first_line):
    if not first_line:
        return ","
    if "\t" in first_line:
        return "\t"
    if first_line.count(";") > first_line.count(","):
        return ";"
    return ","


def read_tsv_or_csv(path):
    path = Path(path)
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        first = f.readline()
        f.seek(0)
        delimiter = detect_delimiter(first)
        reader = csv.DictReader(f, delimiter=delimiter)
        return list(reader)


def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def count_text_rows(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
        # skip header if present
        if lines and any(c.isalpha() for c in lines[0][:20]):
            return max(0, len([l for l in lines[1:] if l.strip()]))
        return len([l for l in lines if l.strip()])
    except Exception:
        return None


def infer_stage(rel_path):
    parts = [p.lower() for p in Path(rel_path).parts]
    if "classification" in parts:
        return "classification"
    if "nextclade" in parts:
        return "nextclade"
    if "nextstrain" in parts or "bioinformatics_analysis" in parts:
        return "bioinformatics"
    if "phenotype_annotation" in parts or "uniprotr" in parts or "uniprotextractr" in parts or "rbioapi" in parts:
        return "phenotype_annotation"
    if "epidemiological_data" in parts or "epi_data" in parts:
        return "epidemiological_data"
    if "figures" in parts:
        return "reporting"
    if "pipeline_info" in parts:
        return "reporting"
    return "unknown"


def get_or_create_location(cur, country, admin1=None, admin2=None, locality=None, location_code=None, location_code_type=None):
    country = normalize_country(country)
    admin1 = normalize_text(admin1)
    admin2 = normalize_text(admin2)
    locality = normalize_text(locality)
    location_code = normalize_text(location_code)
    location_code_type = normalize_text(location_code_type)
    if not country and not location_code:
        return None
    cur.execute(
        """
        SELECT location_id FROM geographic_locations
        WHERE country IS NOT DISTINCT FROM %s
          AND admin1 IS NOT DISTINCT FROM %s
          AND admin2 IS NOT DISTINCT FROM %s
          AND locality IS NOT DISTINCT FROM %s
          AND location_code IS NOT DISTINCT FROM %s
          AND location_code_type IS NOT DISTINCT FROM %s
        """,
        (country, admin1, admin2, locality, location_code, location_code_type),
    )
    row = cur.fetchone()
    if row:
        return row[0]
    cur.execute(
        """
        INSERT INTO geographic_locations (country, admin1, admin2, locality, location_code, location_code_type)
        VALUES (%s, %s, %s, %s, %s, %s)
        RETURNING location_id
        """,
        (country, admin1, admin2, locality, location_code, location_code_type),
    )
    return cur.fetchone()[0]


# ---------------------------------------------------------------------------
# Provenance: register every output file
# ---------------------------------------------------------------------------

_TEXT_EXTS = ("csv", "tsv", "json", "txt", "fasta", "fa", "fna", "nwk", "newick", "gff3", "gff", "md", "yml", "yaml")


def _insert_pipeline_output(cur, run_id, stage, process_name, path):
    ext = path.suffix.lower().lstrip(".") if path.suffix else "none"
    row_count = None
    notes = "registered by reference"
    if ext in _TEXT_EXTS:
        notes = "text file registered"
        row_count = count_text_rows(path)
    cur.execute(
        """
        INSERT INTO pipeline_outputs (run_id, stage, process_name, file_type, file_name, file_path, file_hash, row_count, notes)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        """,
        (run_id, stage, process_name, ext, path.name, str(path), file_hash(path), row_count, notes),
    )


def register_pipeline_outputs(conn, run_id, results_dir, stage=None):
    """Register every file under a staged directory as pipeline provenance.

    `stage` should be passed explicitly whenever the caller knows exactly
    which stage this directory belongs to (e.g. a directly-staged process
    input); it's only inferred from the path via `infer_stage()` when the
    directory layout is expected to mirror the published results/ tree.
    """
    results_dir = Path(results_dir)
    if not results_dir.exists():
        print(f"  SKIP: results directory not found: {results_dir}", file=sys.stderr)
        return 0
    files = sorted(p for p in results_dir.rglob("*") if p.is_file())
    if not files:
        print(f"  No files found under {results_dir}", file=sys.stderr)
        return 0
    inserted = 0
    with conn.cursor() as cur:
        for path in files:
            rel = path.relative_to(results_dir)
            row_stage = stage if stage else infer_stage(rel)
            _insert_pipeline_output(cur, run_id, row_stage, str(rel.parent), path)
            inserted += 1
    conn.commit()
    print(f"  Registered {inserted} pipeline output files (stage={stage or 'inferred'})", file=sys.stderr)
    return inserted


def register_dir_outputs(conn, run_id, stage, base_dir):
    """Register every file under a reliably-staged process input directory
    (always an absolute, complete path as staged by Nextflow) with an
    explicit stage tag, bypassing path-keyword guessing entirely."""
    return register_pipeline_outputs(conn, run_id, base_dir, stage=stage)


def register_file_output(conn, run_id, stage, path):
    """Register a single reliably-staged process input file (e.g. auspice_json)."""
    path = Path(path)
    if not path.exists() or not path.is_file():
        print(f"  SKIP: file not found: {path}", file=sys.stderr)
        return 0
    with conn.cursor() as cur:
        _insert_pipeline_output(cur, run_id, stage, stage, path)
    conn.commit()
    print(f"  Registered 1 pipeline output file (stage={stage})", file=sys.stderr)
    return 1


# ---------------------------------------------------------------------------
# Reference genomes and genes from Nextclade dataset zip
# ---------------------------------------------------------------------------

def load_reference_genomes_and_genes(conn, run_id, results_dir):
    results_dir = Path(results_dir)
    zip_files = sorted(results_dir.rglob("nextclade-dataset.zip"))
    if not zip_files:
        print("  SKIP: no nextclade-dataset.zip found", file=sys.stderr)
        return
    loaded = 0
    with conn.cursor() as cur:
        for zf_path in zip_files:
            species = None
            pathogen = None
            for part in reversed(zf_path.parts):
                if part in ("bdbv", "sudv", "ebov", "tafv", "restv", "zaire", "sudan"):
                    species = part
                    break
            cur.execute("SAVEPOINT load_ref_sp")
            try:
                with zipfile.ZipFile(zf_path) as zf, tempfile.TemporaryDirectory() as tmp:
                    zf.extractall(tmp)
                    tmp_path = Path(tmp)
                    pathogen_json = tmp_path / "pathogen.json"
                    reference_fasta = tmp_path / "reference.fasta"
                    gff3 = tmp_path / "genome_annotation.gff3"
                    if pathogen_json.exists():
                        with open(pathogen_json) as f:
                            pdata = json.load(f)
                        attrs = pdata.get("attributes", {})
                        accession = attrs.get("reference accession") or attrs.get("reference name")
                        source = pdata.get("version", {}).get("tag", "nextclade-dataset")
                        pathogen = attrs.get("name", "orthoebolavirus")
                        if species is None:
                            species = attrs.get("name", "").split()[-1].lower()
                    else:
                        accession = None
                        source = "nextclade-dataset"
                    sequence = None
                    length = None
                    if reference_fasta.exists():
                        lines = reference_fasta.read_text().splitlines()
                        seq = "".join(l for l in lines if not l.startswith(">"))
                        sequence = seq
                        length = len(seq)
                    cur.execute(
                        """
                        INSERT INTO reference_genomes (pathogen, species, accession, source, length, sequence)
                        VALUES (%s, %s, %s, %s, %s, %s)
                        ON CONFLICT (pathogen, species, accession) DO UPDATE SET
                            source = EXCLUDED.source,
                            length = EXCLUDED.length,
                            sequence = EXCLUDED.sequence
                        RETURNING ref_id
                        """,
                        (pathogen, species, accession, source, length, sequence),
                    )
                    ref_id = cur.fetchone()[0]
                    if gff3.exists():
                        with open(gff3) as f:
                            for line in f:
                                if line.startswith("#") or not line.strip():
                                    continue
                                cols = line.rstrip("\n").split("\t")
                                if len(cols) < 9:
                                    continue
                                if cols[2] != "CDS":
                                    continue
                                start = parse_int(cols[3])
                                end = parse_int(cols[4])
                                strand = cols[6]
                                attrs = {}
                                for attr in cols[8].split(";"):
                                    if "=" in attr:
                                        k, v = attr.split("=", 1)
                                        attrs[k.strip()] = v.strip()
                                gene_name = attrs.get("Name") or attrs.get("gene") or None
                                product = attrs.get("product") or None
                                if gene_name:
                                    cur.execute(
                                        """
                                        INSERT INTO genes (ref_id, gene_name, product, start_pos, end_pos, strand)
                                        VALUES (%s, %s, %s, %s, %s, %s)
                                        ON CONFLICT DO NOTHING
                                        RETURNING gene_id
                                        """,
                                        (ref_id, gene_name, product, start, end, strand),
                                    )
                loaded += 1
            except Exception as e:
                cur.execute("ROLLBACK TO SAVEPOINT load_ref_sp")
                print(f"  WARNING: failed to process {zf_path}: {e}", file=sys.stderr)
            else:
                cur.execute("RELEASE SAVEPOINT load_ref_sp")
    conn.commit()
    print(f"  Loaded reference genomes from {loaded}/{len(zip_files)} dataset(s)", file=sys.stderr)


# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

def load_samples(conn, run_id, results_dir, species_assignments=None, metadata_tsv=None, species=None, species_screening_only=False):
    results_dir = Path(results_dir)
    if species_assignments:
        species_assignments = Path(species_assignments)
    else:
        species_assignments = results_dir / "classification" / "species_assignments.tsv"
    if not species_assignments.exists() and metadata_tsv:
        species_assignments = Path(metadata_tsv)
    if not species_assignments.exists():
        print(f"  SKIP: species assignments not found: {species_assignments}", file=sys.stderr)
        return

    rows = read_tsv_or_csv(species_assignments)
    with conn.cursor() as cur:
        for row in rows:
            sample = normalize_text(row.get("sample"))
            if not sample:
                continue
            cur.execute(
                """
                INSERT INTO samples (run_id, sample_name, pathogen, species, is_query, best_dataset_file)
                VALUES (%s, %s, %s, %s, %s, %s)
                ON CONFLICT (run_id, sample_name) DO UPDATE
                SET pathogen = EXCLUDED.pathogen,
                    species = EXCLUDED.species,
                    is_query = EXCLUDED.is_query,
                    best_dataset_file = EXCLUDED.best_dataset_file
                """,
                (
                    run_id,
                    sample,
                    normalize_text(row.get("pathogen")),
                    normalize_text(row.get("species")),
                    True,
                    normalize_text(row.get("best_dataset_file")),
                ),
            )
    conn.commit()
    print(f"  Loaded {len(rows)} samples from species assignments", file=sys.stderr)

    # Update from input/per-species metadata
    if not species_screening_only:
        if metadata_tsv:
            _update_samples_from_metadata(conn, run_id, metadata_tsv)
        else:
            # Look for metadata.tsv sibling of results or under classification
            candidates = [
                results_dir.parent / "input" / "metadata.tsv",
                results_dir / "input" / "metadata.tsv",
            ]
            for cand in candidates:
                if cand.exists():
                    _update_samples_from_metadata(conn, run_id, cand)
                    break

        _update_samples_from_metadata_extended(conn, run_id, results_dir)
    _update_samples_from_nextclade(conn, run_id, results_dir, species=species)


def _update_samples_from_metadata(conn, run_id, metadata_path):
    rows = read_tsv_or_csv(metadata_path)
    if not rows:
        return
    with conn.cursor() as cur:
        for row in rows:
            sample = normalize_text(row.get("accession") or row.get("sample") or row.get("sample_name"))
            if not sample:
                continue
            collection_date = normalize_date(row.get("date") or row.get("collection_date"))
            country = normalize_country(row.get("country"))
            admin1 = normalize_text(row.get("region"))
            admin2 = normalize_text(row.get("division"))
            locality = normalize_text(row.get("location"))
            host = normalize_text(row.get("host"))
            strain = normalize_text(row.get("strain"))
            location_id = get_or_create_location(cur, country, admin1, admin2, locality)
            cur.execute(
                """
                UPDATE samples
                SET collection_date = COALESCE(%s, collection_date),
                    country = COALESCE(%s, country),
                    admin1 = COALESCE(%s, admin1),
                    admin2 = COALESCE(%s, admin2),
                    locality = COALESCE(%s, locality),
                    host = COALESCE(%s, host),
                    strain = COALESCE(%s, strain)
                WHERE run_id = %s AND sample_name = %s
                """,
                (collection_date, country, admin1, admin2, locality, host, strain, run_id, sample),
            )
    conn.commit()
    print(f"  Updated {len(rows)} samples from metadata", file=sys.stderr)


def _update_samples_from_metadata_extended(conn, run_id, results_dir):
    files = sorted(results_dir.rglob("metadata_extended.tsv"))
    if not files:
        return
    with conn.cursor() as cur:
        for path in files:
            rows = read_tsv_or_csv(path)
            for row in rows:
                sample = normalize_text(row.get("accession"))
                if not sample:
                    continue
                cur.execute(
                    """
                    UPDATE samples
                    SET outbreak = COALESCE(%s, outbreak),
                        nextclade_qc = COALESCE(%s, nextclade_qc),
                        genome_coverage = COALESCE(%s, genome_coverage),
                        ppx_accession = COALESCE(%s, ppx_accession),
                        insdc_accession = COALESCE(%s, insdc_accession),
                        country = COALESCE(%s, country),
                        admin1 = COALESCE(%s, admin1),
                        admin2 = COALESCE(%s, admin2),
                        locality = COALESCE(%s, locality),
                        host = COALESCE(%s, host),
                        collection_date = COALESCE(%s, collection_date)
                    WHERE run_id = %s AND sample_name = %s
                    """,
                    (
                        normalize_text(row.get("outbreak")),
                        normalize_text(row.get("nextclade_qc")),
                        parse_float(row.get("genome_coverage")),
                        normalize_text(row.get("PPX_accession")),
                        normalize_text(row.get("INSDC_accession")),
                        normalize_country(row.get("country")),
                        normalize_text(row.get("region")),
                        normalize_text(row.get("division")),
                        normalize_text(row.get("location")),
                        normalize_text(row.get("host")),
                        normalize_date(row.get("date")),
                        run_id,
                        sample,
                    ),
                )
    conn.commit()
    print(f"  Updated samples from {len(files)} metadata_extended file(s)", file=sys.stderr)


def _infer_species_from_path(path):
    name = path.name.lower()
    parts = [p.lower() for p in path.parts]
    for species in ("bdbv", "sudv", "ebov", "zaire", "tafv", "restv", "sudan"):
        if species in name or species in parts:
            return species
    return None


def register_nextstrain_msa(conn, run_id, results_dir, species=None):
    """Register the nextstrain_ebola alignment.fasta + metadata_extended.tsv per
    species so pathogen_identification.R can fall back to a direct MSA-based
    identity/distance computation when a species has no usable Nextclade
    screening tree (e.g. SUDV, whose Nextclade dataset only ships a single
    reference genome with no background isolates baked into the tree).
    file_name is ambiguous across species ("alignment.fasta" everywhere), so
    species is instead encoded in process_name for unambiguous R-side lookup.
    """
    results_dir = Path(results_dir)
    species_filter = normalize_text(species)
    align_files = sorted(results_dir.rglob("alignment.fasta"))
    if species_filter:
        align_files = [p for p in align_files if species_filter in [x.lower() for x in p.parts]]
    if not align_files:
        return
    registered = 0
    with conn.cursor() as cur:
        for path in align_files:
            sp = _infer_species_from_path(path)
            if not sp:
                continue
            _insert_pipeline_output(cur, run_id, "nextstrain_ebola", f"nextstrain_msa_{sp}", path)
            registered += 1
            meta_ext = path.parent / "metadata_extended.tsv"
            if meta_ext.exists():
                _insert_pipeline_output(cur, run_id, "nextstrain_ebola", f"nextstrain_metadata_extended_{sp}", meta_ext)
    conn.commit()
    print(f"  Registered {registered} nextstrain MSA file(s)", file=sys.stderr)


def _update_samples_from_nextclade(conn, run_id, results_dir, species=None):
    results_dir = Path(results_dir)
    species_assignments = results_dir / "classification" / "species_assignments.tsv"
    if not species_assignments.exists():
        print(f"  SKIP: species assignments not found: {species_assignments}", file=sys.stderr)
        return

    assignments = read_tsv_or_csv(species_assignments)
    nextclade_dir = results_dir / "nextclade" / "results"
    if not nextclade_dir.exists():
        print(f"  SKIP: nextclade results directory not found: {nextclade_dir}", file=sys.stderr)
        return

    species_filter = normalize_text(species)
    processed = 0
    processed_species = set()
    with conn.cursor() as cur:
        cur.execute("SELECT sample_id, sample_name FROM samples WHERE run_id = %s", (run_id,))
        sample_id_by_name = {normalize_text(r[1]): r[0] for r in cur.fetchall()}
        for assignment in assignments:
            sample_name = normalize_text(assignment.get("sample"))
            row_species = normalize_text(assignment.get("species"))
            best_dataset_file = normalize_text(assignment.get("best_dataset_file"))
            if not sample_name or not row_species or not best_dataset_file:
                continue
            if species_filter and row_species != species_filter:
                continue

            tsv_path = nextclade_dir / best_dataset_file
            if not tsv_path.exists():
                print(f"  WARNING: best-match Nextclade TSV not found: {tsv_path}", file=sys.stderr)
                continue

            rows = read_tsv_or_csv(tsv_path)
            row = None
            for r in rows:
                if normalize_text(r.get("seqName")) == sample_name:
                    row = r
                    break
            if not row:
                print(f"  WARNING: sample {sample_name} not found in {tsv_path}", file=sys.stderr)
                continue

            pathogen = normalize_text(assignment.get("pathogen")) or "orthoebolavirus"
            clade = normalize_text(row.get("clade"))
            outbreak = normalize_text(row.get("outbreak"))
            nextclade_qc = normalize_text(row.get("qc.overallStatus"))
            coverage = parse_float(row.get("coverage"))
            alignment_score = parse_float(row.get("alignmentScore"))
            divergence = parse_float(row.get("divergence"))
            qc_score = parse_float(row.get("qc.overallScore"))
            nuc_substitutions = parse_int(row.get("totalSubstitutions"))
            aa_substitutions = parse_int(row.get("totalAminoacidSubstitutions"))

            # Nextclade Ebola outputs are missing 'divergence' and 'coverage';
            # derive them from the same TSV row.
            alignment_start = parse_int(row.get("alignmentStart")) or 1
            alignment_end = parse_int(row.get("alignmentEnd"))
            genome_length = parse_int(row.get("genomeLength"))
            if alignment_end and alignment_start:
                alignment_length = alignment_end - alignment_start + 1
            else:
                alignment_length = genome_length or None

            if divergence is None and nuc_substitutions is not None and alignment_length and alignment_length > 0:
                divergence = nuc_substitutions / alignment_length

            if coverage is None and alignment_length and alignment_length > 0:
                total_missing = parse_int(row.get("totalMissing")) or 0
                total_nonacgtn = parse_int(row.get("totalNonACGTNs")) or 0
                coverage = (alignment_length - total_missing - total_nonacgtn) / alignment_length

            # Update the sample with Nextclade data and raw row JSON
            cur.execute(
                """
                UPDATE samples
                SET clade = COALESCE(%s, clade),
                    outbreak = COALESCE(%s, outbreak),
                    nextclade_qc = COALESCE(%s, nextclade_qc),
                    qc_score = COALESCE(%s, qc_score),
                    genome_coverage = COALESCE(%s, genome_coverage),
                    nuc_substitution_count = COALESCE(%s, nuc_substitution_count),
                    aa_mutation_count = COALESCE(%s, aa_mutation_count),
                    alignment_score = COALESCE(%s, alignment_score),
                    divergence = COALESCE(%s, divergence),
                    nextclade_json = %s::jsonb
                WHERE run_id = %s AND sample_name = %s
                """,
                (
                    clade,
                    outbreak,
                    nextclade_qc,
                    qc_score,
                    coverage,
                    nuc_substitutions,
                    aa_substitutions,
                    alignment_score,
                    divergence,
                    json.dumps(row),
                    run_id,
                    sample_name,
                ),
            )

            # Get sample_id for linking
            cur.execute(
                "SELECT sample_id FROM samples WHERE run_id = %s AND sample_name = %s",
                (run_id, sample_name),
            )
            sample_id = cur.fetchone()
            if not sample_id:
                continue
            sample_id = sample_id[0]

            # Upsert clade and link
            if clade:
                cur.execute(
                    """
                    INSERT INTO clades (run_id, pathogen, species, clade_name)
                    VALUES (%s, %s, %s, %s)
                    ON CONFLICT (run_id, species, clade_name) DO UPDATE
                    SET pathogen = EXCLUDED.pathogen
                    RETURNING clade_id
                    """,
                    (run_id, pathogen, species, clade),
                )
                clade_id = cur.fetchone()[0]
                cur.execute(
                    """
                    INSERT INTO sample_clade (sample_id, clade_id)
                    VALUES (%s, %s)
                    ON CONFLICT DO NOTHING
                    """,
                    (sample_id, clade_id),
                )

            # Upsert outbreak and link
            if outbreak:
                cur.execute(
                    """
                    INSERT INTO outbreaks (run_id, pathogen, species, name)
                    VALUES (%s, %s, %s, %s)
                    ON CONFLICT DO NOTHING
                    RETURNING outbreak_id
                    """,
                    (run_id, pathogen, species, outbreak),
                )
                outbreak_id = cur.fetchone()
                if not outbreak_id:
                    cur.execute(
                        "SELECT outbreak_id FROM outbreaks WHERE run_id = %s AND species = %s AND name = %s",
                        (run_id, species, outbreak),
                    )
                    outbreak_id = cur.fetchone()
                if outbreak_id:
                    cur.execute(
                        """
                        INSERT INTO sample_outbreak (sample_id, outbreak_id)
                        VALUES (%s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (sample_id, outbreak_id[0]),
                    )

            # Upsert the Nextclade screening tree and all its tips once per species
            newick_path = tsv_path.with_suffix(".nwk")
            aligned_path = tsv_path.with_suffix(".aligned.fasta")

            if row_species not in processed_species:
                processed_species.add(row_species)

                # Wipe any previous Nextclade tree/tips for this run/species to avoid duplicates
                cur.execute(
                    """
                    DELETE FROM tree_tips
                    WHERE tree_id IN (
                        SELECT tree_id FROM phylogenetic_trees
                        WHERE run_id = %s AND species = %s AND tree_method = %s
                    )
                    """,
                    (run_id, row_species, "nextclade"),
                )
                cur.execute(
                    """
                    DELETE FROM phylogenetic_trees
                    WHERE run_id = %s AND species = %s AND tree_method = %s
                    """,
                    (run_id, row_species, "nextclade"),
                )

                # Prefer the full .auspice.json tree (it has all reference tips and cleaned names).
                # Use the .nwk only as a fallback because its labels are raw and may not match tree_tips.
                newick = None
                tree_source = None
                auspice_path = tsv_path.with_suffix(".auspice.json")
                auspice_tree = None
                if auspice_path.exists():
                    newick, auspice_tree = _newick_from_auspice(auspice_path)
                    tree_source = str(auspice_path)
                if not newick and newick_path.exists():
                    newick = newick_path.read_text().strip() or None
                    tree_source = str(newick_path)

                tree_id = None
                if newick:
                    cur.execute(
                        """
                        INSERT INTO phylogenetic_trees (run_id, pathogen, species, tree_method, tree_source, newick)
                        VALUES (%s, %s, %s, %s, %s, %s)
                        RETURNING tree_id
                        """,
                        (run_id, pathogen, row_species, "nextclade", tree_source, newick),
                    )
                    tree_id = cur.fetchone()[0]

                if auspice_tree and tree_id:
                    for tip in _collect_tips(auspice_tree):
                        raw_name = tip.get("name")
                        label = _clean_newick_name(raw_name)
                        if not label:
                            continue
                        node_attrs = tip.get("node_attrs", {})
                        tip_sample_id = sample_id_by_name.get(normalize_text(raw_name))
                        is_query = bool(tip_sample_id)
                        tip_clade = normalize_text(_get_node_value(node_attrs, "clade"))
                        tip_outbreak = normalize_text(_get_node_value(node_attrs, "outbreak"))
                        tip_div = parse_float(_get_node_value(node_attrs, "div"))
                        tip_cov = parse_float(_get_node_value(node_attrs, "genome_coverage"))
                        tip_qc = normalize_text(_get_node_value(node_attrs, "nextclade_qc"))
                        tip_nuc = parse_int(_get_node_value(node_attrs, "nuc_mutation_count"))
                        tip_aa = parse_int(_get_node_value(node_attrs, "aa_mutation_count"))

                        # For the query tip, use the values from the Nextclade TSV row
                        if is_query and tip_sample_id == sample_id:
                            tip_clade = clade if clade is not None else tip_clade
                            tip_outbreak = outbreak if outbreak is not None else tip_outbreak
                            tip_div = divergence if divergence is not None else tip_div
                            tip_cov = coverage if coverage is not None else tip_cov
                            tip_qc = nextclade_qc if nextclade_qc is not None else tip_qc
                            tip_nuc = nuc_substitutions if nuc_substitutions is not None else tip_nuc
                            tip_aa = aa_substitutions if aa_substitutions is not None else tip_aa

                        cur.execute(
                            """
                            INSERT INTO tree_tips (tree_id, sample_id, label, is_query, clade, outbreak, div, genome_coverage, nextclade_qc, nuc_mutation_count, aa_mutation_count)
                            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                            """,
                            (tree_id, tip_sample_id, label, is_query, tip_clade, tip_outbreak, tip_div, tip_cov, tip_qc, tip_nuc, tip_aa),
                        )
                elif tree_id:
                    # No Auspice tree available: at least record the query tip
                    cur.execute(
                        """
                        INSERT INTO tree_tips (tree_id, sample_id, label, is_query, clade, outbreak, div, genome_coverage, nextclade_qc, nuc_mutation_count, aa_mutation_count)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        """,
                        (tree_id, sample_id, _clean_newick_name(sample_name), True, clade, outbreak, divergence, coverage, nextclade_qc, nuc_substitutions, aa_substitutions),
                    )

            # Register the TSV, NWK, and aligned FASTA in pipeline_outputs
            for p in (tsv_path, newick_path, aligned_path):
                if p.exists():
                    _insert_pipeline_output(cur, run_id, "nextclade", "nextclade_screening", p)

            processed += 1

    conn.commit()
    print(f"  Updated {processed} samples from best-match Nextclade outputs", file=sys.stderr)


# ---------------------------------------------------------------------------
# Trees and tips
# ---------------------------------------------------------------------------

def _get_node_value(node_attrs, key):
    val = node_attrs.get(key)
    if isinstance(val, dict):
        return val.get("value")
    return val


def _clean_newick_name(name):
    name = (name or "").strip()
    # Newick metacharacters and spaces (unquoted names cannot contain these)
    name = re.sub(r"[,;\(\):\[\] ]", "_", name)
    if not name:
        name = "unnamed"
    return name


def _newick_from_auspice_node(node, parent_div=0.0):
    attrs = node.get("node_attrs", {}) if isinstance(node, dict) else {}
    div = parse_float(_get_node_value(attrs, "div")) or 0.0
    length = max(div - parent_div, 0.0)
    name = _clean_newick_name(node.get("name")) if isinstance(node, dict) else "unnamed"
    children = node.get("children") if isinstance(node, dict) else None
    if children:
        subtrees = [_newick_from_auspice_node(child, parent_div=div) for child in children]
        return f'({",".join(subtrees)}){name}:{length}'
    return f'{name}:{length}'


def _newick_from_auspice(path):
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    tree = data.get("tree") if isinstance(data, dict) else None
    if not tree:
        return None, None
    return f'({_newick_from_auspice_node(tree)});', tree


def _collect_tips(tree, tips=None):
    if tips is None:
        tips = []
    children = tree.get("children")
    if not children:
        tips.append(tree)
    else:
        for child in children:
            _collect_tips(child, tips)
    return tips


def load_trees_and_tips(conn, run_id, results_dir, auspice_json=None, iqtree=None, species=None):
    results_dir = Path(results_dir)
    species_filter = normalize_text(species)
    # Auspice JSON files under nextstrain_ebola/<species>/auspice/
    auspice_files = []
    if auspice_json and Path(auspice_json).exists():
        auspice_files = [Path(auspice_json)]
    else:
        auspice_files = sorted(results_dir.rglob("*all-outbreaks.json"))
    if species_filter:
        auspice_files = [p for p in auspice_files if species_filter in [x.lower() for x in p.parts]]
    if not auspice_files:
        print("  SKIP: no Auspice JSON files found", file=sys.stderr)
        return
    loaded = 0
    with conn.cursor() as cur:
        for path in auspice_files:
            species = None
            for part in path.parts:
                if part.lower() in ("bdbv", "sudv", "ebov", "zaire", "tafv", "restv", "sudan"):
                    species = part.lower()
                    break
            newick_candidates = sorted(results_dir.rglob("tree.nwk"))
            if species_filter:
                newick_candidates = [p for p in newick_candidates if species_filter in [x.lower() for x in p.parts]]
            newick = newick_candidates[0].read_text().strip() if newick_candidates else None
            tree_source = str(newick_candidates[0]) if newick_candidates else str(path)
            cur.execute(
                """
                INSERT INTO phylogenetic_trees (run_id, pathogen, species, tree_method, tree_source, newick)
                VALUES (%s, %s, %s, %s, %s, %s)
                RETURNING tree_id
                """,
                (run_id, "orthoebolavirus", species, "nextstrain", tree_source, newick),
            )
            tree_id = cur.fetchone()[0]
            # Load tip metadata if present
            tip_candidates = sorted(results_dir.rglob("*_tip_metadata.tsv"))
            if species_filter:
                tip_candidates = [p for p in tip_candidates if species_filter in [x.lower() for x in p.parts]]
            tip_path = tip_candidates[0] if tip_candidates else None
            if tip_path and tip_path.exists():
                for row in read_tsv_or_csv(tip_path):
                    label = normalize_text(row.get("label"))
                    if not label:
                        continue
                    is_query = str(row.get("is_query", "")).strip().lower() in ("true", "1", "yes")
                    sample_id = None
                    if is_query:
                        cur.execute(
                            "SELECT sample_id FROM samples WHERE run_id = %s AND sample_name = %s",
                            (run_id, label),
                        )
                        r = cur.fetchone()
                        if r:
                            sample_id = r[0]
                    country = normalize_country(row.get("country"))
                    admin1 = normalize_text(row.get("division"))
                    locality = normalize_text(row.get("location"))
                    location_id = get_or_create_location(cur, country, admin1, None, locality)
                    cur.execute(
                        """
                        INSERT INTO tree_tips (
                            tree_id, sample_id, label, is_query, ppx_accession, insdc_accession,
                            country, admin1, admin2, locality, host, outbreak, tip_date, div,
                            genome_coverage, nextclade_qc, aa_mutation_count, nuc_mutation_count
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (
                            tree_id, sample_id, label, is_query,
                            normalize_text(row.get("PPX_accession")),
                            normalize_text(row.get("INSDC_accession")),
                            country,
                            admin1,
                            None,
                            locality,
                            normalize_text(row.get("host")),
                            normalize_text(row.get("outbreak")),
                            normalize_date(row.get("date")),
                            parse_float(row.get("div")),
                            parse_float(row.get("genome_coverage")),
                            normalize_text(row.get("nextclade_qc")),
                            parse_int(row.get("aa_mutation_count")),
                            parse_int(row.get("nuc_mutation_count")),
                        ),
                    )
                    loaded += 1
            else:
                # Fall back to parsing the Auspice JSON for tip attributes
                with open(path) as f:
                    data = json.load(f)
                for tip in _collect_tips(data.get("tree", {})):
                    label = normalize_text(tip.get("name"))
                    if not label:
                        continue
                    node_attrs = tip.get("node_attrs", {})
                    cur.execute(
                        "SELECT sample_id FROM samples WHERE run_id = %s AND sample_name = %s",
                        (run_id, label),
                    )
                    r = cur.fetchone()
                    sample_id = r[0] if r else None
                    is_query = sample_id is not None
                    country = normalize_country(_get_node_value(node_attrs, "country"))
                    admin1 = normalize_text(_get_node_value(node_attrs, "division"))
                    locality = normalize_text(_get_node_value(node_attrs, "location"))
                    location_id = get_or_create_location(cur, country, admin1, None, locality)
                    cur.execute(
                        """
                        INSERT INTO tree_tips (
                            tree_id, sample_id, label, is_query, ppx_accession, insdc_accession,
                            country, admin1, admin2, locality, host, outbreak, tip_date, div,
                            genome_coverage, nextclade_qc, aa_mutation_count, nuc_mutation_count
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (
                            tree_id, sample_id, label, is_query,
                            normalize_text(_get_node_value(node_attrs, "PPX_accession")),
                            normalize_text(_get_node_value(node_attrs, "INSDC_accession")),
                            country,
                            admin1,
                            None,
                            locality,
                            normalize_text(_get_node_value(node_attrs, "host")),
                            normalize_text(_get_node_value(node_attrs, "outbreak")),
                            normalize_date(_get_node_value(node_attrs, "date")),
                            parse_float(_get_node_value(node_attrs, "div")),
                            parse_float(_get_node_value(node_attrs, "genome_coverage")),
                            normalize_text(_get_node_value(node_attrs, "nextclade_qc")),
                            parse_int(_get_node_value(node_attrs, "aa_mutation_count")),
                            parse_int(_get_node_value(node_attrs, "nuc_mutation_count")),
                        ),
                    )
                    loaded += 1
        if iqtree:
            iqtree_path = Path(iqtree)
            if iqtree_path.exists():
                species = _infer_species_from_path(iqtree_path)
                newick = iqtree_path.read_text().strip()
                cur.execute(
                    """
                    INSERT INTO phylogenetic_trees (run_id, pathogen, species, tree_method, tree_source, newick)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    """,
                    (run_id, "orthoebolavirus", species, "iqtree2", str(iqtree_path), newick),
                )
        conn.commit()
    print(f"  Loaded {len(auspice_files)} tree(s) and {loaded} tip(s)", file=sys.stderr)


# ---------------------------------------------------------------------------
# Mutations
# ---------------------------------------------------------------------------

def load_mutations(conn, run_id, results_dir):
    results_dir = Path(results_dir)
    with conn.cursor() as cur:
        # From query_mutations.tsv
        for path in sorted(results_dir.rglob("*_query_mutations.tsv")):
            for row in read_tsv_or_csv(path):
                sample = normalize_text(row.get("sample"))
                gene = normalize_text(row.get("gene"))
                pos = parse_int(row.get("position"))
                ref_aa = normalize_text(row.get("ref_aa"))
                alt_aa = normalize_text(row.get("alt_aa"))
                label = normalize_text(row.get("mutation_label"))
                if not sample or not gene or not pos:
                    continue
                protein_id = _protein_id_for_gene(cur, gene)
                mutation_id = _insert_mutation(cur, protein_id, gene, label, ref_aa, pos, alt_aa)
                cur.execute(
                    "SELECT sample_id FROM samples WHERE run_id = %s AND sample_name = %s",
                    (run_id, sample),
                )
                r = cur.fetchone()
                if r:
                    cur.execute(
                        """
                        INSERT INTO sample_mutation (sample_id, mutation_id) VALUES (%s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (r[0], mutation_id),
                    )

        # From Nextclade TSV aaSubstitutions (filter to matching species from filename)
        for path in sorted(results_dir.rglob("nextclade.tsv")) + sorted(results_dir.rglob("*_nextclade.tsv")):
            file_species = _infer_species_from_path(path)
            for row in read_tsv_or_csv(path):
                sample = normalize_text(row.get("seqName"))
                aa_subs = str(row.get("aaSubstitutions", "")).strip()
                if not sample or not aa_subs or aa_subs in ("NA", ""):
                    continue
                where = "run_id = %s AND sample_name = %s"
                params = [run_id, sample]
                if file_species:
                    where += " AND (species = %s OR species IS NULL)"
                    params.append(file_species)
                cur.execute(
                    f"SELECT sample_id FROM samples WHERE {where}",
                    params,
                )
                r = cur.fetchone()
                if not r:
                    continue
                sample_id = r[0]
                for mut in aa_subs.split(","):
                    mut = mut.strip()
                    if not mut:
                        continue
                    parsed = _parse_aa_mutation(mut)
                    if not parsed:
                        continue
                    gene, ref_aa, pos, alt_aa = parsed
                    protein_id = _protein_id_for_gene(cur, gene)
                    mutation_id = _insert_mutation(cur, protein_id, gene, mut, ref_aa, pos, alt_aa)
                    cur.execute(
                        """
                        INSERT INTO sample_mutation (sample_id, mutation_id) VALUES (%s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (sample_id, mutation_id),
                    )

        # From mutation_matrix.tsv: if a row label matches a sample name and value == 1
        for path in sorted(results_dir.rglob("*_mutation_matrix.tsv")):
            rows = read_tsv_or_csv(path)
            if not rows:
                continue
            label_key = "label" if "label" in rows[0] else None
            if not label_key:
                continue
            for row in rows:
                label = normalize_text(row.get(label_key))
                if not label:
                    continue
                cur.execute(
                    "SELECT sample_id FROM samples WHERE run_id = %s AND sample_name = %s",
                    (run_id, label),
                )
                r = cur.fetchone()
                if not r:
                    continue
                sample_id = r[0]
                for col, val in row.items():
                    if col == label_key:
                        continue
                    if str(val).strip() not in ("1", "True", "true"):
                        continue
                    parsed = _parse_mutation_column(col)
                    if not parsed:
                        continue
                    gene, ref_aa, pos, alt_aa = parsed
                    protein_id = _protein_id_for_gene(cur, gene)
                    mutation_id = _insert_mutation(cur, protein_id, gene, col, ref_aa, pos, alt_aa)
                    cur.execute(
                        """
                        INSERT INTO sample_mutation (sample_id, mutation_id) VALUES (%s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (sample_id, mutation_id),
                    )
    conn.commit()
    print("  Loaded mutations", file=sys.stderr)


def _parse_aa_mutation(mut):
    # e.g. GP:Y151F or L:S85-
    m = re.fullmatch(r"([A-Za-z0-9_]+):([A-Za-z\*])([0-9]+)([A-Za-z\*\-])", mut)
    if not m:
        return None
    return m.group(1), m.group(2), parse_int(m.group(3)), m.group(4)


def _parse_mutation_column(col):
    # e.g. L:K135M or GP:A489V
    m = re.fullmatch(r"([A-Za-z0-9_]+):([A-Za-z\*])([0-9]+)([A-Za-z\*\-])", col)
    if not m:
        return None
    return m.group(1), m.group(2), parse_int(m.group(3)), m.group(4)


def _protein_id_for_gene(cur, gene):
    if not gene:
        return None
    cur.execute(
        "SELECT protein_id FROM proteins WHERE uniprot_accession IS NULL AND protein_name = %s LIMIT 1",
        (gene,),
    )
    r = cur.fetchone()
    if r:
        return r[0]
    # Try gene name match
    cur.execute(
        "SELECT p.protein_id FROM proteins p JOIN genes g ON p.gene_id = g.gene_id WHERE g.gene_name = %s LIMIT 1",
        (gene,),
    )
    r = cur.fetchone()
    if r:
        return r[0]
    return None


def _insert_mutation(cur, protein_id, gene, label, ref_aa, pos, alt_aa):
    if not label:
        return None
    mutation_type = "substitution"
    if alt_aa == "-":
        mutation_type = "deletion"
    if ref_aa == "-":
        mutation_type = "insertion"
    cur.execute(
        """
        INSERT INTO mutations (protein_id, gene_id, mutation_label, ref_aa, position, alt_aa, mutation_type)
        VALUES (%s, (SELECT gene_id FROM genes WHERE gene_name = %s LIMIT 1), %s, %s, %s, %s, %s)
        ON CONFLICT (protein_id, mutation_label) DO UPDATE SET
            ref_aa = EXCLUDED.ref_aa,
            position = EXCLUDED.position,
            alt_aa = EXCLUDED.alt_aa,
            mutation_type = EXCLUDED.mutation_type
        RETURNING mutation_id
        """,
        (protein_id, gene, label, ref_aa, pos, alt_aa, mutation_type),
    )
    return cur.fetchone()[0]


# ---------------------------------------------------------------------------
# Phenotype annotations (proteins, discovery, mutagenesis)
# ---------------------------------------------------------------------------

def _rglob_dirs(dirs, patterns):
    hits = []
    for d in dirs:
        for p in patterns:
            hits.extend(d.rglob(p))
    return sorted(set(hits))


def _first_col(row, *names):
    """Return the first non-empty value among several known column-name
    variants. UniProt annotation outputs go through R (UniprotR/UniProtExtractR),
    which sanitizes column names differently depending on the export step
    (e.g. "Function [CC]" -> "Function..CC." -> "func.Function..CC."), so the
    same semantic field can show up under several different exact header
    strings depending on which output file is being read.
    """
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return value
    return None


def load_phenotype_annotations(conn, run_id, uniprotr_dir=None, extractr_dir=None, rbioapi_dir=None):
    search_dirs = []
    if uniprotr_dir:
        search_dirs.append(Path(uniprotr_dir))
    if extractr_dir:
        search_dirs.append(Path(extractr_dir))
    if rbioapi_dir:
        search_dirs.append(Path(rbioapi_dir))
    if not search_dirs:
        print("  SKIP: no phenotype annotation directories provided", file=sys.stderr)
        return
    with conn.cursor() as cur:
        # Discovery.tsv -> proteins and entity_evidence
        for path in _rglob_dirs(search_dirs, ["*_discovery.tsv"]):
            for row in read_tsv_or_csv(path):
                sample = normalize_text(row.get("sample"))
                gene = normalize_text(row.get("gene"))
                reason = normalize_text(row.get("reason"))
                accessions = str(row.get("uniprot_accessions", "")).strip()
                if not sample or not gene:
                    continue
                if accessions:
                    for acc in re.split(r"[,|]", accessions):
                        acc = acc.strip()
                        if not acc:
                            continue
                        protein_id = _get_or_create_protein(cur, gene, acc)
                        cur.execute(
                            """
                            SELECT sample_id FROM samples WHERE run_id = %s AND sample_name = %s
                            """,
                            (run_id, sample),
                        )
                        r = cur.fetchone()
                        if r:
                            cur.execute(
                                """
                                INSERT INTO entity_evidence (table_name, entity_id, source_id, confidence, reference)
                                VALUES ('proteins', %s, NULL, %s, %s)
                                ON CONFLICT DO NOTHING
                                """,
                                (protein_id, reason, f"discovery:{path.name}"),
                            )

        # UniProt download TSV/CSV -> enrich proteins. Different annotation
        # stages sanitize UniProt's column names differently (raw download vs
        # UniProtExtractR vs the combined UniprotR export), so every field is
        # looked up across all observed variants via _first_col().
        for path in _rglob_dirs(search_dirs, [
            "*_uniprot_download.tsv",
            "*_uniprotextractr_clean.tsv",
            "*_uniprotextractr.tsv",
            "*_uniprotextractr.csv",
            "*_uniprotr_combined.tsv",
            "*_uniprotr_combined.csv",
        ]):
            rows = read_tsv_or_csv(path)
            if not rows:
                continue
            for row in rows:
                acc = normalize_text(_first_col(row, "Entry", "taxa.Entry", "accession"))
                if not acc:
                    continue
                gene = normalize_text(_first_col(
                    row, "Gene.Names", "Gene Names", "Gene Names (primary )", "taxa.Gene.Names",
                ))
                protein_name = normalize_text(_first_col(row, "Protein.names", "Protein names", "taxa.Protein.names"))
                organism = normalize_text(_first_col(row, "Organism", "taxa.Organism"))
                length = parse_int(_first_col(row, "Length"))
                protein_id = _get_or_create_protein(cur, gene, acc, protein_name=protein_name, organism=organism, length=length)
                # Functions
                func_text = _first_col(row, "Function [CC]", "Function..CC.", "func.Function..CC.")
                if func_text:
                    cur.execute(
                        """
                        INSERT INTO protein_functions (protein_id, source, function_text)
                        VALUES (%s, %s, %s) ON CONFLICT DO NOTHING
                        """,
                        (protein_id, "UniProtKB", func_text),
                    )
                # GO terms
                go_ids = _first_col(
                    row, "Gene Ontology (GO)", "Gene Ontology IDs", "Gene.Ontology.IDs", "go.Gene.Ontology.IDs",
                )
                if go_ids:
                    for go_id in re.findall(r"GO:\d+", go_ids):
                        cur.execute(
                            """
                            INSERT INTO go_terms (protein_id, go_id, source)
                            VALUES (%s, %s, %s) ON CONFLICT DO NOTHING
                            """,
                            (protein_id, go_id, "UniProtKB"),
                        )
                # Domains
                domain_text = _first_col(
                    row, "Domain [FT]", "Domain..FT.", "Domain..FT.edit", "family.Domain..FT.",
                )
                if domain_text:
                    for m in re.finditer(r"DOMAIN\s+(\d+)\.\.(\d+);\s*/note=\"([^\"]+)\"", domain_text):
                        cur.execute(
                            """
                            INSERT INTO protein_domains (protein_id, domain_name, start_pos, end_pos, source)
                            VALUES (%s, %s, %s, %s, %s) ON CONFLICT DO NOTHING
                            """,
                            (protein_id, m.group(3), parse_int(m.group(1)), parse_int(m.group(2)), "UniProtKB"),
                        )

        # rbioapi mutagenesis/variation -> mutation_phenotypes
        for path in _rglob_dirs(search_dirs, ["*_mutagenesis.tsv", "*_variation.tsv"]):
            rows = read_tsv_or_csv(path)
            if not rows:
                continue
            for row in rows:
                gene = normalize_text(row.get("gene"))
                pos_start = parse_int(row.get("position_start"))
                pos_end = parse_int(row.get("position_end"))
                ref_aa = normalize_text(row.get("original_aa"))
                alt_aa = normalize_text(row.get("alternative_aa"))
                description = row.get("description") or row.get("consequence") or row.get("notes")
                mutation_id = None
                if gene and pos_start:
                    cur.execute(
                        """
                        SELECT mutation_id FROM mutations
                        WHERE position = %s AND ref_aa = %s AND alt_aa = %s
                          AND (protein_id IN (SELECT protein_id FROM proteins WHERE protein_name = %s)
                               OR gene_id = (SELECT gene_id FROM genes WHERE gene_name = %s LIMIT 1))
                        LIMIT 1
                        """,
                        (pos_start, ref_aa, alt_aa, gene, gene),
                    )
                    r = cur.fetchone()
                    if r:
                        mutation_id = r[0]
                cur.execute(
                    """
                    INSERT INTO mutation_phenotypes (mutation_id, phenotype, effect, evidence, source)
                    VALUES (%s, %s, %s, %s, %s)
                    """,
                    (mutation_id, "mutagenesis" if "mutagenesis" in path.name else "variation", description, str(path.name), "rbioapi"),
                )
    conn.commit()
    print("  Loaded phenotype annotations", file=sys.stderr)


# ---------------------------------------------------------------------------
# HMMER/Pfam domain annotations
# ---------------------------------------------------------------------------

def _split_query_id(query_id):
    if "_" in query_id:
        return query_id.split("_", 1)
    return query_id, ""


def _ensure_hmm_domain_columns(cur):
    """Back-fill HMMER columns on warehouses created before this schema change."""
    cur.execute(
        """
        ALTER TABLE protein_domains
        ADD COLUMN IF NOT EXISTS hmm_accession TEXT,
        ADD COLUMN IF NOT EXISTS evalue DOUBLE PRECISION,
        ADD COLUMN IF NOT EXISTS bit_score DOUBLE PRECISION,
        ADD COLUMN IF NOT EXISTS hmm_start INTEGER,
        ADD COLUMN IF NOT EXISTS hmm_end INTEGER,
        ADD COLUMN IF NOT EXISTS ali_start INTEGER,
        ADD COLUMN IF NOT EXISTS ali_end INTEGER,
        ADD COLUMN IF NOT EXISTS env_start INTEGER,
        ADD COLUMN IF NOT EXISTS env_end INTEGER;

        ALTER TABLE protein_domains
        ALTER COLUMN evalue TYPE DOUBLE PRECISION,
        ALTER COLUMN bit_score TYPE DOUBLE PRECISION;
        """
    )


def load_hmm_annotations(conn, run_id, hmm_dir, species=None):
    hmm_dir = Path(hmm_dir) if hmm_dir else None
    if not hmm_dir or not hmm_dir.exists():
        print("  SKIP: no HMM annotations directory provided", file=sys.stderr)
        return
    paths = sorted(hmm_dir.rglob("*_hmm_domain_table.txt"))
    if not paths:
        print(f"  SKIP: no HMM domain tables found in {hmm_dir}", file=sys.stderr)
        return
    loaded = 0
    with conn.cursor() as cur:
        _ensure_hmm_domain_columns(cur)
        for path in paths:
            for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                # hmmscan --domtblout: 22 fixed fields then the description
                parts = line.split(None, 22)
                if len(parts) < 22:
                    continue
                target_name = normalize_text(parts[0])
                target_accession = normalize_text(parts[1])
                query_id = normalize_text(parts[3])
                evalue = parse_float(parts[12])  # i-Evalue (domain independent)
                bit_score = parse_float(parts[13])
                hmm_from = parse_int(parts[15])
                hmm_to = parse_int(parts[16])
                ali_from = parse_int(parts[17])
                ali_to = parse_int(parts[18])
                env_from = parse_int(parts[19])
                env_to = parse_int(parts[20])
                description = normalize_text(parts[22]) if len(parts) > 22 else ""
                if not query_id or not target_name:
                    continue
                sample, gene = _split_query_id(query_id)
                sample_id = None
                if sample:
                    cur.execute(
                        "SELECT sample_id FROM samples WHERE run_id = %s AND sample_name = %s",
                        (run_id, sample),
                    )
                    r = cur.fetchone()
                    if r:
                        sample_id = r[0]
                protein_id = _get_or_create_protein(
                    cur, gene, None,
                    protein_name=query_id, organism=species, length=parse_int(parts[5]),
                )
                cur.execute(
                    """
                    INSERT INTO protein_domains (
                        protein_id, domain_name, start_pos, end_pos, source,
                        hmm_accession, evalue, bit_score,
                        hmm_start, hmm_end, ali_start, ali_end, env_start, env_end
                    )
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT DO NOTHING
                    """,
                    (
                        protein_id, target_name, env_from, env_to, "HMMER/Pfam",
                        target_accession, evalue, bit_score,
                        hmm_from, hmm_to, ali_from, ali_to, env_from, env_to,
                    ),
                )
                if description:
                    cur.execute(
                        """
                        INSERT INTO protein_functions (protein_id, source, function_text)
                        VALUES (%s, %s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (protein_id, "HMMER/Pfam", description),
                    )
                if sample_id:
                    cur.execute(
                        """
                        INSERT INTO entity_evidence (table_name, entity_id, source_id, confidence, reference)
                        VALUES ('proteins', %s, NULL, %s, %s)
                        ON CONFLICT DO NOTHING
                        """,
                        (protein_id, "HMMER/Pfam", f"hmm:{path.name}"),
                    )
                loaded += 1
    conn.commit()
    print(f"  Loaded {loaded} HMM domain rows", file=sys.stderr)


def _get_or_create_protein(cur, gene, uniprot_accession, protein_name=None, organism=None, length=None):
    if uniprot_accession:
        cur.execute(
            "SELECT protein_id FROM proteins WHERE uniprot_accession = %s LIMIT 1",
            (uniprot_accession,),
        )
        r = cur.fetchone()
        if r:
            return r[0]
    # If no accession, look by gene placeholder
    name = protein_name or gene
    cur.execute(
        "SELECT protein_id FROM proteins WHERE uniprot_accession IS NULL AND protein_name = %s LIMIT 1",
        (name,),
    )
    r = cur.fetchone()
    if r:
        if uniprot_accession:
            cur.execute(
                "UPDATE proteins SET uniprot_accession = %s WHERE protein_id = %s",
                (uniprot_accession, r[0]),
            )
        return r[0]
    # Need a gene_id
    cur.execute(
        "SELECT gene_id FROM genes WHERE gene_name = %s LIMIT 1",
        (gene,),
    )
    r = cur.fetchone()
    gene_id = r[0] if r else None
    cur.execute(
        """
        INSERT INTO proteins (gene_id, uniprot_accession, protein_name, organism, length)
        VALUES (%s, %s, %s, %s, %s)
        RETURNING protein_id
        """,
        (gene_id, uniprot_accession, name, organism, length),
    )
    return cur.fetchone()[0]


# ---------------------------------------------------------------------------
# Epidemiological data
# ---------------------------------------------------------------------------

def load_epi_data(conn, run_id, epi_raw_dir, epi_search_summary=None, species=None):
    if not epi_raw_dir:
        print("  SKIP: no epi raw directory provided", file=sys.stderr)
        return
    epi_dir = Path(epi_raw_dir)
    if not epi_dir.exists():
        print(f"  SKIP: epi data directory not found: {epi_dir}", file=sys.stderr)
        return
    if species:
        species_dir = epi_dir / species
        if not species_dir.exists():
            print(f"  SKIP: no epi data subfolder for species '{species}': {species_dir}", file=sys.stderr)
            return
        epi_dir = species_dir
    csv_files = sorted(epi_dir.glob("*.csv"))
    if not csv_files:
        print(f"  SKIP: no CSV files in {epi_dir}", file=sys.stderr)
        return

    with conn.cursor() as cur:
        for csv_path in csv_files:
            print(f"  Loading epi dataset: {csv_path.name}", file=sys.stderr)
            rows = read_tsv_or_csv(csv_path)
            if not rows:
                continue

            # Detect dataset type from filename/columns
            dataset_type = _infer_dataset_type(csv_path.name, rows[0].keys())
            pathogen = None
            if species:
                pathogen = "orthoebolavirus"
            cur.execute(
                """
                INSERT INTO epidemiological_datasets (run_id, dataset_name, source, row_count, pathogen, species, dataset_type)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
                RETURNING dataset_id
                """,
                (run_id, csv_path.stem, "HDX", len(rows), pathogen, species, dataset_type),
            )
            dataset_id = cur.fetchone()[0]

            for row in rows:
                _load_epi_row(cur, dataset_id, row, dataset_type)
    conn.commit()
    print(f"  Loaded {len(csv_files)} epidemiological datasets", file=sys.stderr)


def _infer_dataset_type(name, columns):
    name_l = name.lower()
    if "subnational" in name_l:
        return "subnational"
    if "summary" in name_l or "outbreaks" in name_l or "before" in name_l:
        return "summary"
    if "cases-and-deaths" in name_l or "cumulative" in name_l or "national" in name_l:
        return "wide"
    if "indicator" in columns or "measure" in columns or "value" in columns:
        return "long"
    return "wide"


def _load_epi_row(cur, dataset_id, row, dataset_type):
    raw = {k: v for k, v in row.items() if v is not None and str(v).strip()}

    # Date extraction
    record_date = None
    report_date = None
    reference_date = None
    for k in ("Date", "date", "DATE", "report_date", "publication_date", "reference_date"):
        if k in row and row[k]:
            d = normalize_date(row[k])
            if d:
                if k.lower() in ("report_date", "publication_date"):
                    report_date = d
                elif k.lower() == "reference_date":
                    reference_date = d
                else:
                    record_date = d

    # Location
    country = normalize_country(row.get("Country") or row.get("country") or row.get("COUNTRY") or row.get("location_country"))
    admin1 = normalize_text(row.get("region") or row.get("division") or row.get("REGION") or row.get("location_name"))
    admin2 = normalize_text(row.get("admin2") or row.get("location_name"))
    locality = normalize_text(row.get("location") or row.get("locality"))
    location_code = normalize_text(row.get("location_code") or row.get("ID_REGION") or row.get("ISO_3"))
    location_code_type = normalize_text(row.get("location_code_type") or ("pcode" if "pcode" in str(row.get("location_code_type", "")) else None))
    location_id = get_or_create_location(cur, country, admin1, admin2, locality, location_code, location_code_type)

    # Default fields
    measure = None
    case_classification = None
    time_period = None
    value = None
    unit = "count"
    cases = None
    deaths = None
    suspected = None
    recovered = None
    gender_breakdown = None
    source_url = normalize_text(row.get("source_url"))
    source_indicator_name = normalize_text(row.get("source_indicator_name"))
    indicator_label = normalize_text(row.get("indicator_label") or row.get("Indicator"))

    if dataset_type == "long" or "measure" in row:
        measure = normalize_text(row.get("measure"))
        case_classification = normalize_text(row.get("case_classification"))
        time_period = normalize_text(row.get("time_period"))
        value = parse_float(row.get("value"))
        if not measure and indicator_label:
            measure, case_classification, time_period = _parse_indicator(indicator_label)
    elif dataset_type == "wide":
        # Try to extract known numeric columns
        cases = _coalesce_int(row, ["total_cases", "EBOLA_CONFIRMED_CASE", "CONFIRMED_EBOLA_CASE", "confirmed_cases", "cases"])
        deaths = _coalesce_int(row, ["total_deaths", "EBOLA_DEATHS", "CONFIRMED_EBOLA_DEATHS", "confirmed_deaths", "deaths"])
        suspected = _coalesce_int(row, ["total_suspected_cases", "EBOLA_SUSPECTED_CASE", "SUSPECTED_EBOLA_CASES", "suspected_cases", "suspected"])
        recovered = _coalesce_int(row, ["total_cured", "RECOVERIES", "recovered", "cured"])
        # Demographics as JSONB
        gender_breakdown = {}
        for k in ["WOMEN", "MEN", "GENDER_NOT_SPECIFIED", "DEATHS_WOMEN", "DEATHS_MEN", "DEATHS_GENDER_NOT_SPECIFIED", "RECOVERIES_WOMEN", "RECOVERIES_MEN"]:
            if k in row and row[k]:
                v = parse_int(row[k])
                if v is not None:
                    gender_breakdown[k] = v
        if not gender_breakdown:
            gender_breakdown = None
    elif dataset_type == "summary":
        cases = parse_int(row.get("Reported number of human cases"))
        deaths = parse_int(row.get("Reported number of deaths among cases"))
        measure = "cases"
        case_classification = "reported"
        subtype = normalize_text(row.get("Ebola subtype"))
        if subtype:
            indicator_label = subtype

    cur.execute(
        """
        INSERT INTO epidemiological_records (
            dataset_id, record_date, report_date, reference_date, location_id,
            country, admin1, admin2, locality, location_code, location_code_type,
            measure, case_classification, time_period, value, unit,
            cases, deaths, suspected, recovered, gender_breakdown,
            source_url, source_indicator_name, indicator_label, raw_data
        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """,
        (
            dataset_id, record_date, report_date, reference_date, location_id,
            country, admin1, admin2, locality, location_code, location_code_type,
            measure, case_classification, time_period, value, unit,
            cases, deaths, suspected, recovered, json.dumps(gender_breakdown) if gender_breakdown else None,
            source_url, source_indicator_name, indicator_label, json.dumps(raw),
        ),
    )


def _coalesce_int(row, keys):
    for k in keys:
        if k in row and row[k]:
            v = parse_int(row[k])
            if v is not None:
                return v
    return None


def _parse_indicator(text):
    text_l = text.lower()
    measure = None
    classification = None
    period = None
    if "case" in text_l:
        measure = "cases"
    elif "death" in text_l:
        measure = "deaths"
    elif "recover" in text_l:
        measure = "recovered"
    if "confirm" in text_l:
        classification = "confirmed"
    elif "probable" in text_l:
        classification = "probable"
    elif "suspect" in text_l:
        classification = "suspected"
    if "cumulative" in text_l:
        period = "cumulative"
    elif "new" in text_l:
        period = "new"
    return measure, classification, period


# ---------------------------------------------------------------------------
# MultiQC custom-content summary
# ---------------------------------------------------------------------------

# Key tables to report row counts for in the MultiQC summary. Kept generic
# (not pathogen-specific) so it works for any future pathogen workflow.
SUMMARY_TABLES = [
    "samples",
    "reference_genomes",
    "genomes",
    "mutations",
    "proteins",
    "phylogenetic_trees",
    "tree_tips",
    "epidemiological_records",
    "pipeline_outputs",
    "literature_domains",
    "literature_papers",
    "literature_extractions",
]


def write_mqc_summary(conn, meta_id, outdir, prefix):
    """Write a MultiQC custom-content TSV reporting row counts per table.

    This also gives BUILD_KNOWLEDGE_DB a real downstream consumer (MultiQC),
    so its outputs aren't a pipeline dead-end.
    """
    counts = {}
    with conn.cursor() as cur:
        for table in SUMMARY_TABLES:
            cur.execute(f"SELECT count(*) FROM {table}")
            counts[table] = cur.fetchone()[0]

    summary_path = Path(outdir) / f"{prefix}_knowledge_warehouse_mqc.tsv"
    with open(summary_path, "w") as f:
        f.write("# id: 'knowledge_warehouse_summary'\n")
        f.write("# section_name: 'Knowledge Warehouse'\n")
        f.write("# description: 'Row counts loaded into the per-run PostgreSQL knowledge warehouse.'\n")
        f.write("# plot_type: 'table'\n")
        f.write("# pconfig:\n")
        f.write("#     id: 'knowledge_warehouse_summary_table'\n")
        f.write("#     namespace: 'Knowledge Warehouse'\n")
        f.write("Sample\t" + "\t".join(SUMMARY_TABLES) + "\n")
        f.write(meta_id + "\t" + "\t".join(str(counts[t]) for t in SUMMARY_TABLES) + "\n")

    print(f"  MultiQC summary written to {summary_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Database dump
# ---------------------------------------------------------------------------

def dump_database(pg, outdir, prefix):
    dump_path = Path(outdir) / f"{prefix}_genomic_intelligence.sql"
    pg_dump = find_executable("pg_dump")
    env = {"PGUSER": pg.superuser}
    with open(dump_path, "w") as f:
        run_cmd([
            pg_dump,
            "-h", pg.host,
            "-p", str(pg.port),
            "-U", pg.superuser,
            "-d", pg.db_name,
            "--clean",
            "--if-exists",
        ], check=True, env=env, stdout=f)
    print(f"  Database dump written to {dump_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, help="Output directory for warehouse files")
    parser.add_argument("--meta-id", required=True, help="Run / group identifier")
    parser.add_argument("--prefix", default="query", help="File prefix for output artifacts")
    parser.add_argument("--results-dir", help="Root results directory (used to discover classification, nextclade, etc.)")
    parser.add_argument("--bioinfo-dir", help="Staged bioinformatics results directory (used to register Nextstrain outputs in pipeline_outputs)")
    parser.add_argument("--species-assignments", help="species_assignments.tsv")
    parser.add_argument("--metadata-tsv", help="Per-species or input metadata TSV")
    parser.add_argument("--epi-raw-dir", help="Directory containing epidemiological CSVs")
    parser.add_argument("--epi-search-summary", help="rhdx_search_results.tsv")
    parser.add_argument("--species", help="Current run's species short code (e.g. bdbv, sudv)")
    parser.add_argument("--uniprotr-dir", help="UniProtR annotation output directory")
    parser.add_argument("--extractr-dir", help="UniProtExtractR annotation output directory")
    parser.add_argument("--rbioapi-dir", help="rbioapi annotation output directory")
    parser.add_argument("--auspice-json", help="Auspice JSON file for the phylogenetic tree")
    parser.add_argument("--iqtree", help="IQ-TREE2 treefile (Newick) to add to phylogenetic_trees")
    parser.add_argument("--query-data-dir", help="Directory of EXTRACT_QUERY_PROTEINS outputs (discovery.tsv, accessions.txt, query_proteins.fasta, etc.)")
    parser.add_argument("--hmm-dir", help="Directory of HMM_ANNOTATE outputs (*_hmm_*.txt)")
    parser.add_argument("--evidence-qc-dir", help="Directory containing evidence_qc/<species>/<domain>/ output tree")
    parser.add_argument("--schema-path", type=Path, default=Path(__file__).parent.parent / "database" / "knowledge_schema.sql", help="SQL schema file to load")
    parser.add_argument("--views-path", type=Path, default=Path(__file__).parent.parent / "database" / "knowledge_views.sql", help="Analytical views SQL file to load")
    parser.add_argument("--db-host", help="Host of an already-running shared PostgreSQL instance (started by START_KNOWLEDGE_DB). When set, this script connects to it instead of managing its own temporary server, and skips the final dump/stop (owned by STOP_KNOWLEDGE_DB).")
    parser.add_argument("--db-port", type=int, help="Port of the shared PostgreSQL instance (required with --db-host)")
    parser.add_argument("--species-screening-only", action="store_true", help="For the Nextclade species-assignment stage: only load species_assignments and the three Nextclade files per species (TSV, NWK, aligned FASTA)")
    return parser.parse_args()


def load_literature_evidence(conn, run_id, evidence_qc_dir):
    """Ingest clean literature evidence JSONs and per-domain QC reports."""
    evidence_qc_dir = Path(evidence_qc_dir)
    if not evidence_qc_dir.exists():
        print(f"  SKIP: evidence QC directory not found: {evidence_qc_dir}", file=sys.stderr)
        return
    report_files = sorted(evidence_qc_dir.rglob("qc_report.json"))
    if not report_files:
        print(f"  SKIP: no evidence QC reports found under {evidence_qc_dir}", file=sys.stderr)
        return
    with conn.cursor() as cur:
        for report_path in report_files:
            with open(report_path, encoding="utf-8") as f:
                report = json.load(f)
            species = normalize_text(report.get("species"))
            domain = normalize_text(report.get("domain"))
            if not species or not domain:
                continue
            expected_fields = report.get("expected_fields")
            summary = report.get("summary") or {}
            confidence_distribution = summary.get("confidence_distribution")
            field_coverage = summary.get("field_coverage")
            quote_presence_rate = parse_float(summary.get("quote_presence_rate"))
            duplicate_pmids = parse_int(summary.get("duplicate_pmids"))
            total_papers = parse_int(report.get("total_papers"))
            clean_count = parse_int(report.get("clean_count"))
            failed_count = parse_int(report.get("failed_count"))
            min_completeness = parse_float(report.get("min_completeness"))
            cur.execute(
                """
                INSERT INTO literature_domains (
                    run_id, species, domain, expected_fields, total_papers, clean_count, failed_count,
                    min_completeness, confidence_distribution, field_coverage, quote_presence_rate, duplicate_pmids
                ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (run_id, species, domain) DO UPDATE SET
                    expected_fields = EXCLUDED.expected_fields,
                    total_papers = EXCLUDED.total_papers,
                    clean_count = EXCLUDED.clean_count,
                    failed_count = EXCLUDED.failed_count,
                    min_completeness = EXCLUDED.min_completeness,
                    confidence_distribution = EXCLUDED.confidence_distribution,
                    field_coverage = EXCLUDED.field_coverage,
                    quote_presence_rate = EXCLUDED.quote_presence_rate,
                    duplicate_pmids = EXCLUDED.duplicate_pmids
                RETURNING domain_id
                """,
                (
                    run_id, species, domain,
                    json.dumps(expected_fields) if expected_fields is not None else None,
                    total_papers, clean_count, failed_count, min_completeness,
                    json.dumps(confidence_distribution) if confidence_distribution is not None else None,
                    json.dumps(field_coverage) if field_coverage is not None else None,
                    quote_presence_rate, duplicate_pmids,
                ),
            )
            domain_id = cur.fetchone()[0]
            for p in (report.get("per_paper") or []):
                if p.get("status") != "clean":
                    continue
                pmid = normalize_text(p.get("pmid"))
                if not pmid:
                    continue
                clean_path = report_path.parent / "clean" / f"{pmid}.json"
                if not clean_path.exists():
                    print(f"  WARNING: clean file not found for {pmid} in {domain}", file=sys.stderr)
                    continue
                with open(clean_path, encoding="utf-8") as f:
                    paper = json.load(f)
                title = normalize_text(paper.get("title"))
                authors = paper.get("authors")
                year = normalize_text(paper.get("year"))
                doi = normalize_text(paper.get("doi"))
                journal = normalize_text(paper.get("journal"))
                publication_date = normalize_text(paper.get("publication_date"))
                keywords = paper.get("keywords")
                pmcid = normalize_text(paper.get("pmcid"))
                cited_by_count = parse_int(paper.get("cited_by_count"))
                is_oa = paper.get("is_oa")
                qc_score = parse_float(paper.get("qc_score"))
                cur.execute(
                    """
                    INSERT INTO literature_papers (
                        domain_id, run_id, pmid, species, domain, title, authors, year, doi, journal,
                        publication_date, keywords, pmcid, cited_by_count, is_oa, qc_score, status
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT (run_id, pmid, species, domain) DO UPDATE SET
                        title = EXCLUDED.title,
                        authors = EXCLUDED.authors,
                        year = EXCLUDED.year,
                        doi = EXCLUDED.doi,
                        journal = EXCLUDED.journal,
                        publication_date = EXCLUDED.publication_date,
                        keywords = EXCLUDED.keywords,
                        pmcid = EXCLUDED.pmcid,
                        cited_by_count = EXCLUDED.cited_by_count,
                        is_oa = EXCLUDED.is_oa,
                        qc_score = EXCLUDED.qc_score,
                        status = EXCLUDED.status
                    RETURNING paper_id
                    """,
                    (
                        domain_id, run_id, pmid, species, domain, title,
                        json.dumps(authors) if authors is not None else None,
                        year, doi, journal, publication_date,
                        json.dumps(keywords) if keywords is not None else None,
                        pmcid, cited_by_count,
                        json.dumps(is_oa) if is_oa is not None else None,
                        qc_score, "clean",
                    ),
                )
                paper_id = cur.fetchone()[0]
                for extraction in paper.get("extraction", []):
                    field = normalize_text(extraction.get("field"))
                    if not field:
                        continue
                    present = bool(extraction.get("present"))
                    value = extraction.get("value")
                    quote = normalize_text(extraction.get("quote"))
                    confidence = normalize_text(extraction.get("confidence"))
                    cur.execute(
                        """
                        INSERT INTO literature_extractions (paper_id, field, present, value, quote, confidence)
                        VALUES (%s, %s, %s, %s, %s, %s)
                        ON CONFLICT (paper_id, field) DO UPDATE SET
                            present = EXCLUDED.present,
                            value = EXCLUDED.value,
                            quote = EXCLUDED.quote,
                            confidence = EXCLUDED.confidence
                        """,
                        (
                            paper_id, field, present,
                            json.dumps(value) if value is not None else None,
                            quote, confidence,
                        ),
                    )
    conn.commit()
    print(f"  Loaded literature evidence from {len(report_files)} domain report(s)", file=sys.stderr)


def main():
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    shared_mode = bool(args.db_host)

    if shared_mode:
        pg = TemporaryPostgres(outdir / "unused_data_dir", host=args.db_host, port=args.db_port)
    else:
        pg_data_dir = outdir / "postgres_data"
        pg_log_file = outdir / "postgres.log"
        pg = TemporaryPostgres(pg_data_dir, pg_log_file)

    try:
        if not shared_mode:
            print("Initializing PostgreSQL data directory...", file=sys.stderr)
            pg.init()
            print("Starting PostgreSQL...", file=sys.stderr)
            pg.start()
            print("Creating database...", file=sys.stderr)
            pg.create_db()

        with pg.connect() as conn:
            schema_path = Path(args.schema_path)
            views_path = Path(args.views_path)

            def _load_schema_and_run():
                print("Loading schema...", file=sys.stderr)
                load_schema(conn, schema_path)
                print("Loading analytical views...", file=sys.stderr)
                load_views(conn, views_path)
                print("Creating run record...", file=sys.stderr)
                create_run(conn, args.meta_id)

            if shared_mode:
                # Multiple species processes may connect concurrently; guard
                # schema creation with an advisory lock so they don't race on
                # CREATE TABLE against the same live database.
                with_advisory_lock(conn, SCHEMA_LOCK_KEY, _load_schema_and_run)
            else:
                _load_schema_and_run()

            results_dir = args.results_dir
            if not results_dir:
                # Fallback to sibling results of outdir
                results_dir = outdir.parent / "results"
            results_dir = Path(results_dir).resolve()
            bioinfo_dir = Path(args.bioinfo_dir).resolve() if args.bioinfo_dir else results_dir

            print("Registering pipeline outputs...", file=sys.stderr)
            register_pipeline_outputs(conn, args.meta_id, bioinfo_dir, stage="bioinformatics")

            if not args.species_screening_only:
                print("Registering phenotype annotation outputs...", file=sys.stderr)
                if args.uniprotr_dir:
                    register_dir_outputs(conn, args.meta_id, "phenotype_annotation", args.uniprotr_dir)
                if args.extractr_dir:
                    register_dir_outputs(conn, args.meta_id, "phenotype_annotation", args.extractr_dir)
                if args.rbioapi_dir:
                    register_dir_outputs(conn, args.meta_id, "phenotype_annotation", args.rbioapi_dir)
                if args.query_data_dir:
                    register_dir_outputs(conn, args.meta_id, "phenotype_annotation", args.query_data_dir)
                if args.hmm_dir:
                    register_dir_outputs(conn, args.meta_id, "phenotype_annotation", args.hmm_dir)

                print("Registering epidemiological data outputs...", file=sys.stderr)
                if args.epi_raw_dir:
                    register_dir_outputs(conn, args.meta_id, "epidemiological_data", args.epi_raw_dir)

                print("Registering auspice tree output...", file=sys.stderr)
                if args.auspice_json:
                    register_file_output(conn, args.meta_id, "bioinformatics", args.auspice_json)
                if args.iqtree:
                    register_file_output(conn, args.meta_id, "bioinformatics", args.iqtree)

                print("Loading reference genomes...", file=sys.stderr)
                load_reference_genomes_and_genes(conn, args.meta_id, results_dir)

                print("Registering nextstrain MSA fallback files...", file=sys.stderr)
                register_nextstrain_msa(conn, args.meta_id, results_dir, species=args.species)

            print("Loading samples...", file=sys.stderr)
            load_samples(
                conn, args.meta_id, results_dir,
                species_assignments=args.species_assignments,
                metadata_tsv=args.metadata_tsv,
                species=args.species,
                species_screening_only=args.species_screening_only,
            )

            if not args.species_screening_only:
                print("Loading phylogenetic trees and tips...", file=sys.stderr)
                load_trees_and_tips(conn, args.meta_id, results_dir, auspice_json=args.auspice_json, iqtree=args.iqtree, species=args.species)

                print("Loading mutations...", file=sys.stderr)
                load_mutations(conn, args.meta_id, results_dir)

                print("Loading phenotype annotations...", file=sys.stderr)
                load_phenotype_annotations(
                    conn, args.meta_id,
                    uniprotr_dir=args.uniprotr_dir,
                    extractr_dir=args.extractr_dir,
                    rbioapi_dir=args.rbioapi_dir,
                )

                print("Loading HMMER/Pfam annotations...", file=sys.stderr)
                load_hmm_annotations(conn, args.meta_id, hmm_dir=args.hmm_dir, species=args.species)

                print("Loading epidemiological data...", file=sys.stderr)
                load_epi_data(
                    conn, args.meta_id,
                    epi_raw_dir=args.epi_raw_dir,
                    epi_search_summary=args.epi_search_summary,
                    species=args.species,
                )

                print("Loading literature evidence...", file=sys.stderr)
                if args.evidence_qc_dir:
                    load_literature_evidence(conn, args.meta_id, args.evidence_qc_dir)

            print("Writing MultiQC summary...", file=sys.stderr)
            write_mqc_summary(conn, args.meta_id, outdir, args.prefix)

        if shared_mode:
            print("Shared DB mode: skipping per-species dump/stop (owned by STOP_KNOWLEDGE_DB).", file=sys.stderr)
        else:
            print("Dumping database...", file=sys.stderr)
            dump_database(pg, outdir, args.prefix)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
    finally:
        if not shared_mode:
            print("Stopping PostgreSQL...", file=sys.stderr)
            pg.stop()


if __name__ == "__main__":
    main()
