# PostgreSQL / pgirl Database Navigation Cheatsheet

Essential commands for inspecting the database without manual browsing — useful when the database grows too large to check by hand.

---

## Connection

```bash
# Use the local PostgreSQL 16 client
/opt/homebrew/opt/postgresql@16/bin/psql pgirl

# Or with full connection string
/opt/homebrew/opt/postgresql@16/bin/psql postgresql://localhost:5432/pgirl

# Exit psql
\q
```

---

## Listing & Discovery

```sql
-- List all tables
\dt

-- List all views
\dv

-- List all materialized views
\dm

-- List all indexes
\di

-- List all indexes for a specific table
\di public.protein_variants

-- Show table / view definition
\d protein_variants
\d+ protein_variants          -- includes comments and storage

-- Show a view's underlying query
\d+ v_mutation_intelligence

-- List all schemas
\dn

-- List all functions
\df

-- List all sequences
\ds
```

---

## Database Size & Table Sizes

```sql
-- Total database size
SELECT pg_size_pretty(pg_database_size('pgirl')) AS db_size;

-- Size of every table + indexes, sorted largest first
SELECT
    schemaname,
    relname                     AS table_name,
    pg_size_pretty(pg_total_relation_size(relid)) AS total_size
FROM pg_stat_user_tables
WHERE schemaname = 'public'
ORDER BY pg_total_relation_size(relid) DESC;

-- Size of each table, toast, and indexes separately
SELECT
    relname AS table_name,
    pg_size_pretty(pg_relation_size(relid))             AS table_size,
    pg_size_pretty(pg_indexes_size(relid))              AS indexes_size,
    pg_size_pretty(pg_total_relation_size(relid))       AS total_size
FROM pg_stat_user_tables
WHERE schemaname = 'public'
ORDER BY pg_total_relation_size(relid) DESC;

-- Index sizes
SELECT
    indexrelname AS index_name,
    relname      AS table_name,
    pg_size_pretty(pg_relation_size(indexrelid)) AS index_size
FROM pg_stat_user_indexes
WHERE schemaname = 'public'
ORDER BY pg_relation_size(indexrelid) DESC;
```

---

## Index Inspection

```sql
-- All indexes in the database
SELECT
    tablename,
    indexname,
    indexdef
FROM pg_indexes
WHERE schemaname = 'public'
ORDER BY tablename, indexname;

-- Indexes for one table
SELECT indexname, indexdef
FROM pg_indexes
WHERE schemaname = 'public' AND tablename = 'protein_variants';

-- Count indexes per table
SELECT tablename, count(*) AS index_count
FROM pg_indexes
WHERE schemaname = 'public'
GROUP BY tablename
ORDER BY index_count DESC;

-- Check whether a specific column is indexed
SELECT tablename, indexname
FROM pg_indexes
WHERE schemaname = 'public'
  AND indexdef ILIKE '%(species_id, gene, "position")%';

-- GIN indexes only
SELECT tablename, indexname, indexdef
FROM pg_indexes
WHERE schemaname = 'public'
  AND indexdef ILIKE '%USING gin%';

-- Partial indexes only (those with a WHERE clause)
SELECT tablename, indexname, indexdef
FROM pg_indexes
WHERE schemaname = 'public'
  AND indexdef ILIKE '%WHERE%';

-- Index usage statistics (how often each index is used)
SELECT
    indexrelname AS index_name,
    relname      AS table_name,
    idx_scan,
    idx_tup_read,
    idx_tup_fetch
FROM pg_stat_user_indexes
WHERE schemaname = 'public'
ORDER BY idx_scan DESC;

-- Tables that may be missing indexes (large tables with lots of sequential scans)
SELECT
    relname AS table_name,
    seq_scan,
    seq_tup_read,
    idx_scan,
    idx_tup_fetch
FROM pg_stat_user_tables
WHERE schemaname = 'public'
ORDER BY seq_scan DESC;
```

---

## Table Row Counts & Structure

```sql
-- Row counts for all tables (approximate, fast)
SELECT
    relname AS table_name,
    n_live_tup AS row_count
FROM pg_stat_user_tables
WHERE schemaname = 'public'
ORDER BY n_live_tup DESC;

-- Exact row count for a specific table
SELECT count(*) FROM protein_variants;

-- Columns of a table
SELECT column_name, data_type, is_nullable
FROM information_schema.columns
WHERE table_schema = 'public' AND table_name = 'protein_variants'
ORDER BY ordinal_position;

-- Foreign keys and their referenced tables/columns
SELECT
    tc.table_name,
    kcu.column_name,
    ccu.table_name  AS foreign_table_name,
    ccu.column_name AS foreign_column_name
FROM information_schema.table_constraints tc
JOIN information_schema.key_column_usage kcu
    ON tc.constraint_name = kcu.constraint_name
JOIN information_schema.constraint_column_usage ccu
    ON tc.constraint_name = ccu.constraint_name
WHERE tc.constraint_type = 'FOREIGN KEY'
  AND tc.table_schema = 'public'
ORDER BY tc.table_name, kcu.column_name;

-- Primary keys
SELECT
    tc.table_name,
    kcu.column_name
FROM information_schema.table_constraints tc
JOIN information_schema.key_column_usage kcu
    ON tc.constraint_name = kcu.constraint_name
WHERE tc.constraint_type = 'PRIMARY KEY'
  AND tc.table_schema = 'public'
ORDER BY tc.table_name;
```

---

## Query Performance

```sql
-- Explain the plan for a query (does not run)
EXPLAIN SELECT * FROM protein_variants WHERE species_id = 'EBOV' AND gene = 'GP' AND "position" = 50;

-- Explain with actual execution times and row counts
EXPLAIN (ANALYZE, BUFFERS, FORMAT TEXT)
SELECT * FROM protein_variants
WHERE species_id = 'EBOV' AND gene = 'GP' AND "position" = 50;

-- Check if the query uses an index
EXPLAIN (ANALYZE, BUFFERS)
SELECT * FROM gene_function WHERE protein_function_tsv @@ plainto_tsquery('receptor binding');

-- Most time-consuming queries since stats reset
SELECT
    substring(query, 1, 80) AS query_preview,
    calls,
    round(total_exec_time::numeric, 2) AS total_ms,
    round(mean_exec_time::numeric, 2)  AS mean_ms,
    round(rows::numeric / calls, 2)    AS avg_rows
FROM pg_stat_statements
ORDER BY total_exec_time DESC
LIMIT 20;
```

---

## Gene Function Intelligence Queries

```sql
-- Basic protein function lookup
SELECT species_id, gene, protein_name, protein_length_aa, protein_function
FROM gene_function
WHERE species_id = 'BDBV' AND gene = 'GP';

-- Search protein function text with full-text search
SELECT species_id, gene, protein_name
FROM gene_function
WHERE protein_function_tsv @@ plainto_tsquery('receptor binding');

-- Find genes that have a domain at a specific position
SELECT species_id, gene, protein_name
FROM gene_function
WHERE key_domains @> '[{"start": 50}]';

-- Find entries that cite a specific PubMed ID
SELECT species_id, gene, protein_name
FROM gene_function
WHERE literature_refs && '{PMID:11836430}';

-- Find entries with a PDB structure
SELECT species_id, gene, pdb_ids
FROM gene_function
WHERE pdb_ids && '{5JQ3}';
```

---

## Mutation / Variant Surveillance Queries

```sql
-- All mutations at a specific position
SELECT * FROM mutations
WHERE species_id = 'EBOV' AND protein = 'GP' AND "position" = 50;

-- Public-health-relevant mutations only
SELECT * FROM mutations
WHERE species_id = 'EBOV' AND protein = 'GP' AND public_health_relevant = true;

-- Variant observation with full context
SELECT * FROM v_genomic_intelligence_snapshot
WHERE species_id = 'EBOV' AND gene = 'GP' AND "position" = 50;

-- Pre-computed intelligence snapshot
SELECT * FROM mv_mutation_intelligence
WHERE species_id = 'EBOV' AND protein = 'GP' AND "position" = 50;

-- Count variants per gene / species
SELECT species_id, gene, count(*) AS variant_count
FROM protein_variants
GROUP BY species_id, gene
ORDER BY species_id, gene;
```

---

## Refreshing Materialized Views

```sql
-- Refresh a single materialized view
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_intelligence;

-- Refresh all intelligence mat views
SELECT refresh_intelligence_matviews();

-- Refresh all pgirl mat views at once
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_summary;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_with_phenotype;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_frequency;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_lineage_breakdown;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_lineage_country;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_phenotype_surveillance;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_phenotype_geo_temporal;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_position_variants;
REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_intelligence;
```

---

## Exporting Data

```sql
-- Export a query to CSV
\COPY (SELECT * FROM gene_function ORDER BY species_id, gene)
TO '/tmp/gene_function_export.csv' CSV HEADER;

-- Export a table to CSV
\COPY gene_function TO '/tmp/gene_function_export.csv' CSV HEADER;

-- From bash (no need to be in psql)
/opt/homebrew/opt/postgresql@16/bin/psql pgirl -P format=unaligned -P tuples_only=on \
  -c "SELECT * FROM gene_function WHERE species_id='EBOV'" > /tmp/ebov_gene_function.txt
```

---

## Maintenance & Health

```sql
-- Long-running queries / locks
SELECT
    pid,
    state,
    now() - query_start AS duration,
    substring(query, 1, 100) AS query_preview
FROM pg_stat_activity
WHERE state != 'idle'
ORDER BY duration DESC;

-- Current locks
SELECT
    l.locktype,
    l.relation::regclass,
    l.mode,
    l.granted,
    a.query
FROM pg_locks l
JOIN pg_stat_activity a ON l.pid = a.pid
WHERE a.state != 'idle';

-- Vacuum / bloat overview
SELECT
    schemaname,
    relname,
    n_dead_tup,
    n_live_tup,
    last_vacuum,
    last_autovacuum,
    last_analyze,
    last_autoanalyze
FROM pg_stat_user_tables
WHERE schemaname = 'public'
ORDER BY n_dead_tup DESC;

-- Update planner statistics for all tables
ANALYZE;
ANALYZE protein_variants;
ANALYZE gene_function;

-- VACUUM (reclaim space) — avoid on production during peak hours
VACUUM ANALYZE protein_variants;
```

---

## Quick One-Liners

```bash
# Database size
/opt/homebrew/opt/postgresql@16/bin/psql pgirl -c "SELECT pg_size_pretty(pg_database_size('pgirl'));"

# Number of indexes
/opt/homebrew/opt/postgresql@16/bin/psql pgirl -c "SELECT count(*) FROM pg_indexes WHERE schemaname='public';"

# Number of materialized views
/opt/homebrew/opt/postgresql@16/bin/psql pgirl -c "SELECT count(*) FROM pg_matviews WHERE schemaname='public';"

# Largest 5 tables
/opt/homebrew/opt/postgresql@16/bin/psql pgirl -P format=unaligned -P tuples_only=on \
  -c "SELECT relname, pg_size_pretty(pg_total_relation_size(relid)) FROM pg_stat_user_tables WHERE schemaname='public' ORDER BY pg_total_relation_size(relid) DESC LIMIT 5;"

# Refresh all mat views
/opt/homebrew/opt/postgresql@16/bin/psql pgirl -c "SELECT refresh_intelligence_matviews();"
```

---

## File Locations

- SQL optimization script: `database/optimizations.sql`
- db_sync script: `scripts/db_sync.py`
- Seeding script: `database/ebola/gene_function/seed_gene_function_from_uniprot.py`
- Config: `config.py`
