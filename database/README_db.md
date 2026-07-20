# PGIRL Reference Database — Setup & Usage

PostgreSQL reference database for the Genomic Epidemic Intelligence System.

## Prerequisites

```bash
pip install psycopg2-binary pyyaml
```

PostgreSQL 14+ must be installed and running.

## One-time setup

```bash
# 1. Create the database
createdb pgirl

# 2. Apply the current schema
psql pgirl < database/schema.sql

# 3. Fetch Ebola reference proteomes
/Users/christianndekezi/anaconda3/bin/python3 database/ebola/protein_variants/fetch_reference_proteomes.py

# 4. Run the variant calling pipeline (all 6 species)
/Users/christianndekezi/anaconda3/bin/python3 -u database/ebola/protein_variants/call_variants.py --batch-size 50

# 5. Assign lineages to Ebola genomes with Nextclade
#    (requires Nextclade CLI installed:
#     /Users/christianndekezi/anaconda3/bin/conda install -c bioconda nextclade -y)
#    Note: Nextclade currently provides clade/lineage information only for
#    Zaire ebolavirus (EBOV). SUDV, BDBV, RESTV, TAFV and BOMV are skipped
#    or reported as unassigned because the datasets lack clade definitions.
/Users/christianndekezi/anaconda3/bin/python3 database/ebola/protein_variants/assign_lineages.py

# 6. Create surveillance views for the intelligence engine
/opt/homebrew/opt/postgresql@16/bin/psql -d pgirl -f database/views/surveillance_views.sql
```

## Reload / reset Ebola variant data (during development)

```bash
# Clear existing Ebola variant data
/opt/homebrew/opt/postgresql@16/bin/psql -d pgirl -c "
DELETE FROM protein_variants WHERE pathogen_id='ebola';
DELETE FROM genome_metadata WHERE pathogen_id='ebola';
"

# Re-run pipeline
/Users/christianndekezi/anaconda3/bin/python3 -u database/ebola/protein_variants/call_variants.py --batch-size 50
```

## Incremental update (only new genomes)

```bash
/Users/christianndekezi/anaconda3/bin/python3 -u database/ebola/protein_variants/call_variants.py --batch-size 50 --skip-existing
```

## Surveillance views for the intelligence engine

After the variant data is loaded, create the surveillance views that the
intelligence engine queries when interpreting a new genome:

```bash
/opt/homebrew/opt/postgresql@16/bin/psql -d pgirl -f database/views/surveillance_views.sql
```

These views turn raw `protein_variants` + `genome_metadata` into epidemiological
intelligence:

```text
┌─────────────────────┐     ┌─────────────────────┐
│   genome_metadata   │     │   protein_variants  │
│  (who, where, when) │     │  (what changed)     │
└──────────┬──────────┘     └──────────┬──────────┘
           │                             │
           └─────────────┬───────────────┘
                         ▼
           ┌─────────────────────────────┐
           │   v_genome_variant_profile    │  ← base view: every variant + context
           └─────────────┬───────────────┘
                         │
     ┌───────────────────┼───────────────────┐
     │                   │                   │
     ▼                   ▼                   ▼
┌────────────┐  ┌──────────────┐  ┌─────────────────┐
│v_mutation_ │  │ v_mutation_  │  │ v_mutation_     │
│surveillance│  │  geography   │  │    trends       │
│(by country, │  │(first/last/  │  │(yearly counts) │
│  year)     │  │  countries)  │  │                 │
└────────────┘  └──────────────┘  └─────────────────┘
     │
     ▼
┌──────────────────────────────────────────────────────┐
│         v_mutation_frequency                         │
│  (country frequency vs global frequency)             │
└──────────────────────────────────────────────────────┘

Other derived views:
  v_mutation_with_phenotype  → joins curated evidence
  v_mutation_co_occurrence   → mutation pairs in same genome
  v_variant_recent_detection → latest observation per mutation
  v_species_summary          → denominator/context counts
```

| View | Purpose |
|---|---|
| `v_genome_variant_profile` | Every variant joined to its genome metadata |
| `v_mutation_surveillance` | Per-mutation counts by species, country, year |
| `v_mutation_geography` | First/last seen date and list of countries per mutation |
| `v_mutation_trends` | Yearly counts per mutation for emergence detection |
| `v_mutation_frequency` | Country-level and global frequency of each mutation |
| `v_mutation_with_phenotype` | Links observed mutations to curated genotype-phenotype evidence |
| `v_mutation_co_occurrence` | Pairs of mutations seen in the same genome |
| `v_species_summary` | Total genomes per species/country/year |
| `v_variant_recent_detection` | Most recent observation of each mutation |
| `v_mutation_lineage_breakdown` | Per-mutation counts split by lineage (requires Nextclade) |
| `v_mutation_lineage_surveillance` | Per-mutation counts split by lineage, country and year |
| `v_phenotype_surveillance` | Aggregate inventory of genotype-phenotype candidates from literature |
| `v_phenotype_geo_temporal` | Where/when each phenotype has been observed, by country, lineage, and year |

### Fast lookup materialized views

For the intelligence engine, the following materialized views pre-aggregate the
above so queries return in milliseconds:

| Materialized view | Source | Purpose |
|---|---|---|
| `mv_mutation_summary` | `v_mutation_geography` | Fast mutation existence + geography |
| `mv_mutation_frequency` | `v_mutation_frequency` | Fast frequency lookups |
| `mv_position_variants` | `v_genome_variant_profile` | All substitutions per position |
| `mv_mutation_lineage_breakdown` | `v_mutation_lineage_breakdown` | Per-mutation counts by lineage |
| `mv_mutation_lineage_country` | `v_mutation_lineage_surveillance` | Per-mutation counts by lineage, country and year |
| `mv_phenotype_surveillance` | `v_phenotype_surveillance` | Fast inventory of genotype-phenotype candidates |
| `mv_mutation_with_phenotype` | `v_mutation_with_phenotype` | Fast mutation → phenotype lookups for observed variants |
| `mv_phenotype_geo_temporal` | `v_phenotype_geo_temporal` | Fast phenotype × region × lineage × year lookups |

Refresh them after every variant or genotype-phenotype load:

```bash
/opt/homebrew/opt/postgresql@16/bin/psql -d pgirl -c "
REFRESH MATERIALIZED VIEW mv_mutation_summary;
REFRESH MATERIALIZED VIEW mv_mutation_frequency;
REFRESH MATERIALIZED VIEW mv_position_variants;
REFRESH MATERIALIZED VIEW mv_mutation_lineage_breakdown;
REFRESH MATERIALIZED VIEW mv_mutation_lineage_country;
REFRESH MATERIALIZED VIEW mv_phenotype_surveillance;
REFRESH MATERIALIZED VIEW mv_mutation_with_phenotype;
REFRESH MATERIALIZED VIEW mv_phenotype_geo_temporal;
"
```

`scripts/db_sync.py` refreshes these automatically after the variants source
runs.

Example queries:

```sql
-- Has GP:A82V been seen before? Where?
SELECT * FROM v_mutation_geography
WHERE pathogen_id = 'ebola' AND species_id = 'EBOV'
  AND gene = 'GP' AND position = 82 AND alt_aa = 'V';

-- Frequency of a mutation in Sierra Leone
SELECT * FROM v_mutation_frequency
WHERE pathogen_id = 'ebola' AND species_id = 'EBOV'
  AND gene = 'GP' AND position = 82 AND alt_aa = 'V'
  AND collection_country_code = 'SLE';

-- Trend of a mutation over time
SELECT collection_year, genome_count FROM v_mutation_trends
WHERE pathogen_id = 'ebola' AND species_id = 'EBOV'
  AND gene = 'GP' AND position = 82 AND alt_aa = 'V'
ORDER BY collection_year;
```

## How the intelligence engine queries the DB

Every query is **pathogen-scoped** then **species-scoped**. Once the upstream
bioinformatics identifies e.g. `Sudan ebolavirus`, all queries add:

```sql
WHERE pathogen_id = 'ebola' AND species_id = 'SUDV'
```

This means a `GP:A82V` mutation in Ebola never collides with any same-position
mutation in Dengue or Influenza — they are completely separate rows in separate
species arms.

### Key queries the intelligence engine uses

**1. Has this species been in this location before?**
```sql
SELECT * FROM v_geographic_lookup
WHERE species_id = 'SUDV' AND country = 'Uganda';
```

**2. Has this specific mutation been reported before? What does it mean?**
```sql
SELECT * FROM v_mutation_intelligence
WHERE pathogen_id = 'ebola' AND species_id = 'SUDV'
  AND hgvs_protein = 'GP:A82V';
```

**3. Is this a novel mutation (not in the DB)?**
```sql
SELECT COUNT(*) FROM mutations
WHERE pathogen_id = 'ebola' AND species_id = 'SUDV'
  AND hgvs_protein = 'GP:T230A';
-- returns 0 → novel, flag for expert review
```

**4. What outbreaks has this lineage caused?**
```sql
SELECT * FROM v_lineage_outbreak_summary
WHERE lineage_id = 'SUDV-MUBENDE';
```

**5. What gene domain does this mutation sit in?**
```sql
SELECT key_domains FROM gene_function
WHERE species_id = 'SUDV' AND gene = 'VP35';
-- Parse JSONB to find which domain covers position 230
```

**6. Full epidemiology for an outbreak (cases, CFR, country health system, disease params):**
```sql
SELECT * FROM v_disease_epidemiology
WHERE pathogen_id = 'ebola' AND species_id = 'EBOV'
ORDER BY start_date DESC LIMIT 10;
```

**7. Aggregate outbreak burden per country:**
```sql
SELECT * FROM v_country_outbreak_burden
ORDER BY total_cases_all_outbreaks DESC NULLS LAST;
```

## Adding a new pathogen (e.g. Dengue)

1. Add a row to `pathogens` (`pathogen_id = 'dengue'`)
2. Add rows to `species` (DENV-1 through DENV-4)
3. Add reference genomes, outbreaks, lineages, mutations, genotype_phenotype
   — all with `pathogen_id = 'dengue'`
4. The schema does not change — just new rows in existing tables

## Genotype-phenotype extraction

Candidate mutation-phenotype associations are extracted from PubMed abstracts
using the NCBI PubTator API. They are inserted as `unverified` rows in
`genotype_phenotype` with `record_flagged = True` for human review.

```bash
# Extract and verify (run directly)
/Users/christianndekezi/anaconda3/bin/python3 database/ebola/genotype_phenotype/extract_from_pubtator.py --max-results 50
/Users/christianndekezi/anaconda3/bin/python3 database/ebola/genotype_phenotype/verify_extracted_pmids.py

# Or via db_sync.py
python scripts/db_sync.py --pathogen ebola --sources genotype_phenotype
```

The `genotype_phenotype` table now captures mutations, motifs, and clades linked
to phenotypes such as immune escape, vaccine effectiveness, disease severity,
virulence, transmission, and drug resistance.

Clade/strain names extracted by PubTator are looked up against `lineages.known_aliases`
(e.g. `Makona` → `EBOV-Ebov-2013`, `Mayinga` → `EBOV-Ebov-1976`, `Kikwit` → `EBOV-Ebov-1995`).
When a match is found, `lineage_id` is filled; otherwise `lineage_id` is left `NULL`
and the annotation stays at the species level (`EBOV`).

## API-driven updates

The scripts in `database/api_sources/` (`outbreak_api.py`, `literature_api.py`)
fetch new records from external sources and INSERT/UPDATE rows in this database.
The `api_refresh_log` table tracks each run.

Use `scripts/db_sync.py` to orchestrate a sync:
```bash
# Default: sync protein variants only (incremental — skips existing genomes)
python scripts/db_sync.py --pathogen ebola

# Sync outbreaks + variants (outbreaks requires WHO DON API access)
python scripts/db_sync.py --pathogen ebola --sources outbreaks variants

# Sync literature-derived genotype-phenotype candidates (PubTator + PMID verification)
python scripts/db_sync.py --pathogen ebola --sources genotype_phenotype

# Sync everything
python scripts/db_sync.py --pathogen ebola --sources all
```

## Schema evolution

`database/schema.sql` is the single source of truth for the current schema.
During development it is regenerated from the working database when the schema
changes. There are no incremental migration files.

## Table overview

| Table | Purpose |
|---|---|
| `pathogens` | One row per pathogen group (ebola, dengue...) |
| `species` | One row per virus species (EBOV, SUDV, DENV-1...) |
| `reference_genomes` | RefSeq canonical references + gene coordinates |
| `outbreaks` | All historical outbreak events (now includes `lineage_id` FK + `country_code`) |
| `lineages` | Named lineages/clades with geographic distribution; `known_aliases` maps common strain names (e.g. `Makona`) to the canonical `lineage_id` |
| `mutations` | All curated mutation observations |
| `genotype_phenotype` | Mutation → phenotype associations with evidence |
| `gene_function` | Gene biology, domains, hotspots per species |
| `disease_epidemiology` | Disease-level parameters (R0, incubation, CFR, transmission routes) per pathogen × species |
| `country_context` | WHO GHO health-system indicators per country (demographics, doctors, beds, mortality) |
| `api_refresh_log` | Tracks API update runs |

**Removed tables** (migration 003):
- `outbreak_mutations` — redundant with `mutations.reported_in_outbreaks` and `outbreaks.key_mutations`
- `lineage_outbreaks` — replaced by `outbreaks.lineage_id` FK
- `geographic_history` — redundant with `outbreaks` (view queries outbreaks directly)

## Views

| View | Returns |
|---|---|
| `v_mutation_intelligence` | Full mutation + all phenotype associations + gene context |
| `v_geographic_lookup` | Aggregated geographic occurrence per species/lineage/country (queries outbreaks directly) |
| `v_lineage_outbreak_summary` | All outbreaks for a lineage with full metadata (joins via outbreaks.lineage_id) |
| `v_disease_epidemiology` | Outbreak enriched with country health-system context + disease parameters + lineage |
| `v_country_outbreak_burden` | Aggregate outbreak burden per country with health-system context |
| `v_country_context_latest` | Latest value per country per WHO GHO indicator |
| `v_curator_queue` | Records awaiting human verification (from migration 001) |
