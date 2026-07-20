/*
 * Reference database update workflow.
 *
 * Wraps scripts/db_sync.py to fetch/update reference data from external APIs
 * (NCBI, Nextstrain) and curated sources into the PostgreSQL PGIRL database.
 *
 * This is intentionally kept separate from SAMPLE_ANALYSIS so DB curation is not
 * rerun for every sample.
 */

include { DB_SYNC } from '../modules/db_sync'

workflow DB_UPDATE {
    take:
    db_url
    pathogens
    sources
    dry_run

    main:
    DB_SYNC(db_url, pathogens, sources, dry_run)
}
