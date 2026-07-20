#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * PGIRL reference database update workflow.
 *
 * Usage:
 *   nextflow run db_update.nf
 *   nextflow run db_update.nf --db_url postgresql://localhost:5432/pgirl
 *   nextflow run db_update.nf --sources variants genotype_phenotype gene_function --dry_run true
 */

include { DB_UPDATE } from './workflows/db_update'

workflow {
    DB_UPDATE(
        params.db_url,
        params.pathogens,
        params.sources,
        params.dry_run,
    )
}
