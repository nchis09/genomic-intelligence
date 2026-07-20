#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * PGIRL setup / environment verification workflow.
 *
 * Usage:
 *   ./setup.sh
 *   # or directly:
 *   ./nextflow run setup.nf
 *
 * This workflow does two things:
 *   1. CHECK_ENV checks that Python packages and external tools are available.
 *   2. SETUP_DATABASE creates the pgirl PostgreSQL database (if needed), loads
 *      the schema, fetches the canonical Ebola reference proteomes, and syncs
 *      curated reference data.
 */

include { SETUP_ENV } from './modules/setup_env'
include { CHECK_ENV } from './modules/check_env'

process SETUP_DATABASE {
    tag "setup_database"

    output:
    stdout

    script:
    """
    bash "${projectDir}/scripts/setup_database.sh"
    """
}

workflow {
    SETUP_ENV().view()
    CHECK_ENV().view()
    SETUP_DATABASE().view()
}
