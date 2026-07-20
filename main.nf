#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * PGIRL end-to-end sample analysis workflow.
 *
 * Usage:
 *   nextflow run main.nf
 *   nextflow run main.nf --input_fasta input/input_FASTA.fasta --input_metadata input/input_metadata.csv
 *   nextflow run main.nf --db_url postgresql://localhost:5432/pgirl --use_llm false
 */

include { SAMPLE_ANALYSIS } from './workflows/sample_analysis'

workflow {
    SAMPLE_ANALYSIS(
        file(params.input_fasta),
        file(params.input_metadata),
        params.outdir,
        params.db_url,
        params.use_llm,
    )
}
