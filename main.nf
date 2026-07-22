#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * PGIRL Genomic Epidemic Intelligence System — Bioinformatics Pipeline
 *
 * Current scope: Nextclade classification + nextstrain/ebola phylogenetics.
 * Downstream stages (DB query, epi query, evidence integration) will be
 * re-enabled when ready.
 *
 * Usage:
 *   nextflow run main.nf
 *   nextflow run main.nf --input_fasta input/input_FASTA.fasta --input_metadata input/metadata.tsv
 *   nextflow run main.nf --outdir output
 */

include { SAMPLE_ANALYSIS } from './workflows/sample_analysis'

workflow {
    SAMPLE_ANALYSIS(
        file(params.input_fasta),
        file(params.input_metadata),
        params.outdir,
    )
}
