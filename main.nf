#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pgirl/genomic-intelligence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/pgirl/genomic-intelligence
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMIC_INTELLIGENCE  } from './workflows/genomic-intelligence'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_genomic-intelligence_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_genomic-intelligence_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow PGIRL_GENOMIC_INTELLIGENCE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GENOMIC_INTELLIGENCE (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = channel.empty() // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    println(new File('assets/pipeline_banner.txt').text)

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    PGIRL_GENOMIC_INTELLIGENCE (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        PGIRL_GENOMIC_INTELLIGENCE.out.multiqc_report
    )

    // Safety-net stop/dump for the shared knowledge-warehouse PostgreSQL
    // instance in case the pipeline failed or was interrupted before the
    // in-DAG EXPORT_KNOWLEDGE_DB -> STOP_KNOWLEDGE_DB steps (see
    // subworkflows/local/knowledge_warehouse/main.nf) could run. On a normal
    // successful completion the server is already stopped by then, so this
    // checks pg_ctl status first and exits quietly instead of logging a
    // false-failure warning.
    workflow.onComplete = {
        if (params.skip_knowledge_warehouse) {
            return
        }
        def data_dir = file("${params.kw_data_dir}")
        if (!data_dir.exists()) {
            return
        }
        def dump_dir = file("${params.outdir}/knowledge_warehouse")
        dump_dir.mkdirs()
        def cmd = "set -e; for d in ${projectDir}/work/conda/env-*/; do if [ -f \"\${d}bin/pg_ctl\" ]; then export PATH=\"\${d}bin:\$PATH\"; if ! pg_ctl -D \"${data_dir}\" status >/dev/null 2>&1; then exit 0; fi; pg_dump -h ${params.kw_db_host} -p ${params.kw_db_port} -U postgres -f \"${dump_dir}/run_genomic_intelligence.sql\" genomic_intelligence; pg_ctl -D \"${data_dir}\" -m fast stop || true; break; fi; done"
        def proc = ["/bin/bash", "-c", cmd.toString()].execute()
        proc.waitFor()
        if (proc.exitValue() != 0) {
            log.warn "Knowledge warehouse safety-net stop hook failed: ${proc.err.text}"
        } else {
            log.info "Knowledge warehouse safety-net check complete."
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
