/*
 * Local module: LITERATURE_PDF
 *
 * Downloads full-text PDFs for the screened paper JSONs.
 * Tries Europe PMC OA, Unpaywall, and DOI resolver in order.
 */

process LITERATURE_PDF {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_pdf.yml"

    input:
    tuple val(meta), path("*.json")

    output:
    tuple val(meta), path("*.pdf"), optional: true, emit: pdfs
    path "pdf_download_summary.json", emit: summary

    when:
    !params.skip_literature_pdf && (task.ext.when == null || task.ext.when)

    script:
    def email_arg = params.pdf_email ? "--email \"${params.pdf_email}\"" : ""
    def timeout_arg = params.pdf_timeout ? "--timeout ${params.pdf_timeout}" : ""
    def sleep_arg = params.pdf_sleep ? "--sleep ${params.pdf_sleep}" : ""
    def retries_arg = params.pdf_retries ? "--retries ${params.pdf_retries}" : ""
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    python ${projectDir}/bin/fetch_literature_pdfs.py \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        ${email_arg} \
        ${timeout_arg} \
        ${sleep_arg} \
        ${retries_arg}
    """
}
