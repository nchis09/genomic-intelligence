/*
 * Local module: LITERATURE_GROBID
 *
 * Converts downloaded PDFs to TEI XML using GROBID.
 */

process LITERATURE_GROBID {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_grobid.yml"

    input:
    tuple val(meta), path("*.pdf")

    output:
    tuple val(meta), path("*.tei.xml"), optional: true, emit: xml
    path "grobid_summary.json", emit: summary

    when:
    !params.skip_literature_grobid && (task.ext.when == null || task.ext.when)

    script:
    def url_arg = params.grobid_url ? "--grobid-url ${params.grobid_url}" : ""
    def timeout_arg = params.grobid_timeout ? "--timeout ${params.grobid_timeout}" : ""
    def sleep_arg = params.grobid_sleep ? "--sleep ${params.grobid_sleep}" : ""
    def retries_arg = params.grobid_retries ? "--retries ${params.grobid_retries}" : ""
    def batch_arg = params.grobid_batch_size ? "--batch-size ${params.grobid_batch_size}" : ""
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    python ${projectDir}/bin/convert_literature_pdfs_to_xml.py \
        --input-dir . \
        --outdir . \
        ${url_arg} \
        ${timeout_arg} \
        ${sleep_arg} \
        ${retries_arg} \
        ${batch_arg}
    """
}
