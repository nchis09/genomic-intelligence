/*
 * Local module: PUBMED_METADATA
 *
 * Fetches PubMed metadata (title, abstract, authors, journal, year, DOI,
 * keywords, PMCID) for each paper in a LITERATURE_SEARCH results.json.
 * Writes one JSON file per PMID.
 */

process PUBMED_METADATA {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    // Use a dedicated conda env; Nextflow builds and activates it.
    conda "${projectDir}/envs/pgirl_pubmed.yml"

    input:
    tuple val(meta), path("results.json")

    output:
    tuple val(meta), path("*.json"), emit: metadata

    when:
    !params.skip_pubmed_metadata && (task.ext.when == null || task.ext.when)

    script:
    def email = params.pubmed_email ? "--email ${params.pubmed_email}" : ""
    def batch_size = params.pubmed_batch_size ?: 100
    def sleep = params.pubmed_sleep ?: 0.5
    """
    python ${projectDir}/bin/fetch_pubmed_metadata.py \
        --results-json results.json \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --outdir . \
        ${email} \
        --batch-size ${batch_size} \
        --sleep ${sleep}
    """
}
