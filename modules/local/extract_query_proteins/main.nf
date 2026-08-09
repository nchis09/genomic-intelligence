/*
 * Local module: EXTRACT_QUERY_PROTEINS
 *
 * Dynamically discover UniProt accessions from the Auspice tree using:
 *   - Phylogenetic neighbors
 *   - Outbreak/strain matching
 *   - Shared AA mutations (per gene)
 *   - GenBank → UniProt cross-reference mapping
 * Downloads a UniProtKB TSV for downstream annotation.
 */

process EXTRACT_QUERY_PROTEINS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container null

    input:
    tuple val(meta), path(auspice_json), path(results_dir)

    output:
    tuple val(meta), path("*_discovery.tsv")            , emit: discovery
    tuple val(meta), path("*_all_accessions.txt")       , emit: accessions
    tuple val(meta), path("*_uniprot_download.tsv")     , emit: uniprot_tsv
    tuple val(meta), path("*_query_mutations.tsv")      , emit: mutations
    tuple val(meta), path("*_query_proteins.fasta")     , emit: proteins
    tuple val(meta), path("*_discovery_summary.json")   , emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def species = meta.species ?: meta.id.replaceAll(/^.*_/, '')
    def query   = params.query_samples ?: (meta.query_samples ?: '')
    """
    python3 ${projectDir}/bin/extract_query_proteins.py \\
        --auspice ${auspice_json} \\
        --results_dir ${results_dir} \\
        --query_samples "${query}" \\
        --species ${species} \\
        --prefix ${prefix}
    """
}
