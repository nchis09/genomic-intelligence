/*
 * Local module: BUILD_KNOWLEDGE_DB
 *
 * Build a temporary PostgreSQL knowledge warehouse from key pipeline outputs.
 */

process BUILD_KNOWLEDGE_DB {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_knowledge.yml"

    input:
    tuple val(meta), path(species_assignments), path(metadata_tsv), path(epi_raw_dir), path(epi_search_summary)

    output:
    tuple val(meta), path("knowledge_warehouse"), emit: knowledge_db

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = ["--outdir knowledge_warehouse", "--meta-id ${meta.id}", "--prefix ${prefix}"]
    if (species_assignments.name != 'NO_FILE') args << "--species-assignments ${species_assignments}"
    if (metadata_tsv.name != 'NO_FILE')          args << "--metadata-tsv ${metadata_tsv}"
    if (epi_raw_dir.name != 'NO_FILE')           args << "--epi-raw-dir ${epi_raw_dir}"
    if (epi_search_summary.name != 'NO_FILE')    args << "--epi-search-summary ${epi_search_summary}"
    // Register every published pipeline output in the provenance table without staging the whole tree.
    args << "--results-dir ${params.outdir}"
    """
    mkdir -p knowledge_warehouse

    python3 ${projectDir}/bin/build_knowledge_db.py ${args.join(' ')}
    """
}
