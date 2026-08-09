/*
 * Local module: ROUTE_PATHOGEN
 *
 * No-op passthrough process. Its only purpose is to give the pathogen
 * dispatch decision (which happens via a `.branch{}` channel operator in
 * PATHOGEN_ROUTER, and therefore leaves no trace in Nextflow's DAG export)
 * a real node in the execution graph. Every species group passes through
 * this process before being routed to its pathogen-specific workflow, so
 * diagram tools (e.g. nf-metro) render one "Route Pathogen" node with
 * edges fanning out to each pathogen workflow (Ebola now, others later).
 */

process ROUTE_PATHOGEN {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils"
    container null

    input:
    tuple val(meta), path(fasta), path(metadata)

    output:
    tuple val(meta), path(fasta), path(metadata), emit: routed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Routing '${meta.id}' (pathogen: ${meta.pathogen}, species: ${meta.species}) to its pathogen-specific workflow" >&2
    """
}
