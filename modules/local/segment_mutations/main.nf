/*
 * Local module: SEGMENT_MUTATIONS
 *
 * Produce mutations_{segment}.tsv for every genome segment a flu build
 * actually produced — query AND background tips — so phenotype annotation
 * (EXTRACT_QUERY_PROTEINS) gets reference-relative AA mutations without
 * any modification to the vendored Nextstrain workflows under data/
 * (gitignored; fresh clones arrive clean).
 *
 * Per segment, bin/extract_segment_mutations.py prefers `nextclade run`
 * on the build's own aligned.fasta (seasonal: staged reference.fasta +
 * genemap.gff, or the ds_{segment} datasets fetched for internals) and
 * falls back to deriving ref-relative mutations from the segment's
 * Auspice JSON branch attributes (avian, augur-based).
 */

process SEGMENT_MUTATIONS {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"
    container null

    input:
    tuple val(meta), path(results_dir), path(auspice_jsons)

    output:
    tuple val(meta), path("mutations_*.tsv", optional: true), emit: mutations_tsvs
    tuple val(meta), path("segment_mutations_report.txt")   , emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def jsons   = auspice_jsons instanceof List ? auspice_jsons : [ auspice_jsons ]
    def js_args = jsons
        .findAll { it && !it.name.startsWith('NO_FILE') }
        .collect { it.toString() }
        .join(' ')
    """
    python3 ${projectDir}/bin/extract_segment_mutations.py \\
        --results_dir "${results_dir}" \\
        --pathogen "${meta.pathogen ?: ''}" \\
        --lineage "${meta.species ?: ''}" \\
        --outdir . \\
        --report segment_mutations_report.txt \\
        ${js_args ? "--auspice ${js_args}" : ''}
    """
}
