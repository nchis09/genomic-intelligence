/*
 * Local module: GEOGRAPHIC_SUMMARY
 *
 * Pre-generate the per-species geographic spread narrative ("how the strain
 * spread, where the query genomes most likely came from") from the
 * transmission-context TSVs, so the dashboard renders text instantly instead
 * of calling Ollama at view time. Falls back to a deterministic template —
 * geo_summary.json is always produced.
 */

process GEOGRAPHIC_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_pathogen_genomics.yml"

    input:
    tuple val(meta), path(tc_dir)
    path llm_note_r

    output:
    tuple val(meta), path("geo_summary.json"), emit: json

    when:
    !params.skip_transmission_context && (task.ext.when == null || task.ext.when)

    script:
    def species = meta.species ?: meta.id
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"
    export PG_OLLAMA_HOST="${params.ollama_host}"
    export PG_OLLAMA_MODEL="${params.ollama_model}"

    Rscript ${projectDir}/bin/summarize_geographic.R \\
        --tc_dir ${tc_dir} \\
        --species ${species} \\
        --outdir . \\
        --llm_note_r ${llm_note_r}
    """
}
