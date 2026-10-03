/*
 * Local module: GENERATE_SPREAD_ASSESSMENT
 *
 * Pre-generate the per-query "Transmission & Spread Assessment" narrative for
 * one species. Reads the transmission-context TSVs produced by
 * COMPUTE_TRANSMISSION_CONTEXT, assembles a fact list per query, calls the
 * local Ollama model (dashboard/modules/llm_note.R), and writes
 * spread_assessment.json. Falls back to a deterministic template when Ollama
 * is unreachable — the file is always produced.
 */

process GENERATE_SPREAD_ASSESSMENT {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_pathogen_genomics.yml"

    input:
    tuple val(meta), path(tc_dir)
    path llm_note_r
    path module_r

    output:
    tuple val(meta), path("spread_assessment.json"), emit: json

    script:
    def species = meta.species ?: meta.id
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"
    export PG_OLLAMA_HOST="${params.ollama_host}"
    export PG_OLLAMA_MODEL="${params.ollama_model}"

    Rscript ${projectDir}/bin/generate_spread_assessment.R \\
        --tc_dir ${tc_dir} \\
        --species ${species} \\
        --outdir . \\
        --llm_note_r ${llm_note_r} \\
        --module_r ${module_r}
    """
}
