/*
 * Local module: GENERATE_TREE_NOTES
 *
 * Pre-generate the plain-language phylogenetic-tree interpretation note for
 * one species. Reads the exported DuckDB knowledge warehouse, builds the
 * evolutionary summary (dashboard/modules/llm_note.R), calls the local Ollama
 * model, and writes tree_note.json. Falls back to a deterministic template
 * note when Ollama is unreachable — the file is always produced.
 */

process GENERATE_TREE_NOTES {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_pathogen_genomics.yml"

    input:
    tuple val(meta), path(duckdb_file)
    path llm_note_r

    output:
    tuple val(meta), path("tree_note.json"), emit: json

    script:
    def species = meta.species ?: meta.id
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"
    export PG_OLLAMA_HOST="${params.ollama_host}"
    export PG_OLLAMA_MODEL="${params.ollama_model}"

    Rscript ${projectDir}/bin/generate_tree_notes.R \\
        --duckdb ${duckdb_file} \\
        --species ${species} \\
        --outdir . \\
        --llm-note-r ${llm_note_r}
    """
}
