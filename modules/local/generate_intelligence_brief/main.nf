/*
 * Local module: GENERATE_INTELLIGENCE_BRIEF
 *
 * Pre-generate the per-species "Intelligence Brief" — the front-page
 * situation report — for one species. Reads the staged pipeline outputs
 * (species-identification TSVs, transmission-context dir, mutation-profile
 * TSVs, tree_note.json, spread_assessment.json, DuckDB warehouse), assembles
 * a structured fact list, calls the local Ollama model via
 * dashboard/modules/intelligence_brief.R + llm_note.R, and writes
 * intelligence_brief.json. Falls back to a deterministic template when
 * Ollama is unreachable — the file is always produced.
 *
 * Any input may be an empty file list when its upstream stage was skipped
 * (--skip_* params); the R script treats those sections as data gaps.
 */
process GENERATE_INTELLIGENCE_BRIEF {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_pathogen_genomics.yml"

    input:
    tuple val(meta), path(duckdb_file), path(pi_dir), path(tc_dir), path(mp_dir), path(tree_note), path(spread_json)
    path llm_note_r
    path module_r

    output:
    tuple val(meta), path("intelligence_brief.json"), emit: json

    script:
    def species = meta.species ?: meta.id
    // Per-domain literature metadata JSONs live under the published outdir;
    // read by path (not staged) as a fallback when the warehouse's literature
    // tables are empty. Resolved to an absolute path like the knowledge
    // warehouse lookups in main.nf.
    def lit_dir = file("${params.outdir}/literature_retrieval/literature_metadata/${species}").toAbsolutePath()
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"
    export PG_OLLAMA_HOST="${params.ollama_host}"
    export PG_OLLAMA_MODEL="${params.ollama_model}"

    # Optional inputs: only pass flags for staged inputs that actually exist
    # (skipped upstream stages arrive as empty file lists).
    EXTRA_ARGS=""
    if [ -e "${pi_dir}" ];      then EXTRA_ARGS="\$EXTRA_ARGS --pi_dir ${pi_dir}"; fi
    if [ -e "${tc_dir}" ];      then EXTRA_ARGS="\$EXTRA_ARGS --tc_dir ${tc_dir}"; fi
    if [ -e "${mp_dir}" ];      then EXTRA_ARGS="\$EXTRA_ARGS --mp_dir ${mp_dir}"; fi
    if [ -e "${tree_note}" ];   then EXTRA_ARGS="\$EXTRA_ARGS --tree_note ${tree_note}"; fi
    if [ -e "${spread_json}" ]; then EXTRA_ARGS="\$EXTRA_ARGS --spread_json ${spread_json}"; fi
    if [ -e "${duckdb_file}" ]; then EXTRA_ARGS="\$EXTRA_ARGS --duckdb ${duckdb_file}"; fi
    if [ -d "${lit_dir}" ];     then EXTRA_ARGS="\$EXTRA_ARGS --lit_dir ${lit_dir}"; fi

    Rscript ${projectDir}/bin/generate_intelligence_brief.R \\
        --species ${species} \\
        --outdir . \\
        --llm_note_r ${llm_note_r} \\
        --module_r ${module_r} \\
        \$EXTRA_ARGS
    """
}
