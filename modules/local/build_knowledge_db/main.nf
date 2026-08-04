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
    tuple val(meta), path(species_assignments), path(metadata_tsv), path(epi_raw_dir), path(epi_search_summary), path(bioinformatics_results), path(uniprotr_results), path(extractr_results), path(rbioapi_results), path(auspice_json), path(query_data_files, stageAs: 'query_data/*'), path(hmm_files, stageAs: 'hmm/*')
    val(db_host)
    val(db_port)

    output:
    tuple val(meta), path("knowledge_warehouse"), emit: knowledge_db
    tuple val(meta), path("knowledge_warehouse/*_mqc.tsv"), emit: mqc_summary

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = ["--outdir knowledge_warehouse", "--meta-id ${meta.id}", "--prefix ${prefix}"]
    if (meta.species) args << "--species ${meta.species}"
    if (!species_assignments.name.startsWith('NO_FILE')) args << "--species-assignments ${species_assignments}"
    if (!metadata_tsv.name.startsWith('NO_FILE'))          args << "--metadata-tsv ${metadata_tsv}"
    if (!epi_raw_dir.name.startsWith('NO_FILE'))           args << "--epi-raw-dir ${epi_raw_dir}"
    if (!epi_search_summary.name.startsWith('NO_FILE'))    args << "--epi-search-summary ${epi_search_summary}"
    if (!uniprotr_results.name.startsWith('NO_FILE'))      args << "--uniprotr-dir ${uniprotr_results}"
    if (!extractr_results.name.startsWith('NO_FILE'))      args << "--extractr-dir ${extractr_results}"
    if (!rbioapi_results.name.startsWith('NO_FILE'))       args << "--rbioapi-dir ${rbioapi_results}"
    if (!auspice_json.name.startsWith('NO_FILE'))          args << "--auspice-json ${auspice_json}"
    def qd_list  = query_data_files instanceof List ? query_data_files : [query_data_files]
    def qd_real  = qd_list.findAll { !it.name.startsWith('NO_FILE') }
    if (qd_real)                                           args << "--query-data-dir query_data"
    def hmm_list = hmm_files instanceof List ? hmm_files : [hmm_files]
    def hmm_real = hmm_list.findAll { !it.name.startsWith('NO_FILE') }
    if (hmm_real)                                          args << "--hmm-dir hmm"
    // `bioinformatics_results` is the guaranteed-staged copy of NEXTSTRAIN_EBOLA's
    // own "results/" output dir (already an input above) -- reuse it directly
    // instead of params.outdir, which is a relative string that would otherwise
    // resolve inside this task's own work directory, not the real output tree.
    args << "--results-dir ${bioinformatics_results}"
    // Connect to the pipeline-lifetime shared Postgres (started by
    // START_KNOWLEDGE_DB before any species processing begins) instead of
    // spinning up a throwaway per-species instance -- lets downstream figures
    // stages query this same live database after ingestion completes.
    args << "--db-host ${db_host} --db-port ${db_port}"
    """
    mkdir -p knowledge_warehouse

    # Ensure conda env's Python is used (old-style source activate may not update PATH)
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"

    python3 ${projectDir}/bin/build_knowledge_db.py ${args.join(' ')}
    """
}
