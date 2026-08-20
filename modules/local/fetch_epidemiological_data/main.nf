/*
 * Local module: FETCH_EPIDEMIOLOGICAL_DATA
 *
 * Search HDX for a disease term using the rhdx R package and download
 * the first CSV/XLSX/XLS/JSON resource found for the species group.
 */

process FETCH_EPIDEMIOLOGICAL_DATA {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_rhdx.yml"
    container null

    input:
    tuple val(meta), val(search_term), val(species)

    output:
    tuple val(meta), path("epi_data"),               emit: epi_raw
    tuple val(meta), path("rhdx_search_results.tsv"), emit: search_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def rows = params.epi_hdx_rows ?: 20
    def species_arg = species ? "--species \"${species}\"" : ""
    def hdx_key_arg = params.hdx_key ? "--hdx-key \"${params.hdx_key}\"" : ""
    def user_agent_arg = params.hdx_user_agent ? "--user-agent \"${params.hdx_user_agent}\"" : ""
    def ignore_hdx_errors_arg = params.ignore_hdx_errors ? "--ignore-hdx-errors" : ""
    """
    # Optional HTTP(S) proxy for compute nodes without direct outbound access.
    [ -n "${params.http_proxy}" ] && [ "${params.http_proxy}" != "null" ] && export http_proxy="${params.http_proxy}"
    [ -n "${params.https_proxy}" ] && [ "${params.https_proxy}" != "null" ] && export https_proxy="${params.https_proxy}"

    # Ensure rhdx is cloned and installed once in a concurrency-safe way.
    RHDX_DIR="\$(${projectDir}/bin/setup_tool.sh rhdx https://github.com/dickoa/rhdx.git --install)"

    # Point R at the shared tool library where rhdx was installed.
    export R_LIBS_USER="\${RHDX_DIR}/r_library:\${R_LIBS_USER:-}"

    \$CONDA_PREFIX/bin/Rscript ${projectDir}/bin/fetch_rhdx.R \
        --disease "${search_term}" \
        --rhdx_dir "${projectDir}/tools/rhdx" \
        --mapping "${projectDir}/database/hdx_ebola_datasets.yml" \
        ${species_arg} \
        --outdir . \
        --rows ${rows} \
        ${hdx_key_arg} \
        ${user_agent_arg} \
        ${ignore_hdx_errors_arg}
    """
}
