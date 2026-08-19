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
    """
    # Auto-download the rhdx tool on first run if the local clone is empty/missing.
    if [ ! -f "${projectDir}/tools/rhdx/DESCRIPTION" ]; then
        rm -rf "${projectDir}/tools/rhdx"
        git clone https://github.com/dickoa/rhdx.git "${projectDir}/tools/rhdx"
    fi

    \$CONDA_PREFIX/bin/Rscript ${projectDir}/bin/fetch_rhdx.R \
        --disease "${search_term}" \
        --rhdx_dir "${projectDir}/tools/rhdx" \
        --mapping "${projectDir}/database/hdx_ebola_datasets.yml" \
        ${species_arg} \
        --outdir . \
        --rows ${rows}
    """
}
