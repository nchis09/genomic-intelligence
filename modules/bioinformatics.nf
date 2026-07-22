/*
 * Bioinformatics module: Nextclade classification + nextstrain/ebola phylogenetics.
 *
 * Split into 4 stages so Nextflow shows per-step progress:
 *   1. NEXTCLADE_SCREEN   — screen all sequences against 15 datasets (1 task)
 *   2. NEXTCLADE_ASSIGN   — assign species per sample (1 task)
 *   3. NEXTCLADE_ANALYZE  — full Nextclade analysis (1 task per sample)
 *   4. NEXTSTRAIN_PHYLO   — nextstrain/ebola phylogenetics (1 task)
 *
 * Outputs land in:
 *   <outdir>/nextclade_classification/   (Nextclade results)
 *   <outdir>/nextstrain_ebola/           (phylogenetic results)
 */

def runner = "${projectDir}/intelligence_engine/bioinformatics/run_nextclade.sh"

process NEXTCLADE_SCREEN {
    tag "screen_datasets"

    input:
    path input_fasta
    path input_metadata
    val outdir

    output:
    val "${projectDir}/${outdir}", emit: out_dir

    script:
    def out = "${projectDir}/${outdir}"
    """
    bash "${runner}" "${input_fasta}" "${input_metadata}" "${out}" --step screen
    """
}

process NEXTCLADE_ASSIGN {
    tag "assign_species"

    input:
    path input_fasta
    path input_metadata
    val out_dir

    output:
    val out_dir, emit: out_dir

    script:
    """
    bash "${runner}" "${input_fasta}" "${input_metadata}" "${out_dir}" --step assign
    """
}

process NEXTCLADE_ANALYZE {
    tag "${sample_id}"

    input:
    path input_fasta
    path input_metadata
    val out_dir
    val sample_id

    output:
    val out_dir, emit: out_dir

    script:
    """
    bash "${runner}" "${input_fasta}" "${input_metadata}" "${out_dir}" --step analyze --sample "${sample_id}"
    """
}

process NEXTSTRAIN_CONFIG {
    tag "phylo_config"

    input:
    path input_fasta
    path input_metadata
    val out_dir

    output:
    val out_dir, emit: out_dir

    script:
    """
    bash "${runner}" "${input_fasta}" "${input_metadata}" "${out_dir}" --step phylo_config
    """
}

process NEXTSTRAIN_BUILD {
    tag "${species}"
    maxForks 1
    errorStrategy 'ignore'

    input:
    path input_fasta
    path input_metadata
    val out_dir
    val species

    output:
    val out_dir, emit: out_dir

    script:
    """
    bash "${runner}" "${input_fasta}" "${input_metadata}" "${out_dir}" --step phylo_build --species "${species}"
    """
}
