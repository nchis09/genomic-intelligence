/*
 * Bioinformatics module: Nextclade classification + nextstrain/ebola phylogenetics.
 *
 * Takes the full multi-sample FASTA (no per-sample splitting needed).
 * run_nextclade.sh handles:
 *   1. Screening against all Nextclade datasets
 *   2. Species/pathogen assignment per sample
 *   3. Full Nextclade analysis per sample (sample-centric folders)
 *   4. Auto-routing to nextstrain/ebola for ebolavirus samples
 *
 * Outputs land in:
 *   <outdir>/nextclade_classification/   (Nextclade results)
 *   <outdir>/nextstrain_ebola/           (phylogenetic results)
 */
process BIOINFORMATICS {
    tag "nextclade_classification"

    input:
    path input_fasta
    path input_metadata
    val outdir

    output:
    val "${projectDir}/${outdir}/nextclade_classification", emit: nextclade_dir
    val "${projectDir}/${outdir}/nextstrain_ebola",         emit: nextstrain_dir

    script:
    def runner = "${projectDir}/intelligence_engine/bioinformatics/run_nextclade.sh"
    def out   = "${projectDir}/${outdir}"
    """
    bash "${runner}" "${input_fasta}" "${input_metadata}" "${out}"
    """
}
