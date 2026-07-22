/*
 * PGIRL Bioinformatics Analysis Workflow
 *
 * Current scope: bioinformatics only (classification + phylogenetics).
 * Downstream stages (DB query, epi query, evidence integration, intelligence
 * synthesis) are disabled and will be re-enabled when ready.
 *
 * Stages (currently active):
 *   1. Nextclade classification (screen → assign species → full analysis)
 *   2. Nextstrain/ebola phylogenetics (auto-routed for ebolavirus samples)
 *
 * No per-sample splitting is needed — run_nextclade.sh handles multi-sample
 * FASTA natively and creates sample-centric output folders.
 */

include { BIOINFORMATICS } from '../modules/bioinformatics'

workflow SAMPLE_ANALYSIS {
    take:
    input_fasta
    input_metadata
    outdir

    main:
    BIOINFORMATICS(input_fasta, input_metadata, outdir)

    emit:
    nextclade_dir  = BIOINFORMATICS.out.nextclade_dir
    nextstrain_dir = BIOINFORMATICS.out.nextstrain_dir
}
