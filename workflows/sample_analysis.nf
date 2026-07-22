/*
 * PGIRL Bioinformatics Analysis Workflow
 *
 * Current scope: bioinformatics only (classification + phylogenetics).
 * Downstream stages (DB query, epi query, evidence integration, intelligence
 * synthesis) are disabled and will be re-enabled when ready.
 *
 * Stages (currently active — each shown as separate Nextflow progress):
 *   1. NEXTCLADE_SCREEN   — screen sequences against 15 pathogen datasets
 *   2. NEXTCLADE_ASSIGN   — assign species per sample based on best QC
 *   3. NEXTCLADE_ANALYZE  — full Nextclade analysis per sample
 *   4. NEXTSTRAIN_PHYLO   — nextstrain/ebola phylogenetic tree building
 *
 * No per-sample splitting is needed — run_nextclade.sh handles multi-sample
 * FASTA natively and creates sample-centric output folders.
 */

include { NEXTCLADE_SCREEN  } from '../modules/bioinformatics'
include { NEXTCLADE_ASSIGN  } from '../modules/bioinformatics'
include { NEXTCLADE_ANALYZE } from '../modules/bioinformatics'
include { NEXTSTRAIN_CONFIG } from '../modules/bioinformatics'
include { NEXTSTRAIN_BUILD  } from '../modules/bioinformatics'

workflow SAMPLE_ANALYSIS {
    take:
    input_fasta
    input_metadata
    outdir

    main:
    // Step 1: Screen against all datasets
    screen_ch = NEXTCLADE_SCREEN(input_fasta, input_metadata, outdir)

    // Step 2: Assign species per sample
    assign_ch = NEXTCLADE_ASSIGN(input_fasta, input_metadata, screen_ch.out_dir)

    // Read assignments.tsv to get sample IDs, then fan out per sample
    sample_ids = assign_ch.out_dir
        .map { out_dir ->
            def tsv = new File("${out_dir}/nextclade_classification/assignments.tsv")
            tsv.readLines().drop(1)  // skip header
                .collect { line -> line.split('\t')[0] }  // seqName column
        }
        .flatMap { it }

    // Step 3: Full Nextclade analysis — one task per sample
    analyze_ch = NEXTCLADE_ANALYZE(
        input_fasta,
        input_metadata,
        assign_ch.out_dir,
        sample_ids,
    )

    // Step 4a: Generate nextstrain/ebola config (after all samples done)
    config_ch = NEXTSTRAIN_CONFIG(
        input_fasta,
        input_metadata,
        analyze_ch.out_dir.collect().map { it[0] },
    )

    // Read ebola species from assignments.tsv, fan out per species
    ebola_species = config_ch.out_dir
        .map { out_dir ->
            def tsv = new File("${out_dir}/nextclade_classification/assignments.tsv")
            tsv.readLines().drop(1)
                .findAll { line -> line.split('\t')[3] == 'orthoebolavirus' }
                .collect { line -> line.split('\t')[2] }  // species_label column
                .unique()
        }
        .flatMap { it }

    // Step 4b: Build nextstrain/ebola — one task per species
    build_ch = NEXTSTRAIN_BUILD(
        input_fasta,
        input_metadata,
        config_ch.out_dir,
        ebola_species,
    )

    emit:
    nextstrain_dir = build_ch.out_dir
}
