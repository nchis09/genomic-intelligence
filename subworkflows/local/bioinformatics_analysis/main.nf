/*
 * Subworkflow: BIOINFORMATICS_ANALYSIS
 *
 * Pathogen-specific bioinformatics stage. Currently wraps the Nextstrain
 * Ebola Snakemake workflow. Additional tools can be added here in future
 * pathogen workflows.
 *
 * Input:  ch_species_data - channel of [ meta, fasta, metadata ]
 * Output: auspice, results
 */

include { NEXTSTRAIN_EBOLA             } from '../../../modules/local/nextstrain_ebola/main'
include { NEXTSTRAIN_EBOLA_INGEST      } from '../../../modules/local/nextstrain_ebola_ingest/main'
include { MAFFT_ALIGN                  } from '../../../modules/local/mafft_align/main'
include { IQTREE2                      } from '../../../modules/local/iqtree2/main'

workflow BIOINFORMATICS_ANALYSIS {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    //
    // Fetch/update the Nextstrain Ebola background sequences and metadata.
    // Skipped when --skip_nextstrain_ingest is true, in which case the
    // pipeline expects pre-generated files in data/nextstrain_ebola/ingest/data/.
    //
    ch_ingest_bg = params.skip_nextstrain_ingest
        ? ch_species_data.map { meta, fasta, metadata ->
            def species = (meta.species ?: meta.id).replaceAll(/^.*_/, '')
            [ meta,
              file("${projectDir}/data/nextstrain_ebola/ingest/data/${species}/sequences.fasta"),
              file("${projectDir}/data/nextstrain_ebola/ingest/data/${species}/metadata.tsv") ]
        }
        : NEXTSTRAIN_EBOLA_INGEST(ch_species_data.map { it[0] }).background

    ch_ebola_input = ch_species_data.join(ch_ingest_bg)

    //
    // Real bioinformatics workhorse: Nextstrain Ebola Snakemake workflow
    //
    NEXTSTRAIN_EBOLA(ch_ebola_input)

    //
    // Build a model-aware ML tree from the Nextstrain subsampled sequences
    //
    MAFFT_ALIGN(NEXTSTRAIN_EBOLA.out.results_dir)
    IQTREE2(MAFFT_ALIGN.out.alignment)

    emit:
    auspice   = NEXTSTRAIN_EBOLA.out.auspice      // channel: [ meta, json ]
    results   = NEXTSTRAIN_EBOLA.out.results_dir  // channel: [ meta, dir ]
    alignment = MAFFT_ALIGN.out.alignment         // channel: [ meta, fasta ]
    tree      = IQTREE2.out.tree                  // channel: [ meta, newick ]
}
