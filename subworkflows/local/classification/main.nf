/*
 * Subworkflow: CLASSIFICATION
 *
 * Run Nextclade classification on input sequences:
 *   1. Download Nextclade dataset(s)
 *   2. Run Nextclade on each sample against each dataset
 *   3. Assign species based on best Nextclade QC score
 *   4. Group samples by species and emit per-species FASTA + metadata
 *
 * Input:  ch_samplesheet - channel of [ meta, fasta, metadata ]
 *         ch_datasets    - channel of dataset names to screen against
 * Output: species_groups - channel of [ meta(pathogen, species), fasta, metadata ] per species
 */

include { NEXTCLADE_DATASETGET    } from '../../../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN           } from '../../../modules/nf-core/nextclade/run/main'
include { SPECIES_ASSIGN          } from '../../../modules/local/species_assign/main'

workflow CLASSIFICATION {
    take:
    ch_samplesheet // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_datasets    // channel: val(dataset_name)

    main:
    //
    // MODULE: Download Nextclade datasets
    //
    NEXTCLADE_DATASETGET(
        ch_datasets,
        channel.value([])  // no specific tag → latest
    )

    //
    // MODULE: Run Nextclade classification
    //
    // Build fasta-only channel for Nextclade
    ch_fasta = ch_samplesheet.map { meta, fasta, _metadata -> [ meta, fasta ] }

    // Combine each sample FASTA with each downloaded dataset
    // so Nextclade runs once per sample × dataset combination
    // Add dataset name to meta.id to avoid filename collisions
    ch_nextclade_input = ch_fasta
        .combine(NEXTCLADE_DATASETGET.out.dataset)
        .map { meta, fasta, dataset ->
            def dataset_suffix = dataset.name  // directory name (e.g., "bdbv", "sudan")
            def new_meta = meta + [id: "${meta.id}_${dataset_suffix}"]
            [ new_meta, fasta, dataset ]
        }

    NEXTCLADE_RUN(
        ch_nextclade_input.map { meta, fasta, dataset -> [ meta, fasta ] },
        ch_nextclade_input.map { meta, fasta, dataset -> dataset }
    )

    //
    // MODULE: Assign species based on Nextclade QC scores
    //
    ch_all_tsvs = NEXTCLADE_RUN.out.tsv
        .map { meta, tsv -> tsv }
        .collect()

    SPECIES_ASSIGN(
        ch_samplesheet.first(),
        ch_all_tsvs
    )

    //
    // Build per-species channel from SPECIES_ASSIGN outputs
    // Each species gets its own [meta, fasta, metadata] tuple
    //
    ch_species_groups = SPECIES_ASSIGN.out.species_groups_json
        .splitJson()
        .map { group ->
            def meta = [
                id: "${group.pathogen}_${group.species}",
                pathogen: group.pathogen,
                species: group.species,
                query_samples: group.samples.join(',')
            ]
            def fasta = file(group.fasta)
            def metadata = file(group.metadata)
            [ meta, fasta, metadata ]
        }

    emit:
    species_groups = ch_species_groups             // channel: [ meta(pathogen, species), fasta, metadata ]
    tsv            = NEXTCLADE_RUN.out.tsv         // channel: [ meta, tsv ]
    json           = NEXTCLADE_RUN.out.json        // channel: [ meta, json ]
    assignments    = SPECIES_ASSIGN.out.assignments // path: species_assignments.tsv
}
