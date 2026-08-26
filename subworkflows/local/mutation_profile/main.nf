/*
 * Subworkflow: MUTATION_PROFILE_WF
 *
 * Runs the overall mutation profile analysis (analysis #1) per detected
 * species from the exported knowledge warehouse DuckDB.
 *
 * Input:  ch_mutation_profile - channel: [ val(meta), path(duckdb_file) ]
 *         ch_translations     - channel: path(translations_dir)
 * Output: tsv     - channel: path(*.tsv)
 *         mqc_tsv - channel: path(*_mqc.tsv)
 */

include { PATHOGEN_MUTATION_PROFILE } from '../../../modules/local/pathogen_mutation_profile/main'

workflow MUTATION_PROFILE_WF {
    take:
    ch_mutation_profile // channel: [ meta, duckdb_file ]
    ch_translations     // channel: path(translations_dir)

    main:
    PATHOGEN_MUTATION_PROFILE(ch_mutation_profile, ch_translations)

    emit:
    tsv     = PATHOGEN_MUTATION_PROFILE.out.tsv
    mqc_tsv = PATHOGEN_MUTATION_PROFILE.out.mqc_tsv
}
