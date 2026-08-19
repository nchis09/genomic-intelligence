/*
 * Subworkflow: REPORTING
 *
 * Final reporting stage. Currently contains only MultiQC.
 *
 * Input:  ch_multiqc_files - channel of MultiQC input files
 *         multiqc_config   - path to MultiQC config file
 *         multiqc_logo     - path to MultiQC logo file
 * Output: multiqc_report
 */

include { MULTIQC } from '../../../modules/nf-core/multiqc/main'

workflow REPORTING {
    take:
    ch_multiqc_files  // channel: MultiQC input files
    multiqc_config
    multiqc_logo

    main:
    def ch_multiqc_config = multiqc_config
        ? file(multiqc_config, checkIfExists: true)
        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_logo = multiqc_logo
        ? file(multiqc_logo, checkIfExists: true)
        : file("${projectDir}/assets/logo.png", checkIfExists: true)

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'genomic-intelligence'],
                files,
                ch_multiqc_config,
                ch_multiqc_logo,
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report  // channel: [ meta, html ]
}
