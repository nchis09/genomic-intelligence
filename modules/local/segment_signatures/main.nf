/*
 * Local module: SEGMENT_SIGNATURES
 *
 * Influenza-only downstream signature tables (bin/flu_segment_signatures.py),
 * run once per flu species group after SEGMENT_MUTATIONS:
 *
 *   constellation.tsv      isolate x segment clade/genotype labels, the
 *                          PB2|PB1|PA|HA|NP|NA|MP|NS signature string and a
 *                          novelty flag vs the group modal signature.
 *   reassortment_flags.tsv isolates flagged by species_assignments or by
 *                          per-segment label/record-species discordance.
 *   di_candidates.tsv      consensus-level defective-interfering screen
 *                          from each record's winning Nextclade QC row.
 *   markers.tsv            query-isolate mutations joined against the
 *                          curated database/flu_markers.yml catalogue
 *                          (resistance / virulence / host adaptation).
 *   group_diversity.tsv    shared/private mutation counts per group_id x
 *                          segment (intra-host / same-site proxy).
 *
 * Inputs beyond the group's own outputs are broadcast: species_assignments
 * and every classification Nextclade TSV (per-record QC + clade) — the
 * script filters them to this group's pathogen/species and prefers each
 * record's best_dataset_file. `genoflu` is a real TSV for avian groups and
 * a NO_FILE placeholder for seasonal ones.
 */

process SEGMENT_SIGNATURES {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"
    container null

    input:
    tuple val(meta), path(mutations_tsvs), path(build_dir), path(metadata), path(genoflu)
    path(assignments)
    path(nextclade_tsvs, stageAs: 'nc_inputs/*')
    path(markers)

    output:
    tuple val(meta), path("*.tsv")                   , emit: tables
    tuple val(meta), path("signatures_report.txt")   , emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def tsvs    = mutations_tsvs instanceof List ? mutations_tsvs : [ mutations_tsvs ]
    def mut_args = tsvs
        .findAll { it && !it.name.startsWith('NO_FILE') }
        .collect { it.toString() }
        .join(' ')
    def genoflu_arg = (genoflu && !genoflu.name.startsWith('NO_FILE')) ? "--genoflu ${genoflu}" : ''
    def metadata_arg = (metadata && !metadata.name.startsWith('NO_FILE')) ? "--metadata ${metadata}" : ''
    """
    python3 ${projectDir}/bin/flu_segment_signatures.py \\
        ${mut_args ? "--mutations ${mut_args}" : ''} \\
        --assignments "${assignments}" \\
        --nextclade nc_inputs/*.tsv \\
        ${metadata_arg} \\
        ${genoflu_arg} \\
        --build-dir "${build_dir}" \\
        --markers "${markers}" \\
        --species "${meta.species ?: ''}" \\
        --pathogen "${meta.pathogen ?: ''}" \\
        --outdir . \\
        --report signatures_report.txt
    """
}
