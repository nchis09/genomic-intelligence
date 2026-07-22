process PUBLISH_REPORTS {
    tag "${genomic_dir}"

    input:
    tuple val(genomic_dir), val(bioinformatics_dir)

    output:
    val("${genomic_dir}")

    script:
    """
    sample_id=\$(basename "${genomic_dir}")
    report_dir="${projectDir}/${params.outdir}/reports/\${sample_id}"
    mkdir -p "\${report_dir}"
    ${params.python} "${projectDir}/scripts/publish_reports.py" \
        --genomic-dir "${genomic_dir}" \
        --bioinformatics-dir "${bioinformatics_dir}" \
        --output-dir "\${report_dir}"
    """
}
