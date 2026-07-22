process RUN_SUMMARY {
    tag "run_summary"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val _trigger  // signal that all per-sample stages have finished

    output:
    path "run_summary.json"

    script:
    """
    ${params.python} "${projectDir}/scripts/run_summary.py" \
        --outdir "${projectDir}/${params.outdir}" \
        --output run_summary.json
    """
}
