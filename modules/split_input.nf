process SPLIT_INPUT {
    tag "split_input"

    input:
    path fasta
    path metadata
    val outdir

    output:
    path "${outdir}/split_input", emit: split_dir

    script:
    """
    export PYTHONPATH="${projectDir}"
    python3 "${projectDir}/scripts/split_samples.py" \
        --fasta "${fasta}" \
        --metadata "${metadata}" \
        --output-dir "${outdir}/split_input"
    """
}
