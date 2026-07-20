process BIOINFORMATICS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(sample_fasta), path(sample_metadata)
    val outdir
    val db_url

    output:
    tuple val(sample_id), val("${projectDir}/${outdir}/bioinformatics/${sample_id}"), emit: bio_tuple

    script:
    """
    export PYTHONPATH="${projectDir}"
    python3 -m intelligence_engine.bioinformatics.pipeline \
        --fasta "${sample_fasta}" \
        --metadata "${sample_metadata}" \
        --output-dir "${projectDir}/${outdir}/bioinformatics" \
        --db-url "${db_url}"
    """
}
