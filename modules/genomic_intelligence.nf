process GENOMIC_INTELLIGENCE {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(evidence_dir), val(bio_dir), val(data_query_dir), val(outdir)

    output:
    val "${projectDir}/${outdir}/genomic_intelligence/${sample_id}", emit: genomic_intelligence_dir

    script:
    """
    export PYTHONPATH="${projectDir}"
    python3 -m intelligence_engine.genomic_intelligence.synthesize \
        --evidence-integration-dir "${evidence_dir}" \
        --output-dir "${projectDir}/${outdir}/genomic_intelligence/${sample_id}" \
        --bioinformatics-dir "${bio_dir}" \
        --data-query-dir "${data_query_dir}"
    """
}
