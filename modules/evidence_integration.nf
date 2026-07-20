process EVIDENCE_INTEGRATION {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(bio_dir), val(data_query_dir), val(outdir), val(db_url)

    output:
    tuple val(sample_id), val("${projectDir}/${outdir}/evidence_integration/${sample_id}"), emit: evidence_tuple

    script:
    """
    export PYTHONPATH="${projectDir}"
    python3 -m intelligence_engine.evidence_integration.pipeline.intelligence_pipeline \
        --bio-output "${bio_dir}/bio_output.json" \
        --epi-output "${data_query_dir}/epi_output.json" \
        --tree-file "${bio_dir}/tree.nwk" \
        --output-dir "${projectDir}/${outdir}/evidence_integration/${sample_id}" \
        --associations "database/exports/genotype_phenotype.csv" \
        --variants "database/exports/protein_variants.csv" \
        --lineages "database/exports/lineages.csv" \
        --genome-metadata "database/exports/genome_metadata.csv" \
        --db-url "${db_url}"
    """
}
