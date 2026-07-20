/*
 * End-to-end per-sample analysis workflow.
 *
 * Each sample is processed as a separate Nextflow task in every stage, so the
 * executor progress shows e.g. "3 of 5 done" rather than waiting for a whole
 * batch to finish.
 *
 * Stages:
 *   1. Split multi-sample FASTA + metadata into per-sample inputs
 *   2. Bioinformatics per sample
 *   3. Local DB query per sample
 *   4. Online epi query per sample
 *   5. Evidence integration per sample
 *   6. Genomic intelligence synthesis per sample
 */

include { SPLIT_INPUT } from '../modules/split_input'
include { BIOINFORMATICS } from '../modules/bioinformatics'
include { DB_QUERY } from '../modules/db_query'
include { EPI_QUERY } from '../modules/epi_query'
include { EVIDENCE_INTEGRATION } from '../modules/evidence_integration'
include { GENOMIC_INTELLIGENCE } from '../modules/genomic_intelligence'

workflow SAMPLE_ANALYSIS {
    take:
    input_fasta
    input_metadata
    outdir
    db_url
    use_llm

    main:
    split_dir = SPLIT_INPUT(input_fasta, input_metadata, outdir)

    // Build one channel item per sample: (sample_id, sample_fasta, sample_metadata)
    sample_ch = split_dir
        | flatMap { dir ->
            dir.listFiles()
                .findAll { it.isDirectory() }
                .collect { sample_dir ->
                    tuple(
                        sample_dir.name,
                        sample_dir / "sample.fasta",
                        sample_dir / "metadata.csv",
                    )
                }
        }

    // Per-sample bioinformatics, DB query, and epi query.
    bio_ch  = BIOINFORMATICS(sample_ch, outdir, db_url)
    db_ch   = DB_QUERY(bio_ch.map { sample_id, bio_dir ->
        tuple(sample_id, new File(bio_dir).getParent())
    }, outdir, db_url)
    epi_ch  = EPI_QUERY(bio_ch.combine(db_ch, by: 0), outdir, use_llm)

    // Evidence integration needs (sample_id, bio_dir, data_query_dir, outdir, db_url)
    evidence_ch = epi_ch
        .combine(bio_ch, by: 0)
        | map { sample_id, data_query_dir, bio_dir ->
            tuple(sample_id, bio_dir, data_query_dir, outdir, db_url)
        }
        | EVIDENCE_INTEGRATION

    // Genomic intelligence needs (sample_id, evidence_dir, bio_dir, data_query_dir, outdir)
    genomic_ch = evidence_ch
        .combine(bio_ch, by: 0)
        .combine(db_ch, by: 0)
        | map { sample_id, evidence_dir, bio_dir, data_query_dir ->
            tuple(sample_id, evidence_dir, bio_dir, data_query_dir, outdir)
        }
        | GENOMIC_INTELLIGENCE

    emit:
    genomic_intelligence_dir = genomic_ch
}
