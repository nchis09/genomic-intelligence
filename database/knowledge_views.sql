-- Genomic Intelligence Knowledge Warehouse analytical views
-- These views flatten and connect tables for downstream statistical / BI tools.

-- Per-run high-level counts useful for MultiQC-style dashboards and QC checks.
CREATE OR REPLACE VIEW v_run_overview AS
SELECT
    ar.run_id,
    (SELECT COUNT(*) FROM samples s WHERE s.run_id = ar.run_id) AS n_samples,
    (SELECT COUNT(*) FROM genomes g JOIN samples s ON g.sample_id = s.sample_id WHERE s.run_id = ar.run_id) AS n_genomes,
    (SELECT COUNT(DISTINCT m.mutation_id)
     FROM mutations m
     WHERE m.protein_id IN (
        SELECT p.protein_id
        FROM proteins p
        JOIN genes ge ON p.gene_id = ge.gene_id
        JOIN reference_genomes rg ON ge.ref_id = rg.ref_id
        WHERE rg.species IN (SELECT DISTINCT s2.species FROM samples s2 WHERE s2.run_id = ar.run_id)
     )) AS n_mutations,
    (SELECT COUNT(*) FROM phylogenetic_trees pt WHERE pt.run_id = ar.run_id) AS n_trees,
    (SELECT COUNT(*) FROM epidemiological_datasets ed WHERE ed.run_id = ar.run_id) AS n_epi_datasets,
    (SELECT COUNT(*) FROM epidemiological_records er
     JOIN epidemiological_datasets ed ON er.dataset_id = ed.dataset_id
     WHERE ed.run_id = ar.run_id) AS n_epi_records,
    (SELECT COUNT(*) FROM literature_domains ld WHERE ld.run_id = ar.run_id) AS n_literature_domains,
    (SELECT COUNT(*) FROM literature_papers lp WHERE lp.run_id = ar.run_id) AS n_literature_papers,
    (SELECT COUNT(*) FROM literature_extractions le
     WHERE le.paper_id IN (SELECT paper_id FROM literature_papers WHERE run_id = ar.run_id)) AS n_literature_extractions,
    (SELECT COUNT(*) FROM pipeline_outputs po WHERE po.run_id = ar.run_id) AS n_pipeline_outputs
FROM analysis_runs ar;

-- Flatten sample -> mutation -> protein -> phenotype for variant effect analysis.
CREATE OR REPLACE VIEW v_sample_mutation_phenotype AS
SELECT
    s.run_id,
    s.sample_id,
    s.sample_name,
    s.species,
    s.pathogen,
    m.mutation_id,
    m.mutation_label,
    m.position,
    m.mutation_type,
    p.protein_id,
    p.uniprot_accession,
    p.protein_name,
    g.gene_name,
    mp.phenotype_id,
    mp.phenotype,
    mp.effect,
    mp.evidence
FROM samples s
LEFT JOIN sample_mutation sm ON s.sample_id = sm.sample_id
LEFT JOIN mutations m ON sm.mutation_id = m.mutation_id
LEFT JOIN proteins p ON m.protein_id = p.protein_id
LEFT JOIN genes g ON p.gene_id = g.gene_id
LEFT JOIN mutation_phenotypes mp ON m.mutation_id = mp.mutation_id;

-- Sample context with geography and temporal information for outbreak mapping.
CREATE OR REPLACE VIEW v_sample_geo_temporal AS
SELECT
    s.run_id,
    s.sample_id,
    s.sample_name,
    s.species,
    s.pathogen,
    s.collection_date,
    s.country,
    s.admin1,
    s.admin2,
    s.locality,
    s.outbreak,
    s.lineage,
    s.clade,
    gl.location_id,
    gl.latitude,
    gl.longitude,
    gl.location_code,
    tt.tip_id,
    tt.tip_date,
    tt.div,
    tt.aa_mutation_count,
    tt.nuc_mutation_count
FROM samples s
LEFT JOIN sample_geo_location sgl ON s.sample_id = sgl.sample_id
LEFT JOIN geographic_locations gl ON sgl.location_id = gl.location_id
LEFT JOIN tree_tips tt ON s.sample_id = tt.sample_id;

-- Per-species, per-domain literature summary for quick risk-factor overviews.
CREATE OR REPLACE VIEW v_species_literature_summary AS
SELECT
    ld.run_id,
    ld.species,
    ld.domain,
    ld.total_papers,
    ld.clean_count,
    ld.failed_count,
    ld.min_completeness,
    ld.quote_presence_rate,
    ld.duplicate_pmids,
    COUNT(DISTINCT lp.paper_id) AS clean_papers,
    COUNT(DISTINCT le.extraction_id) AS extractions
FROM literature_domains ld
LEFT JOIN literature_papers lp ON lp.domain_id = ld.domain_id AND lp.status = 'clean'
LEFT JOIN literature_extractions le ON le.paper_id = lp.paper_id
GROUP BY ld.run_id, ld.species, ld.domain, ld.total_papers, ld.clean_count,
         ld.failed_count, ld.min_completeness, ld.quote_presence_rate, ld.duplicate_pmids;

-- Long-format extraction table joined back to paper metadata for R/Python import.
CREATE OR REPLACE VIEW v_paper_extraction_flat AS
SELECT
    ld.run_id,
    ld.species,
    ld.domain,
    lp.paper_id,
    lp.pmid,
    lp.pmcid,
    lp.doi,
    lp.year,
    lp.title,
    lp.journal,
    lp.publication_date,
    lp.cited_by_count,
    lp.is_oa,
    lp.qc_score,
    lp.status,
    le.field,
    le.present,
    le.value,
    le.quote,
    le.confidence
FROM literature_extractions le
JOIN literature_papers lp ON le.paper_id = lp.paper_id
JOIN literature_domains ld ON lp.domain_id = ld.domain_id;

-- Link each unknown/query sample to the literature evidence available for its species.
CREATE OR REPLACE VIEW v_sample_literature_risk AS
SELECT
    s.run_id,
    s.sample_id,
    s.sample_name,
    s.species,
    s.pathogen,
    s.country,
    s.collection_date,
    lp.paper_id,
    lp.pmid,
    lp.domain,
    lp.qc_score,
    lp.status,
    le.field,
    le.present,
    le.value,
    le.confidence
FROM samples s
LEFT JOIN literature_papers lp ON s.run_id = lp.run_id AND s.species = lp.species AND lp.status = 'clean'
LEFT JOIN literature_extractions le ON lp.paper_id = le.paper_id
WHERE s.is_query = TRUE OR s.is_query IS NULL;
