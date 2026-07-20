-- ============================================================================
-- Database Optimizations: Indexes, GIN, Full-Text Search, Materialized Views
-- Run: psql pgirl -f database/optimizations.sql
-- ============================================================================

-- ============================================================================
-- 1. FOREIGN KEY INDEXES (missing ones only)
--    PostgreSQL does NOT auto-index FK columns; JOINs need these.
-- ============================================================================

-- genotype_phenotype FK indexes (mutation_id already indexed, add the rest)
CREATE INDEX IF NOT EXISTS idx_gp_lineage_id      ON genotype_phenotype(lineage_id);
CREATE INDEX IF NOT EXISTS idx_gp_species_id      ON genotype_phenotype(species_id);
CREATE INDEX IF NOT EXISTS idx_gp_pathogen_id     ON genotype_phenotype(pathogen_id);

-- genome_metadata FK indexes (most already exist, ensure lineage_id)
CREATE INDEX IF NOT EXISTS idx_genomemeta_lineage_id ON genome_metadata(lineage_id);

-- gene_function FK indexes (species_id and pathogen_id already indexed)

-- outbreaks FK indexes (lineage_id missing)
CREATE INDEX IF NOT EXISTS idx_outbreaks_lineage_id ON outbreaks(lineage_id);

-- lineages FK indexes
CREATE INDEX IF NOT EXISTS idx_lineages_pathogen_id ON lineages(pathogen_id);
CREATE INDEX IF NOT EXISTS idx_lineages_species_id  ON lineages(species_id);

-- disease_epidemiology FK indexes
CREATE INDEX IF NOT EXISTS idx_disepi_pathogen_id ON disease_epidemiology(pathogen_id);
CREATE INDEX IF NOT EXISTS idx_disepi_species_id  ON disease_epidemiology(species_id);

-- api_refresh_log FK index
CREATE INDEX IF NOT EXISTS idx_apirefresh_pathogen_id ON api_refresh_log(pathogen_id);

-- protein_variants FK indexes (genome_accession already indexed)
-- reference_accession already indexed, species_id already indexed, pathogen_id already indexed

-- ============================================================================
-- 2. COMPOSITE INDEXES (hot query paths)
-- ============================================================================

-- mutations: (species_id, protein, position) — the #1 intelligence lookup
CREATE INDEX IF NOT EXISTS idx_mutations_scope_position
    ON mutations(species_id, protein, "position");

-- mutations: (pathogen_id, species_id, protein, position) — full scope
CREATE INDEX IF NOT EXISTS idx_mutations_full_scope
    ON mutations(pathogen_id, species_id, protein, "position");

-- genotype_phenotype: (species_id, protein, position) — phenotype by location
CREATE INDEX IF NOT EXISTS idx_gp_scope_position
    ON genotype_phenotype(species_id, protein, "position");

-- genotype_phenotype: (pathogen_id, species_id, protein, position)
CREATE INDEX IF NOT EXISTS idx_gp_full_scope
    ON genotype_phenotype(pathogen_id, species_id, protein, "position");

-- protein_variants: (species_id, gene, position, alt_aa) — variant lookup
CREATE INDEX IF NOT EXISTS idx_proteinvars_variant_lookup
    ON protein_variants(species_id, gene, "position", alt_aa);

-- genome_metadata: (species_id, collection_year, collection_country)
CREATE INDEX IF NOT EXISTS idx_genomemeta_temporal_geo
    ON genome_metadata(species_id, collection_year, collection_country);

-- gene_function: (pathogen_id, species_id, gene) — function lookup
CREATE INDEX IF NOT EXISTS idx_genefunc_full_scope
    ON gene_function(pathogen_id, species_id, gene);

-- ============================================================================
-- 3. GIN INDEXES ON JSONB COLUMNS
--    Enables @> containment queries (e.g., "domain at position 50")
-- ============================================================================

CREATE INDEX IF NOT EXISTS idx_genefunc_domains_gin
    ON gene_function USING GIN (key_domains);

CREATE INDEX IF NOT EXISTS idx_genefunc_sites_gin
    ON gene_function USING GIN (functional_sites);

CREATE INDEX IF NOT EXISTS idx_genefunc_hotspots_gin
    ON gene_function USING GIN (known_hotspots);

CREATE INDEX IF NOT EXISTS idx_genefunc_conserved_gin
    ON gene_function USING GIN (conserved_regions);

-- ============================================================================
-- 4. GIN INDEXES ON ARRAY COLUMNS
--    Enables && overlap and @> containment queries
-- ============================================================================

CREATE INDEX IF NOT EXISTS idx_genefunc_pdb_ids_gin
    ON gene_function USING GIN (pdb_ids);

CREATE INDEX IF NOT EXISTS idx_genefunc_literature_gin
    ON gene_function USING GIN (literature_refs);

CREATE INDEX IF NOT EXISTS idx_gp_literature_gin
    ON genotype_phenotype USING GIN (literature_refs);

CREATE INDEX IF NOT EXISTS idx_mutations_outbreaks_gin
    ON mutations USING GIN (reported_in_outbreaks);

CREATE INDEX IF NOT EXISTS idx_mutations_genbank_gin
    ON mutations USING GIN (genbank_accessions);

-- ============================================================================
-- 5. PARTIAL INDEXES (surveillance dashboards)
-- ============================================================================

-- Only public-health-relevant mutations (small, fast subset)
CREATE INDEX IF NOT EXISTS idx_mutations_ph_relevant
    ON mutations(species_id, protein, "position")
    WHERE public_health_relevant = true;

-- Only flagged genotype-phenotype records
CREATE INDEX IF NOT EXISTS idx_gp_flagged
    ON genotype_phenotype(species_id, protein, "position")
    WHERE record_flagged = true;

-- Only verified mutations
CREATE INDEX IF NOT EXISTS idx_mutations_verified
    ON mutations(species_id, protein, "position")
    WHERE verification_status = 'verified';

-- Only high-evidence genotype-phenotype
CREATE INDEX IF NOT EXISTS idx_gp_high_evidence
    ON genotype_phenotype(species_id, protein, "position")
    WHERE evidence_strength IN ('strong', 'moderate');

-- ============================================================================
-- 6. FULL-TEXT SEARCH ON protein_function
--    Enables: SELECT ... FROM gene_function
--             WHERE protein_function_tsv @@ plainto_tsquery('receptor binding')
-- ============================================================================

-- Add a generated tsvector column (no triggers needed)
ALTER TABLE gene_function
    ADD COLUMN IF NOT EXISTS protein_function_tsv tsvector
    GENERATED ALWAYS AS (to_tsvector('english', coalesce(protein_function, ''))) STORED;

CREATE INDEX IF NOT EXISTS idx_genefunc_fts
    ON gene_function USING GIN (protein_function_tsv);

-- ============================================================================
-- 7. MATERIALIALIZED VIEW: v_mutation_intelligence (snapshot)
--    Replaces the regular view for dashboard / repeated queries.
--    Refresh with: REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_intelligence;
-- ============================================================================

DROP MATERIALIZED VIEW IF EXISTS mv_mutation_intelligence CASCADE;

CREATE MATERIALIZED VIEW mv_mutation_intelligence AS
SELECT
    m.mutation_id,
    m.pathogen_id,
    m.species_id,
    m.hgvs_protein,
    m.protein,
    m."position",
    m.ref_aa,
    m.alt_aa,
    m.frequency,
    m.first_reported_year,
    m.reported_in_outbreaks,
    m.public_health_relevant,
    gp.association_id,
    gp.phenotype_category,
    gp.phenotype_specific,
    gp.effect_size,
    gp.study_type,
    gp.model_system,
    gp.evidence_strength,
    gp.literature_refs   AS gp_literature_refs,
    gp.record_flagged,
    gf.protein_name,
    gf.protein_function,
    gf.key_domains,
    gf.functional_sites,
    gf.pdb_ids           AS gf_pdb_ids,
    gf.literature_refs   AS gf_literature_refs
FROM mutations m
LEFT JOIN genotype_phenotype gp
    ON m.mutation_id = gp.mutation_id
LEFT JOIN gene_function gf
    ON m.species_id = gf.species_id
    AND m.protein = gf.gene
WITH DATA;

CREATE UNIQUE INDEX idx_mv_mut_intelligence_pk
    ON mv_mutation_intelligence(mutation_id, species_id, protein, "position");

CREATE INDEX idx_mv_mut_intelligence_scope
    ON mv_mutation_intelligence(species_id, protein, "position");

CREATE INDEX idx_mv_mut_intelligence_pathogen
    ON mv_mutation_intelligence(pathogen_id);

CREATE INDEX idx_mv_mut_intelligence_ph
    ON mv_mutation_intelligence(species_id, protein, "position")
    WHERE public_health_relevant = true;

-- ============================================================================
-- 8. CROSS-TABLE LINKING VIEW: Genomic Intelligence Snapshot
--    genome_metadata → protein_variants → gene_function → genotype_phenotype
--    One row per variant observation with full biological context.
-- ============================================================================

CREATE OR REPLACE VIEW v_genomic_intelligence_snapshot AS
SELECT
    pv.variant_id,
    pv.pathogen_id,
    pv.species_id,
    pv.gene,
    pv."position",
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    pv.is_synonymous,
    pv.is_stop,
    pv.is_frameshift,
    pv.genome_accession,
    gm.strain,
    gm.isolate,
    gm.collection_date,
    gm.collection_year,
    gm.collection_country,
    gm.collection_country_code,
    gm.host,
    gm.lineage_id,
    gf.protein_name,
    gf.protein_function,
    gf.protein_length_aa,
    gf.key_domains,
    gf.functional_sites,
    gf.pdb_ids,
    gp.association_id,
    gp.phenotype_category,
    gp.phenotype_specific,
    gp.effect_size,
    gp.evidence_strength,
    gp.study_type,
    gp.record_flagged
FROM protein_variants pv
JOIN genome_metadata gm
    ON pv.genome_accession = gm.genome_accession
LEFT JOIN gene_function gf
    ON pv.species_id = gf.species_id
    AND pv.gene = gf.gene
LEFT JOIN genotype_phenotype gp
    ON pv.species_id = gp.species_id
    AND pv.gene = gp.protein
    AND pv."position" = gp."position"
    AND pv.ref_aa = gp.ref_aa
    AND pv.alt_aa = gp.alt_aa;

-- ============================================================================
-- 9. REFRESH FUNCTION (concurrent refresh for mat views)
-- ============================================================================

CREATE OR REPLACE FUNCTION refresh_intelligence_matviews()
RETURNS void AS $$
BEGIN
    -- Refresh existing surveillance mat views
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_with_phenotype;
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_frequency;
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_phenotype_surveillance;
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_phenotype_geo_temporal;
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_position_variants;
    -- New intelligence mat view
    REFRESH MATERIALIZED VIEW CONCURRENTLY mv_mutation_intelligence;
EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Some mat views could not be refreshed concurrently: %', SQLERRM;
END;
$$ LANGUAGE plpgsql;

-- ============================================================================
-- DONE
-- ============================================================================
