-- =============================================================================
-- PGIRL Aggregation Migration
-- =============================================================================
-- Replaces per-genome variant storage with pre-aggregated variant counts.
-- This reduces Ebola from 300K rows / 1.3 GB to ~7K rows / ~15 MB.
-- Scales linearly with distinct variants, not quadratically with genomes.
-- =============================================================================

BEGIN;

-- =============================================================================
-- 1. Drop ALL views and materialized views that depend on protein_variants
-- =============================================================================

DROP VIEW IF EXISTS v_genome_variant_profile CASCADE;
DROP VIEW IF EXISTS v_mutation_frequency CASCADE;
DROP VIEW IF EXISTS v_mutation_lineage_breakdown CASCADE;
DROP VIEW IF EXISTS v_mutation_lineage_surveillance CASCADE;
DROP VIEW IF EXISTS v_mutation_lineage_country CASCADE;
DROP VIEW IF EXISTS v_mutation_geography CASCADE;
DROP VIEW IF EXISTS v_mutation_with_phenotype CASCADE;
DROP VIEW IF EXISTS v_mutation_surveillance CASCADE;
DROP VIEW IF EXISTS v_mutation_trends CASCADE;
DROP VIEW IF EXISTS v_mutation_co_occurrence CASCADE;
DROP VIEW IF EXISTS v_genomic_intelligence_snapshot CASCADE;
DROP VIEW IF EXISTS v_variant_recent_detection CASCADE;
DROP VIEW IF EXISTS v_phenotype_geo_temporal CASCADE;
DROP VIEW IF EXISTS v_phenotype_surveillance CASCADE;
DROP VIEW IF EXISTS v_species_summary CASCADE;

DROP MATERIALIZED VIEW IF EXISTS mv_mutation_frequency CASCADE;
DROP MATERIALIZED VIEW IF EXISTS mv_mutation_lineage_breakdown CASCADE;
DROP MATERIALIZED VIEW IF EXISTS mv_mutation_lineage_country CASCADE;
DROP MATERIALIZED VIEW IF EXISTS mv_mutation_summary CASCADE;
DROP MATERIALIZED VIEW IF EXISTS mv_mutation_with_phenotype CASCADE;
DROP MATERIALIZED VIEW IF EXISTS mv_phenotype_geo_temporal CASCADE;
DROP MATERIALIZED VIEW IF EXISTS mv_phenotype_surveillance CASCADE;
DROP MATERIALIZED VIEW IF EXISTS mv_position_variants CASCADE;

-- =============================================================================
-- 2. Drop the create_variant_view function (per-species/per-gene views)
-- =============================================================================

DROP FUNCTION IF EXISTS create_variant_view(VARCHAR, VARCHAR) CASCADE;
DROP FUNCTION IF EXISTS refresh_intelligence_matviews() CASCADE;

-- =============================================================================
-- 3. Drop indexes on protein_variants (they reference genome_accession)
-- =============================================================================

DROP INDEX IF EXISTS idx_proteinvars_gene;
DROP INDEX IF EXISTS idx_proteinvars_genome;
DROP INDEX IF EXISTS idx_proteinvars_hgvs;
DROP INDEX IF EXISTS idx_proteinvars_pathogen;
DROP INDEX IF EXISTS idx_proteinvars_position;
DROP INDEX IF EXISTS idx_proteinvars_reference;
DROP INDEX IF EXISTS idx_proteinvars_scope;
DROP INDEX IF EXISTS idx_proteinvars_scope_position;
DROP INDEX IF EXISTS idx_proteinvars_species;
DROP INDEX IF EXISTS idx_proteinvars_variant_lookup;
DROP INDEX IF EXISTS uq_proteinvars;

-- =============================================================================
-- 4. Drop FK constraints from protein_variants (genome_accession FK will go)
-- =============================================================================

ALTER TABLE protein_variants DROP CONSTRAINT IF EXISTS protein_variants_genome_accession_fkey;
ALTER TABLE protein_variants DROP CONSTRAINT IF EXISTS protein_variants_pathogen_id_fkey;
ALTER TABLE protein_variants DROP CONSTRAINT IF EXISTS protein_variants_reference_accession_fkey;
ALTER TABLE protein_variants DROP CONSTRAINT IF EXISTS protein_variants_species_id_fkey;

-- =============================================================================
-- 5. Truncate protein_variants (old per-genome data is being replaced)
-- =============================================================================

TRUNCATE TABLE protein_variants;

-- =============================================================================
-- 6. Alter protein_variants: remove genome_accession, add aggregated columns
-- =============================================================================

ALTER TABLE protein_variants DROP COLUMN IF EXISTS genome_accession;
ALTER TABLE protein_variants DROP COLUMN IF EXISTS nucleotide_ref;
ALTER TABLE protein_variants DROP COLUMN IF EXISTS nucleotide_alt;
ALTER TABLE protein_variants DROP COLUMN IF EXISTS is_synonymous;
ALTER TABLE protein_variants DROP COLUMN IF EXISTS is_frameshift;

-- Add aggregated columns
ALTER TABLE protein_variants ADD COLUMN genome_count INT NOT NULL DEFAULT 1;
ALTER TABLE protein_variants ADD COLUMN first_seen_date DATE;
ALTER TABLE protein_variants ADD COLUMN last_seen_date DATE;
ALTER TABLE protein_variants ADD COLUMN first_seen_year SMALLINT;
ALTER TABLE protein_variants ADD COLUMN last_seen_year SMALLINT;
ALTER TABLE protein_variants ADD COLUMN countries_seen TEXT[] DEFAULT '{}';
ALTER TABLE protein_variants ADD COLUMN country_codes TEXT[] DEFAULT '{}';
ALTER TABLE protein_variants ADD COLUMN lineage_ids TEXT[] DEFAULT '{}';
ALTER TABLE protein_variants ADD COLUMN species_total_genomes INT;

-- =============================================================================
-- 7. Add new unique constraint (one row per unique variant per species)
-- =============================================================================

CREATE UNIQUE INDEX uq_proteinvars_variant
    ON protein_variants (pathogen_id, species_id, gene, position, ref_aa, alt_aa, variant_type);

-- =============================================================================
-- 8. Add query-optimized indexes
-- =============================================================================

CREATE INDEX idx_proteinvars_scope ON protein_variants (pathogen_id, species_id, gene);
CREATE INDEX idx_proteinvars_variant_lookup ON protein_variants (species_id, gene, position, alt_aa);
CREATE INDEX idx_proteinvars_hgvs ON protein_variants (hgvs_p);
CREATE INDEX idx_proteinvars_frequency ON protein_variants (pathogen_id, species_id, genome_count DESC);

-- =============================================================================
-- 9. Create new lightweight views (no materialized views needed)
-- =============================================================================

-- v_variant_summary: One row per variant with frequency and geography
CREATE OR REPLACE VIEW v_variant_summary AS
SELECT
    pv.pathogen_id,
    pv.species_id,
    pv.gene,
    pv.position,
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    pv.is_stop,
    pv.genome_count,
    pv.species_total_genomes,
    round(pv.genome_count::numeric / NULLIF(pv.species_total_genomes, 0), 4) AS global_frequency,
    pv.first_seen_date,
    pv.last_seen_date,
    pv.first_seen_year,
    pv.last_seen_year,
    pv.countries_seen,
    pv.country_codes,
    pv.lineage_ids,
    gf.protein_name,
    gf.protein_function,
    gf.key_domains,
    gf.functional_sites,
    gf.known_hotspots
FROM protein_variants pv
LEFT JOIN gene_function gf
    ON pv.species_id = gf.species_id AND pv.gene = gf.gene;

-- v_variant_with_phenotype: Join variants to phenotype associations (at query time)
CREATE OR REPLACE VIEW v_variant_with_phenotype AS
SELECT
    pv.pathogen_id,
    pv.species_id,
    pv.gene,
    pv.position,
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    pv.genome_count,
    pv.species_total_genomes,
    round(pv.genome_count::numeric / NULLIF(pv.species_total_genomes, 0), 4) AS global_frequency,
    pv.first_seen_date,
    pv.last_seen_date,
    pv.first_seen_year,
    pv.last_seen_year,
    pv.countries_seen,
    pv.country_codes,
    pv.lineage_ids,
    gp.association_id,
    gp.phenotype_category,
    gp.phenotype_specific,
    gp.evidence_strength,
    gp.effect_size,
    gp.genotype_description,
    gp.literature_refs,
    gp.study_type,
    gp.verification_status
FROM protein_variants pv
LEFT JOIN genotype_phenotype gp ON
    pv.pathogen_id = gp.pathogen_id AND
    pv.species_id = gp.species_id AND
    (
        -- Exact variant match
        (pv.gene = gp.protein AND pv.position = gp.position AND pv.alt_aa = gp.alt_aa)
        -- Lineage-level association
        OR (gp.lineage_id IS NOT NULL AND gp.lineage_id = ANY(pv.lineage_ids))
        -- Gene-level association (position is NULL, protein is set)
        OR (gp.position IS NULL AND gp.protein IS NOT NULL AND gp.protein = pv.gene)
    );

-- v_species_summary: Genome counts per species/year/country (from genome_metadata)
CREATE OR REPLACE VIEW v_species_summary AS
SELECT
    pathogen_id,
    species_id,
    collection_year,
    collection_country,
    collection_country_code,
    count(*) AS total_genomes,
    min(collection_date) AS earliest_date,
    max(collection_date) AS latest_date
FROM genome_metadata
GROUP BY pathogen_id, species_id, collection_year, collection_country, collection_country_code;

-- v_phenotype_surveillance: Phenotype counts (directly from genotype_phenotype)
CREATE OR REPLACE VIEW v_phenotype_surveillance AS
SELECT
    pathogen_id,
    species_id,
    phenotype_category,
    evidence_strength,
    verification_status,
    CASE
        WHEN position IS NOT NULL THEN
            COALESCE((protein || ':' || ref_aa || position || alt_aa), genotype_description)
        ELSE genotype_description
    END AS genotype_label,
    protein,
    position,
    ref_aa,
    alt_aa,
    genotype_description,
    count(*) AS candidate_count
FROM genotype_phenotype
GROUP BY pathogen_id, species_id, phenotype_category, evidence_strength,
         verification_status, protein, position, ref_aa, alt_aa, genotype_description;

COMMIT;
