-- =============================================================================
-- PGIRL Database Cleanup Migration
-- =============================================================================
-- Drops tables not needed for the focused biological reference catalog:
--   mutations, outbreaks, country_context, disease_epidemiology, api_refresh_log
-- Drops all dependent views, materialized views, indexes, and FK constraints.
-- Trims operational/provenance columns from kept tables.
-- Removes the FK from reference_genomes.outbreak_id -> outbreaks.
-- =============================================================================

BEGIN;

-- =============================================================================
-- 1. Drop views that depend on tables being removed
-- =============================================================================

DROP VIEW IF EXISTS v_mutation_intelligence CASCADE;
DROP VIEW IF EXISTS v_disease_epidemiology CASCADE;
DROP VIEW IF EXISTS v_country_outbreak_burden CASCADE;
DROP VIEW IF EXISTS v_country_context_latest CASCADE;
DROP VIEW IF EXISTS v_curation_queue CASCADE;
DROP VIEW IF EXISTS v_geographic_lookup CASCADE;
DROP VIEW IF EXISTS v_lineage_outbreak_summary CASCADE;

-- =============================================================================
-- 2. Drop all 42 auto-generated protein_variants_{species}_{gene} views
-- =============================================================================

DROP VIEW IF EXISTS protein_variants_bdbv_gp CASCADE;
DROP VIEW IF EXISTS protein_variants_bdbv_l CASCADE;
DROP VIEW IF EXISTS protein_variants_bdbv_np CASCADE;
DROP VIEW IF EXISTS protein_variants_bdbv_vp24 CASCADE;
DROP VIEW IF EXISTS protein_variants_bdbv_vp30 CASCADE;
DROP VIEW IF EXISTS protein_variants_bdbv_vp35 CASCADE;
DROP VIEW IF EXISTS protein_variants_bdbv_vp40 CASCADE;
DROP VIEW IF EXISTS protein_variants_bomv_gp CASCADE;
DROP VIEW IF EXISTS protein_variants_bomv_l CASCADE;
DROP VIEW IF EXISTS protein_variants_bomv_np CASCADE;
DROP VIEW IF EXISTS protein_variants_bomv_vp24 CASCADE;
DROP VIEW IF EXISTS protein_variants_bomv_vp30 CASCADE;
DROP VIEW IF EXISTS protein_variants_bomv_vp35 CASCADE;
DROP VIEW IF EXISTS protein_variants_bomv_vp40 CASCADE;
DROP VIEW IF EXISTS protein_variants_ebov_gp CASCADE;
DROP VIEW IF EXISTS protein_variants_ebov_l CASCADE;
DROP VIEW IF EXISTS protein_variants_ebov_np CASCADE;
DROP VIEW IF EXISTS protein_variants_ebov_vp24 CASCADE;
DROP VIEW IF EXISTS protein_variants_ebov_vp30 CASCADE;
DROP VIEW IF EXISTS protein_variants_ebov_vp35 CASCADE;
DROP VIEW IF EXISTS protein_variants_ebov_vp40 CASCADE;
DROP VIEW IF EXISTS protein_variants_restv_gp CASCADE;
DROP VIEW IF EXISTS protein_variants_restv_l CASCADE;
DROP VIEW IF EXISTS protein_variants_restv_np CASCADE;
DROP VIEW IF EXISTS protein_variants_restv_vp24 CASCADE;
DROP VIEW IF EXISTS protein_variants_restv_vp30 CASCADE;
DROP VIEW IF EXISTS protein_variants_restv_vp35 CASCADE;
DROP VIEW IF EXISTS protein_variants_restv_vp40 CASCADE;
DROP VIEW IF EXISTS protein_variants_sudv_gp CASCADE;
DROP VIEW IF EXISTS protein_variants_sudv_l CASCADE;
DROP VIEW IF EXISTS protein_variants_sudv_np CASCADE;
DROP VIEW IF EXISTS protein_variants_sudv_vp24 CASCADE;
DROP VIEW IF EXISTS protein_variants_sudv_vp30 CASCADE;
DROP VIEW IF EXISTS protein_variants_sudv_vp35 CASCADE;
DROP VIEW IF EXISTS protein_variants_sudv_vp40 CASCADE;
DROP VIEW IF EXISTS protein_variants_tafv_gp CASCADE;
DROP VIEW IF EXISTS protein_variants_tafv_l CASCADE;
DROP VIEW IF EXISTS protein_variants_tafv_np CASCADE;
DROP VIEW IF EXISTS protein_variants_tafv_vp24 CASCADE;
DROP VIEW IF EXISTS protein_variants_tafv_vp30 CASCADE;
DROP VIEW IF EXISTS protein_variants_tafv_vp35 CASCADE;
DROP VIEW IF EXISTS protein_variants_tafv_vp40 CASCADE;

-- =============================================================================
-- 3. Drop tables no longer needed (CASCADE handles remaining dependents)
-- =============================================================================

DROP TABLE IF EXISTS api_refresh_log CASCADE;
DROP TABLE IF EXISTS country_context CASCADE;
DROP TABLE IF EXISTS disease_epidemiology CASCADE;
DROP TABLE IF EXISTS outbreaks CASCADE;
DROP TABLE IF EXISTS mutations CASCADE;

-- =============================================================================
-- 4. Remove FK constraint from reference_genomes -> outbreaks
--    (keep the outbreak_id column as plain text for reference)
-- =============================================================================

ALTER TABLE reference_genomes DROP CONSTRAINT IF EXISTS fk_refgenomes_outbreak;

-- =============================================================================
-- 5. Trim operational/provenance columns from kept tables
-- =============================================================================

-- genome_metadata: remove operational columns
ALTER TABLE genome_metadata DROP COLUMN IF EXISTS bioproject;
ALTER TABLE genome_metadata DROP COLUMN IF EXISTS biosample;
ALTER TABLE genome_metadata DROP COLUMN IF EXISTS submitter;
ALTER TABLE genome_metadata DROP COLUMN IF EXISTS source_url;
ALTER TABLE genome_metadata DROP COLUMN IF EXISTS data_source;
ALTER TABLE genome_metadata DROP COLUMN IF EXISTS notes;

-- genotype_phenotype: remove operational columns
ALTER TABLE genotype_phenotype DROP COLUMN IF EXISTS source_url;
ALTER TABLE genotype_phenotype DROP COLUMN IF EXISTS data_source;
ALTER TABLE genotype_phenotype DROP COLUMN IF EXISTS verified_by;
ALTER TABLE genotype_phenotype DROP COLUMN IF EXISTS verified_at;
ALTER TABLE genotype_phenotype DROP COLUMN IF EXISTS ingested_at;
ALTER TABLE genotype_phenotype DROP COLUMN IF EXISTS notes;

-- lineages: remove operational columns
ALTER TABLE lineages DROP COLUMN IF EXISTS source_url;
ALTER TABLE lineages DROP COLUMN IF EXISTS data_source;
ALTER TABLE lineages DROP COLUMN IF EXISTS verified_by;
ALTER TABLE lineages DROP COLUMN IF EXISTS verified_at;
ALTER TABLE lineages DROP COLUMN IF EXISTS ingested_at;
ALTER TABLE lineages DROP COLUMN IF EXISTS notes;
ALTER TABLE lineages DROP COLUMN IF EXISTS diagnostic_availability;
ALTER TABLE lineages DROP COLUMN IF EXISTS vaccine_availability;
ALTER TABLE lineages DROP COLUMN IF EXISTS therapeutics_available;

-- protein_variants: remove operational columns
ALTER TABLE protein_variants DROP COLUMN IF EXISTS data_source;
ALTER TABLE protein_variants DROP COLUMN IF EXISTS notes;

-- reference_genomes: remove operational columns
ALTER TABLE reference_genomes DROP COLUMN IF EXISTS source_url;
ALTER TABLE reference_genomes DROP COLUMN IF EXISTS data_source;
ALTER TABLE reference_genomes DROP COLUMN IF EXISTS verified_by;
ALTER TABLE reference_genomes DROP COLUMN IF EXISTS verified_at;
ALTER TABLE reference_genomes DROP COLUMN IF EXISTS ingested_at;
ALTER TABLE reference_genomes DROP COLUMN IF EXISTS notes;

-- =============================================================================
-- 6. Drop unused enum types
-- =============================================================================

DROP TYPE IF EXISTS circulation_status_enum CASCADE;
DROP TYPE IF EXISTS risk_tier_enum CASCADE;

-- verification_status_enum and evidence_strength_enum and frequency_enum
-- and genome_quality_enum are still used by kept tables — do NOT drop them.

COMMIT;
