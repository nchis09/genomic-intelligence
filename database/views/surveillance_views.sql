-- Surveillance views for the Genomic Epidemic Intelligence Engine.
-- All views are scoped by pathogen_id and species_id so they can be reused
-- for Ebola, Dengue, Influenza, etc. once those pipelines are built.
--
-- Load with:
--   psql pgirl -f database/views/surveillance_views.sql

-- =============================================================================
-- v_genome_variant_profile
-- Every variant call joined to its genome metadata. Use this as the base view
-- for matching a newly sequenced genome against the historical record.
-- =============================================================================
DROP VIEW IF EXISTS v_genome_variant_profile CASCADE;
CREATE VIEW v_genome_variant_profile AS
SELECT
    pv.pathogen_id,
    pv.species_id,
    pv.genome_accession,
    pv.gene,
    pv.position,
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    gm.collection_date,
    gm.collection_year,
    gm.collection_country,
    gm.collection_country_code,
    gm.collection_region,
    gm.lineage_id,
    gm.host,
    gm.completeness
FROM protein_variants pv
JOIN genome_metadata gm USING (pathogen_id, species_id, genome_accession);

-- =============================================================================
-- v_mutation_surveillance
-- Per-mutation, per-country, per-year counts. This is the core view for
-- frequency calculations and trend detection.
-- =============================================================================
CREATE OR REPLACE VIEW v_mutation_surveillance AS
SELECT
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    alt_aa,
    variant_type,
    hgvs_p,
    collection_year,
    collection_country,
    collection_country_code,
    COUNT(DISTINCT genome_accession) AS genome_count,
    COUNT(*) AS variant_observations
FROM v_genome_variant_profile
GROUP BY
    pathogen_id, species_id, gene, position, ref_aa, alt_aa, variant_type, hgvs_p,
    collection_year, collection_country, collection_country_code;

-- =============================================================================
-- v_mutation_geography
-- Where has each mutation been seen? First detection, last detection, and
-- list of countries. Also includes total genomes for that species for context.
-- =============================================================================
CREATE OR REPLACE VIEW v_mutation_geography AS
SELECT
    gvp.pathogen_id,
    gvp.species_id,
    gvp.gene,
    gvp.position,
    gvp.ref_aa,
    gvp.alt_aa,
    gvp.variant_type,
    gvp.hgvs_p,
    COUNT(DISTINCT gvp.genome_accession) AS genome_count,
    COUNT(*) AS variant_observations,
    MIN(gvp.collection_date) AS first_seen_date,
    MAX(gvp.collection_date) AS last_seen_date,
    MIN(gvp.collection_year) AS first_seen_year,
    MAX(gvp.collection_year) AS last_seen_year,
    ARRAY_AGG(DISTINCT gvp.collection_country ORDER BY gvp.collection_country) FILTER (WHERE gvp.collection_country IS NOT NULL) AS countries,
    ARRAY_AGG(DISTINCT gvp.collection_country_code ORDER BY gvp.collection_country_code) FILTER (WHERE gvp.collection_country_code IS NOT NULL) AS country_codes
FROM v_genome_variant_profile gvp
GROUP BY
    gvp.pathogen_id, gvp.species_id, gvp.gene, gvp.position, gvp.ref_aa, gvp.alt_aa, gvp.variant_type, gvp.hgvs_p;

-- =============================================================================
-- v_mutation_trends
-- Yearly counts per mutation. Useful for detecting emergence or expansion.
-- =============================================================================
CREATE OR REPLACE VIEW v_mutation_trends AS
SELECT
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    alt_aa,
    variant_type,
    hgvs_p,
    collection_year,
    COUNT(DISTINCT genome_accession) AS genome_count,
    COUNT(*) AS variant_observations
FROM v_genome_variant_profile
GROUP BY
    pathogen_id, species_id, gene, position, ref_aa, alt_aa, variant_type, hgvs_p, collection_year;

-- =============================================================================
-- v_mutation_frequency
-- Frequency of each mutation within its species and within each country.
-- country_frequency = genomes with mutation in country / total genomes in country.
-- global_frequency = genomes with mutation globally / total genomes for species.
-- =============================================================================
CREATE OR REPLACE VIEW v_mutation_frequency AS
WITH species_total AS (
    SELECT pathogen_id, species_id, COUNT(*) AS total_genomes
    FROM genome_metadata
    GROUP BY pathogen_id, species_id
),
country_total AS (
    SELECT pathogen_id, species_id, collection_country, collection_country_code, COUNT(*) AS total_genomes
    FROM genome_metadata
    GROUP BY pathogen_id, species_id, collection_country, collection_country_code
),
mutation_country AS (
    SELECT
        pathogen_id,
        species_id,
        gene,
        position,
        ref_aa,
        alt_aa,
        variant_type,
        hgvs_p,
        collection_country,
        collection_country_code,
        COUNT(DISTINCT genome_accession) AS genome_count
    FROM v_genome_variant_profile
    GROUP BY
        pathogen_id, species_id, gene, position, ref_aa, alt_aa, variant_type, hgvs_p,
        collection_country, collection_country_code
),
mutation_global AS (
    SELECT
        pathogen_id,
        species_id,
        gene,
        position,
        ref_aa,
        alt_aa,
        variant_type,
        hgvs_p,
        COUNT(DISTINCT genome_accession) AS genome_count
    FROM v_genome_variant_profile
    GROUP BY
        pathogen_id, species_id, gene, position, ref_aa, alt_aa, variant_type, hgvs_p
)
SELECT
    mc.pathogen_id,
    mc.species_id,
    mc.gene,
    mc.position,
    mc.ref_aa,
    mc.alt_aa,
    mc.variant_type,
    mc.hgvs_p,
    mc.collection_country,
    mc.collection_country_code,
    mc.genome_count AS country_genome_count,
    ct.total_genomes AS country_total_genomes,
    ROUND(mc.genome_count::numeric / NULLIF(ct.total_genomes, 0), 4) AS country_frequency,
    mg.genome_count AS global_genome_count,
    st.total_genomes AS species_total_genomes,
    ROUND(mg.genome_count::numeric / NULLIF(st.total_genomes, 0), 4) AS global_frequency
FROM mutation_country mc
JOIN country_total ct USING (pathogen_id, species_id, collection_country, collection_country_code)
JOIN mutation_global mg USING (pathogen_id, species_id, gene, position, ref_aa, alt_aa, variant_type, hgvs_p)
JOIN species_total st USING (pathogen_id, species_id);

-- =============================================================================
-- v_mutation_with_phenotype
-- Join every observed mutation to curated genotype-phenotype associations.
-- Matches on explicit mutation (protein/position/alt_aa), on lineage, or on
-- a protein-level/genotype_description pattern, so associations from literature
-- extraction that are not residue-specific are still surfaced.
-- =============================================================================
DROP VIEW IF EXISTS v_mutation_with_phenotype CASCADE;
CREATE VIEW v_mutation_with_phenotype AS
SELECT
    gvp.pathogen_id,
    gvp.species_id,
    gvp.genome_accession,
    gvp.gene,
    gvp.position,
    gvp.ref_aa,
    gvp.alt_aa,
    gvp.variant_type,
    gvp.hgvs_p,
    gvp.collection_date,
    gvp.collection_year,
    gvp.collection_country,
    gvp.collection_country_code,
    gvp.collection_region,
    gvp.lineage_id,
    gp.association_id,
    gp.phenotype_category,
    gp.phenotype_specific,
    gp.evidence_strength,
    gp.effect_size,
    gp.genotype_description,
    gp.literature_refs
FROM v_genome_variant_profile gvp
LEFT JOIN genotype_phenotype gp
    ON gvp.pathogen_id = gp.pathogen_id
    AND gvp.species_id = gp.species_id
    AND (
        -- Exact mutation match
        (gvp.gene = gp.protein AND gvp.position = gp.position AND gvp.alt_aa = gp.alt_aa)
        -- Lineage-level match (lineage_id is already resolved from strain/clade aliases
        -- and joined into v_genome_variant_profile)
        OR (gp.lineage_id IS NOT NULL AND gp.lineage_id = gvp.lineage_id)
        -- Protein/genotype description match (when position is null)
        OR (gp.position IS NULL AND gp.protein IS NOT NULL AND gp.protein = gvp.gene)
    );

-- Indexes that support the join strategies above
CREATE INDEX IF NOT EXISTS idx_gp_mutation_match
    ON genotype_phenotype (pathogen_id, species_id, protein, position, alt_aa);
CREATE INDEX IF NOT EXISTS idx_gp_lineage_match
    ON genotype_phenotype (pathogen_id, species_id, lineage_id)
    WHERE lineage_id IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_gp_data_source
    ON genotype_phenotype (data_source);

CREATE INDEX IF NOT EXISTS idx_genomemeta_lineage
    ON genome_metadata (pathogen_id, species_id, lineage_id)
    WHERE lineage_id IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_genomemeta_strain
    ON genome_metadata (strain text_pattern_ops);
CREATE INDEX IF NOT EXISTS idx_genomemeta_isolate
    ON genome_metadata (isolate text_pattern_ops);

CREATE INDEX IF NOT EXISTS idx_lineages_aliases
    ON lineages USING GIN (known_aliases);

-- =============================================================================
-- v_mutation_co_occurrence
-- Pairs of mutations observed in the same genome. Helps detect combinations
-- of potential public-health concern (e.g. GP + VP24 changes together).
-- =============================================================================
CREATE OR REPLACE VIEW v_mutation_co_occurrence AS
SELECT
    a.pathogen_id,
    a.species_id,
    a.genome_accession,
    a.gene AS gene_a,
    a.position AS position_a,
    a.alt_aa AS alt_aa_a,
    a.hgvs_p AS hgvs_p_a,
    b.gene AS gene_b,
    b.position AS position_b,
    b.alt_aa AS alt_aa_b,
    b.hgvs_p AS hgvs_p_b,
    gm.collection_year,
    gm.collection_country
FROM protein_variants a
JOIN protein_variants b
    ON a.pathogen_id = b.pathogen_id
    AND a.species_id = b.species_id
    AND a.genome_accession = b.genome_accession
    AND (a.gene, a.position, a.alt_aa) < (b.gene, b.position, b.alt_aa)
JOIN genome_metadata gm
    ON a.pathogen_id = gm.pathogen_id
    AND a.species_id = gm.species_id
    AND a.genome_accession = gm.genome_accession;

-- =============================================================================
-- v_species_summary
-- High-level counts per species/country/year for denominator/context queries.
-- =============================================================================
CREATE OR REPLACE VIEW v_species_summary AS
SELECT
    pathogen_id,
    species_id,
    collection_year,
    collection_country,
    collection_country_code,
    COUNT(*) AS total_genomes,
    MIN(collection_date) AS earliest_date,
    MAX(collection_date) AS latest_date
FROM genome_metadata
GROUP BY
    pathogen_id, species_id, collection_year, collection_country, collection_country_code;

-- =============================================================================
-- v_variant_recent_detection
-- Most recent observation of each mutation. Useful for flagging re-emergence.
-- =============================================================================
CREATE OR REPLACE VIEW v_variant_recent_detection AS
SELECT DISTINCT ON (pathogen_id, species_id, gene, position, ref_aa, alt_aa)
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    alt_aa,
    variant_type,
    hgvs_p,
    genome_accession,
    collection_date,
    collection_year,
    collection_country
FROM v_genome_variant_profile
ORDER BY pathogen_id, species_id, gene, position, ref_aa, alt_aa, collection_date DESC;

-- =============================================================================
-- MATERIALIZED VIEWS for fast intelligence-engine lookup
-- These pre-aggregate the surveillance data so the engine can query a mutation
-- in milliseconds instead of scanning the full variant table.
-- Refresh after every variant pipeline run:
--   REFRESH MATERIALIZED VIEW mv_mutation_summary;
--   REFRESH MATERIALIZED VIEW mv_mutation_frequency;
--   REFRESH MATERIALIZED VIEW mv_position_variants;
-- =============================================================================

-- =============================================================================
-- v_mutation_lineage_breakdown
-- Per-mutation counts split by lineage. This requires genome_metadata.lineage_id
-- to be populated (e.g. via Nextclade).
-- =============================================================================
CREATE OR REPLACE VIEW v_mutation_lineage_breakdown AS
SELECT
    pv.pathogen_id,
    pv.species_id,
    pv.gene,
    pv.position,
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    gm.lineage_id,
    COUNT(DISTINCT pv.genome_accession) AS genome_count,
    COUNT(*) AS variant_observations
FROM protein_variants pv
JOIN genome_metadata gm USING (pathogen_id, species_id, genome_accession)
WHERE gm.lineage_id IS NOT NULL
GROUP BY
    pv.pathogen_id, pv.species_id, pv.gene, pv.position, pv.ref_aa, pv.alt_aa,
    pv.variant_type, pv.hgvs_p, gm.lineage_id;


DROP MATERIALIZED VIEW IF EXISTS mv_mutation_lineage_breakdown CASCADE;
CREATE MATERIALIZED VIEW mv_mutation_lineage_breakdown AS
SELECT
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    alt_aa,
    variant_type,
    hgvs_p,
    lineage_id,
    genome_count,
    variant_observations
FROM v_mutation_lineage_breakdown;

CREATE UNIQUE INDEX idx_mv_mutation_lineage_breakdown_pk
    ON mv_mutation_lineage_breakdown (pathogen_id, species_id, gene, position, ref_aa, alt_aa, lineage_id);
CREATE INDEX idx_mv_mutation_lineage_breakdown_lookup
    ON mv_mutation_lineage_breakdown (species_id, gene, position, alt_aa, lineage_id);

-- =============================================================================
-- v_mutation_lineage_surveillance
-- Per-mutation counts by lineage, country and year. This is the drill-down view
-- for queries like "Which lineages/countries carry GP:A82V?".
-- =============================================================================
CREATE OR REPLACE VIEW v_mutation_lineage_surveillance AS
SELECT
    pv.pathogen_id,
    pv.species_id,
    pv.gene,
    pv.position,
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    gm.lineage_id,
    gm.collection_year,
    gm.collection_country,
    gm.collection_country_code,
    COUNT(DISTINCT pv.genome_accession) AS genome_count,
    COUNT(*) AS variant_observations
FROM protein_variants pv
JOIN genome_metadata gm USING (pathogen_id, species_id, genome_accession)
WHERE gm.lineage_id IS NOT NULL
GROUP BY
    pv.pathogen_id, pv.species_id, pv.gene, pv.position, pv.ref_aa, pv.alt_aa,
    pv.variant_type, pv.hgvs_p, gm.lineage_id,
    gm.collection_year, gm.collection_country, gm.collection_country_code;

DROP MATERIALIZED VIEW IF EXISTS mv_mutation_lineage_country CASCADE;
CREATE MATERIALIZED VIEW mv_mutation_lineage_country AS
SELECT
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    alt_aa,
    variant_type,
    hgvs_p,
    lineage_id,
    collection_year,
    collection_country,
    collection_country_code,
    genome_count,
    variant_observations
FROM v_mutation_lineage_surveillance;

CREATE UNIQUE INDEX idx_mv_mutation_lineage_country_pk
    ON mv_mutation_lineage_country (pathogen_id, species_id, gene, position, ref_aa, alt_aa, lineage_id, collection_country_code, collection_year);
CREATE INDEX idx_mv_mutation_lineage_country_lookup
    ON mv_mutation_lineage_country (species_id, gene, position, alt_aa, lineage_id, collection_country_code);

DROP MATERIALIZED VIEW IF EXISTS mv_mutation_summary CASCADE;
CREATE MATERIALIZED VIEW mv_mutation_summary AS
SELECT
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    alt_aa,
    variant_type,
    hgvs_p,
    genome_count,
    variant_observations,
    first_seen_date,
    last_seen_date,
    first_seen_year,
    last_seen_year,
    countries,
    country_codes
FROM v_mutation_geography;

CREATE UNIQUE INDEX idx_mv_mutation_summary_pk
    ON mv_mutation_summary (pathogen_id, species_id, gene, position, ref_aa, alt_aa);
CREATE INDEX idx_mv_mutation_summary_hgvs
    ON mv_mutation_summary (hgvs_p);
CREATE INDEX idx_mv_mutation_summary_species
    ON mv_mutation_summary (species_id, gene, position, alt_aa);

DROP MATERIALIZED VIEW IF EXISTS mv_mutation_frequency CASCADE;
CREATE MATERIALIZED VIEW mv_mutation_frequency AS
SELECT
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    alt_aa,
    variant_type,
    hgvs_p,
    collection_country,
    collection_country_code,
    country_genome_count,
    country_total_genomes,
    country_frequency,
    global_genome_count,
    species_total_genomes,
    global_frequency
FROM v_mutation_frequency;

CREATE UNIQUE INDEX idx_mv_mutation_frequency_pk
    ON mv_mutation_frequency (pathogen_id, species_id, gene, position, ref_aa, alt_aa, collection_country_code);
CREATE INDEX idx_mv_mutation_frequency_country
    ON mv_mutation_frequency (collection_country_code, species_id, gene, position, alt_aa);
CREATE INDEX idx_mv_mutation_frequency_global
    ON mv_mutation_frequency (species_id, gene, position, alt_aa, global_frequency);

DROP MATERIALIZED VIEW IF EXISTS mv_position_variants CASCADE;
CREATE MATERIALIZED VIEW mv_position_variants AS
SELECT
    pathogen_id,
    species_id,
    gene,
    position,
    ref_aa,
    COUNT(DISTINCT alt_aa) AS distinct_alt_count,
    ARRAY_AGG(DISTINCT alt_aa ORDER BY alt_aa) AS observed_alt_aas,
    SUM(CASE WHEN alt_aa = ref_aa THEN 0 ELSE 1 END) AS non_ref_variant_observations
FROM v_genome_variant_profile
GROUP BY pathogen_id, species_id, gene, position, ref_aa;

CREATE UNIQUE INDEX idx_mv_position_variants_pk
    ON mv_position_variants (pathogen_id, species_id, gene, position, ref_aa);
CREATE INDEX idx_mv_position_variants_species_pos
    ON mv_position_variants (species_id, gene, position);

-- =============================================================================
-- v_phenotype_surveillance
-- Aggregate of candidate genotype-phenotype associations from literature.
-- Supports high-level queries like "which phenotypes are linked to GP mutations"
-- or "how many unverified candidates mention immune escape".
-- =============================================================================
DROP VIEW IF EXISTS v_phenotype_surveillance CASCADE;
CREATE VIEW v_phenotype_surveillance AS
SELECT
    pathogen_id,
    species_id,
    phenotype_category,
    evidence_strength,
    verification_status,
    CASE
        WHEN position IS NOT NULL THEN COALESCE(protein || ':' || ref_aa || position || alt_aa, genotype_description)
        ELSE genotype_description
    END AS genotype_label,
    protein,
    position,
    ref_aa,
    alt_aa,
    genotype_description,
    COUNT(*) AS candidate_count
FROM genotype_phenotype
GROUP BY
    pathogen_id, species_id, phenotype_category,
    evidence_strength, verification_status, protein, position, ref_aa, alt_aa,
    genotype_description;

DROP MATERIALIZED VIEW IF EXISTS mv_phenotype_surveillance CASCADE;
CREATE MATERIALIZED VIEW mv_phenotype_surveillance AS
SELECT * FROM v_phenotype_surveillance;

CREATE UNIQUE INDEX idx_mv_phenotype_surveillance_pk
    ON mv_phenotype_surveillance (pathogen_id, species_id, phenotype_category, evidence_strength, verification_status, genotype_label);
CREATE INDEX idx_mv_phenotype_surveillance_category
    ON mv_phenotype_surveillance (species_id, phenotype_category, verification_status);
CREATE INDEX idx_mv_phenotype_surveillance_genotype
    ON mv_phenotype_surveillance (species_id, protein, position, alt_aa, genotype_description);

-- =============================================================================
-- v_mutation_with_phenotype already exists above; this is its fast lookup
-- materialized version for the intelligence engine.
-- =============================================================================
DROP MATERIALIZED VIEW IF EXISTS mv_mutation_with_phenotype CASCADE;
CREATE MATERIALIZED VIEW mv_mutation_with_phenotype AS
SELECT * FROM v_mutation_with_phenotype;

CREATE UNIQUE INDEX idx_mv_mutation_with_phenotype_pk
    ON mv_mutation_with_phenotype (pathogen_id, species_id, genome_accession, gene, position, alt_aa, association_id);
CREATE INDEX idx_mv_mutation_with_phenotype_lookup
    ON mv_mutation_with_phenotype (species_id, gene, position, alt_aa, phenotype_category);

-- =============================================================================
-- v_phenotype_geo_temporal
-- Where and when has each phenotype-genotype association been observed?
-- Joins observed variants (and lineage/motif matches) to genome metadata so
-- you can ask "which phenotypes are common in this region/clade and how long
-- have they been present?".
-- =============================================================================
DROP VIEW IF EXISTS v_phenotype_geo_temporal CASCADE;
CREATE VIEW v_phenotype_geo_temporal AS
SELECT
    mwp.pathogen_id,
    mwp.species_id,
    mwp.association_id,
    mwp.phenotype_category,
    mwp.phenotype_specific,
    mwp.evidence_strength,
    mwp.gene AS protein,
    mwp.genotype_description,
    mwp.lineage_id,
    mwp.collection_year,
    mwp.collection_country,
    mwp.collection_country_code,
    MAX(mwp.collection_region) AS collection_region,
    COUNT(DISTINCT mwp.genome_accession) AS genome_count,
    COUNT(*) AS variant_observations,
    MIN(mwp.collection_date) AS first_seen_date,
    MAX(mwp.collection_date) AS last_seen_date,
    MIN(mwp.collection_year) AS first_seen_year,
    MAX(mwp.collection_year) AS last_seen_year
FROM v_mutation_with_phenotype mwp
WHERE mwp.association_id IS NOT NULL
GROUP BY
    mwp.pathogen_id, mwp.species_id, mwp.association_id, mwp.phenotype_category,
    mwp.phenotype_specific, mwp.evidence_strength,
    mwp.gene, mwp.genotype_description, mwp.lineage_id, mwp.collection_year,
    mwp.collection_country, mwp.collection_country_code;

CREATE MATERIALIZED VIEW mv_phenotype_geo_temporal AS
SELECT * FROM v_phenotype_geo_temporal;

CREATE UNIQUE INDEX idx_mv_phenotype_geo_temporal_pk
    ON mv_phenotype_geo_temporal (pathogen_id, species_id, association_id, lineage_id, collection_country, collection_country_code, collection_year, protein);
CREATE INDEX idx_mv_phenotype_geo_temporal_phenotype
    ON mv_phenotype_geo_temporal (species_id, phenotype_category, collection_country_code);
CREATE INDEX idx_mv_phenotype_geo_temporal_lineage
    ON mv_phenotype_geo_temporal (species_id, lineage_id, phenotype_category);
CREATE INDEX idx_mv_phenotype_geo_temporal_genotype
    ON mv_phenotype_geo_temporal (species_id, protein, genotype_description, phenotype_category);
