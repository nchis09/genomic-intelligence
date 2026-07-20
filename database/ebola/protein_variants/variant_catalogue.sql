-- =============================================================================
-- VARIANT CATALOGUE — per-species reference-based protein variant calling.
--
-- ISOLATION MODEL:
--   Every table has pathogen_id as the FIRST column.
--   The intelligence engine ALWAYS queries with WHERE pathogen_id = 'X'.
--   This guarantees virus X never touches virus Y data.
--   Per-species/per-gene views add a second isolation layer.
--
--   When adding a new virus, duplicate this file under database/<virus>/mutation/
--   and change the pathogen_id values. The tables are shared but pathogen-scoped.
-- =============================================================================

-- ---------------------------------------------------------------------------
-- TABLE: reference_proteomes
-- One row per protein/gene product per reference genome.
-- Used as the mapping target for all variant calling.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS reference_proteomes (
    reference_accession     VARCHAR(30) NOT NULL REFERENCES reference_genomes(accession),
    species_id              VARCHAR(50) NOT NULL REFERENCES species(species_id),
    pathogen_id             VARCHAR(20) NOT NULL REFERENCES pathogens(pathogen_id),
    gene                    VARCHAR(30) NOT NULL,       -- e.g. "NP", "GP", "VP35"
    protein_name            VARCHAR(100),               -- full name
    protein_sequence        TEXT NOT NULL,              -- amino acid sequence of reference
    genome_start            INT NOT NULL,
    genome_end              INT NOT NULL,
    strand                  SMALLINT DEFAULT 1,
    protein_length          INT GENERATED ALWAYS AS (LENGTH(protein_sequence)) STORED,
    notes                   TEXT,
    last_curated            DATE,
    PRIMARY KEY (reference_accession, gene)
);

-- Composite index: pathogen first for scoping, then species, then gene
CREATE INDEX IF NOT EXISTS idx_refproteomes_scope
    ON reference_proteomes(pathogen_id, species_id, gene);

-- ---------------------------------------------------------------------------
-- TABLE: genome_metadata
-- One row per submitted/sequenced genome.
-- Links to the species reference it was aligned against.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS genome_metadata (
    genome_accession        VARCHAR(30) PRIMARY KEY,
    pathogen_id             VARCHAR(20) NOT NULL REFERENCES pathogens(pathogen_id),
    species_id              VARCHAR(50) NOT NULL REFERENCES species(species_id),
    reference_accession     VARCHAR(30) NOT NULL REFERENCES reference_genomes(accession),
    -- NCBI metadata
    ncbi_taxonomy_id        INT,
    bioproject              VARCHAR(30),
    biosample               VARCHAR(30),
    strain                  VARCHAR(255),
    isolate                 VARCHAR(255),
    -- Sample metadata
    collection_date         DATE,
    collection_year         SMALLINT,
    collection_country      VARCHAR(255),
    collection_country_code VARCHAR(10),
    collection_region       VARCHAR(255),
    host                    VARCHAR(255),
    isolation_source        VARCHAR(255),
    -- Sequence quality
    genome_length           INT,
    completeness            VARCHAR(30),    -- "complete" | "partial" | "nearly complete"
    genome_quality          genome_quality_enum DEFAULT 'MODERATE',
    -- Source/provenance
    source_db               VARCHAR(50) DEFAULT 'NCBI',
    source_url              TEXT,
    submitter               VARCHAR(255),
    release_date            DATE,
    data_source             VARCHAR(50) DEFAULT 'api_auto',
    last_updated            DATE,
    notes                   TEXT
);

-- Composite indexes: pathogen first for scoping
CREATE INDEX IF NOT EXISTS idx_genomemeta_scope
    ON genome_metadata(pathogen_id, species_id);
CREATE INDEX IF NOT EXISTS idx_genomemeta_pathogen_date
    ON genome_metadata(pathogen_id, collection_date);
CREATE INDEX IF NOT EXISTS idx_genomemeta_pathogen_country
    ON genome_metadata(pathogen_id, collection_country);

-- ---------------------------------------------------------------------------
-- TABLE: protein_variants
-- One row per amino acid variant per genome per gene.
-- This is the single source of truth for reference-based variant calling.
-- Per-species/per-gene views can be created on top.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS protein_variants (
    variant_id              BIGSERIAL PRIMARY KEY,
    pathogen_id             VARCHAR(20) NOT NULL REFERENCES pathogens(pathogen_id),
    species_id              VARCHAR(50) NOT NULL REFERENCES species(species_id),
    gene                    VARCHAR(30) NOT NULL,
    reference_accession     VARCHAR(30) NOT NULL REFERENCES reference_genomes(accession),
    genome_accession        VARCHAR(30) NOT NULL REFERENCES genome_metadata(genome_accession),
    position                INT NOT NULL,               -- 1-based position in reference protein
    ref_aa                  VARCHAR(5) NOT NULL,
    alt_aa                  VARCHAR(5) NOT NULL,
    variant_type            VARCHAR(30) NOT NULL,       -- "substitution" | "insertion" | "deletion" | "stop_gained" | "stop_lost" | "frameshift"
    hgvs_p                  VARCHAR(50),                -- e.g. "p.Ala82Val"
    nucleotide_ref          VARCHAR(100),               -- optional underlying nt context
    nucleotide_alt          VARCHAR(100),
    is_frameshift           BOOLEAN DEFAULT FALSE,
    is_stop                 BOOLEAN DEFAULT FALSE,
    is_synonymous           BOOLEAN DEFAULT FALSE,
    notes                   TEXT,
    data_source             VARCHAR(50) DEFAULT 'api_auto',
    last_updated            DATE
);

-- Composite indexes: pathogen first for scoping, then species + gene
CREATE INDEX IF NOT EXISTS idx_proteinvars_scope
    ON protein_variants(pathogen_id, species_id, gene);
CREATE INDEX IF NOT EXISTS idx_proteinvars_scope_position
    ON protein_variants(pathogen_id, species_id, gene, position);
CREATE INDEX IF NOT EXISTS idx_proteinvars_genome
    ON protein_variants(pathogen_id, genome_accession);
CREATE INDEX IF NOT EXISTS idx_proteinvars_hgvs
    ON protein_variants(pathogen_id, hgvs_p);

-- Unique constraint: pathogen-scoped to avoid duplicate variant calls per genome
CREATE UNIQUE INDEX IF NOT EXISTS uq_proteinvars
    ON protein_variants(pathogen_id, species_id, gene, genome_accession, position, variant_type, alt_aa);

-- ---------------------------------------------------------------------------
-- FUNCTION: create_variant_view(species_id, gene)
-- Creates a per-species/per-gene view that behaves like a separate table.
-- Example usage: SELECT create_variant_view('EBOV', 'GP');
-- ---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION create_variant_view(p_species_id VARCHAR(50), p_gene VARCHAR(30))
RETURNS TEXT AS $$
DECLARE
    view_name TEXT;
BEGIN
    view_name := 'protein_variants_' || LOWER(p_species_id) || '_' || LOWER(p_gene);
    EXECUTE format(
        'CREATE OR REPLACE VIEW %I AS
         SELECT *
         FROM protein_variants
         WHERE species_id = %L AND gene = %L',
        view_name, p_species_id, p_gene
    );
    RETURN view_name;
END;
$$ LANGUAGE plpgsql;
