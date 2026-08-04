-- Genomic Intelligence Knowledge Warehouse schema
--
-- Designed to be pathogen-agnostic. Tables model biological and
-- epidemiological entities, their relationships, annotations, and provenance.
--
-- Layer 1: core biological entities
-- Layer 2: knowledge annotations
-- Layer 3: epidemiological context
-- Layer 4: evidence and provenance

-- Layer 1: core biological entities

CREATE TABLE IF NOT EXISTS analysis_runs (
    run_id TEXT PRIMARY KEY,
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id SERIAL PRIMARY KEY,
    run_id TEXT REFERENCES analysis_runs(run_id),
    sample_name TEXT NOT NULL,
    pathogen TEXT,
    species TEXT,
    qc_score REAL,
    best_dataset_file TEXT,
    collection_date DATE,
    country TEXT,
    admin1 TEXT,
    admin2 TEXT,
    locality TEXT,
    host TEXT,
    outbreak TEXT,
    strain TEXT,
    lineage TEXT,
    clade TEXT,
    ppx_accession TEXT,
    insdc_accession TEXT,
    genome_coverage REAL,
    nuc_substitution_count INTEGER,
    aa_mutation_count INTEGER,
    nextclade_qc TEXT,
    is_query BOOLEAN DEFAULT FALSE,
    UNIQUE(run_id, sample_name)
);

CREATE TABLE IF NOT EXISTS genomes (
    genome_id SERIAL PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(sample_id),
    accession TEXT,
    length INTEGER,
    sequence TEXT,
    UNIQUE(sample_id, accession)
);

CREATE TABLE IF NOT EXISTS reference_genomes (
    ref_id SERIAL PRIMARY KEY,
    pathogen TEXT,
    species TEXT,
    accession TEXT,
    source TEXT,
    length INTEGER,
    sequence TEXT,
    UNIQUE(pathogen, species, accession)
);

CREATE TABLE IF NOT EXISTS genes (
    gene_id SERIAL PRIMARY KEY,
    ref_id INTEGER REFERENCES reference_genomes(ref_id),
    gene_name TEXT,
    product TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    strand TEXT
);

CREATE TABLE IF NOT EXISTS proteins (
    protein_id SERIAL PRIMARY KEY,
    gene_id INTEGER REFERENCES genes(gene_id),
    uniprot_accession TEXT,
    protein_name TEXT,
    organism TEXT,
    length INTEGER,
    UNIQUE(gene_id, uniprot_accession)
);

CREATE TABLE IF NOT EXISTS mutations (
    mutation_id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES proteins(protein_id),
    gene_id INTEGER REFERENCES genes(gene_id),
    mutation_label TEXT NOT NULL,
    ref_aa TEXT,
    position INTEGER,
    alt_aa TEXT,
    mutation_type TEXT,
    UNIQUE(protein_id, mutation_label)
);

CREATE TABLE IF NOT EXISTS phylogenetic_trees (
    tree_id SERIAL PRIMARY KEY,
    run_id TEXT REFERENCES analysis_runs(run_id),
    species TEXT,
    pathogen TEXT,
    tree_source TEXT,
    newick TEXT
);

CREATE TABLE IF NOT EXISTS tree_tips (
    tip_id SERIAL PRIMARY KEY,
    tree_id INTEGER REFERENCES phylogenetic_trees(tree_id),
    sample_id INTEGER REFERENCES samples(sample_id),
    label TEXT NOT NULL,
    is_query BOOLEAN DEFAULT FALSE,
    ppx_accession TEXT,
    insdc_accession TEXT,
    country TEXT,
    admin1 TEXT,
    admin2 TEXT,
    locality TEXT,
    host TEXT,
    outbreak TEXT,
    tip_date DATE,
    div REAL,
    genome_coverage REAL,
    nextclade_qc TEXT,
    aa_mutation_count INTEGER,
    nuc_mutation_count INTEGER
);

CREATE TABLE IF NOT EXISTS outbreaks (
    outbreak_id SERIAL PRIMARY KEY,
    run_id TEXT REFERENCES analysis_runs(run_id),
    name TEXT,
    pathogen TEXT,
    species TEXT,
    location TEXT,
    date_start DATE,
    date_end DATE
);

-- Layer 2: knowledge annotations

CREATE TABLE IF NOT EXISTS protein_functions (
    function_id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES proteins(protein_id),
    source TEXT,
    function_text TEXT
);

CREATE TABLE IF NOT EXISTS protein_domains (
    domain_id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES proteins(protein_id),
    domain_name TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    source TEXT,
    hmm_accession TEXT,
    evalue DOUBLE PRECISION,
    bit_score DOUBLE PRECISION,
    hmm_start INTEGER,
    hmm_end INTEGER,
    ali_start INTEGER,
    ali_end INTEGER,
    env_start INTEGER,
    env_end INTEGER
);

CREATE TABLE IF NOT EXISTS go_terms (
    go_term_id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES proteins(protein_id),
    go_id TEXT,
    go_term TEXT,
    go_category TEXT,
    source TEXT
);

CREATE TABLE IF NOT EXISTS biological_pathways (
    pathway_id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES proteins(protein_id),
    pathway_name TEXT,
    source TEXT
);

CREATE TABLE IF NOT EXISTS protein_interactions (
    interaction_id SERIAL PRIMARY KEY,
    protein_a_id INTEGER REFERENCES proteins(protein_id),
    protein_b_id INTEGER REFERENCES proteins(protein_id),
    score REAL,
    source TEXT
);

CREATE TABLE IF NOT EXISTS mutation_phenotypes (
    phenotype_id SERIAL PRIMARY KEY,
    mutation_id INTEGER REFERENCES mutations(mutation_id),
    phenotype TEXT,
    effect TEXT,
    evidence TEXT,
    source TEXT
);

CREATE TABLE IF NOT EXISTS structural_features (
    feature_id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES proteins(protein_id),
    feature_type TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    description TEXT,
    source TEXT
);

-- Layer 3: epidemiological context

CREATE TABLE IF NOT EXISTS geographic_locations (
    location_id SERIAL PRIMARY KEY,
    country TEXT,
    admin1 TEXT,
    admin2 TEXT,
    locality TEXT,
    location_code TEXT,
    location_code_type TEXT,
    latitude REAL,
    longitude REAL,
    UNIQUE(country, admin1, admin2, locality, location_code)
);

CREATE TABLE IF NOT EXISTS epidemiological_datasets (
    dataset_id SERIAL PRIMARY KEY,
    run_id TEXT REFERENCES analysis_runs(run_id),
    dataset_name TEXT,
    source TEXT,
    retrieval_date DATE,
    row_count INTEGER,
    search_rank INTEGER,
    pathogen TEXT,
    species TEXT,
    dataset_type TEXT
);

CREATE TABLE IF NOT EXISTS epidemiological_records (
    record_id SERIAL PRIMARY KEY,
    dataset_id INTEGER REFERENCES epidemiological_datasets(dataset_id),
    record_date DATE,
    report_date DATE,
    reference_date DATE,
    location_id INTEGER REFERENCES geographic_locations(location_id),
    country TEXT,
    admin1 TEXT,
    admin2 TEXT,
    locality TEXT,
    location_code TEXT,
    location_code_type TEXT,
    measure TEXT,
    case_classification TEXT,
    time_period TEXT,
    value REAL,
    unit TEXT,
    cases INTEGER,
    deaths INTEGER,
    suspected INTEGER,
    recovered INTEGER,
    gender_breakdown JSONB,
    age_breakdown JSONB,
    source_url TEXT,
    source_indicator_name TEXT,
    indicator_label TEXT,
    raw_data JSONB
);

CREATE TABLE IF NOT EXISTS surveillance_records (
    surveillance_id SERIAL PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(sample_id),
    location_id INTEGER REFERENCES geographic_locations(location_id),
    record_date DATE
);

CREATE TABLE IF NOT EXISTS transmission_events (
    event_id SERIAL PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(sample_id),
    related_sample_id INTEGER REFERENCES samples(sample_id),
    event_date DATE,
    location_id INTEGER REFERENCES geographic_locations(location_id),
    confidence TEXT
);

-- Layer 4: evidence and provenance

CREATE TABLE IF NOT EXISTS evidence_sources (
    source_id SERIAL PRIMARY KEY,
    source_name TEXT NOT NULL,
    version TEXT,
    url TEXT,
    retrieval_date DATE
);

CREATE TABLE IF NOT EXISTS entity_evidence (
    entity_evidence_id SERIAL PRIMARY KEY,
    table_name TEXT NOT NULL,
    entity_id INTEGER NOT NULL,
    source_id INTEGER REFERENCES evidence_sources(source_id),
    confidence TEXT,
    reference TEXT,
    retrieved_date DATE
);

CREATE TABLE IF NOT EXISTS database_versions (
    db_id SERIAL PRIMARY KEY,
    database_name TEXT,
    version TEXT,
    retrieval_date DATE
);

-- Layer 4b: pipeline output provenance

CREATE TABLE IF NOT EXISTS pipeline_outputs (
    output_id SERIAL PRIMARY KEY,
    run_id TEXT REFERENCES analysis_runs(run_id),
    stage TEXT,
    process_name TEXT,
    file_type TEXT,
    file_name TEXT,
    file_path TEXT,
    file_hash TEXT,
    row_count INTEGER,
    loaded_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    notes TEXT
);

-- Relationship tables linking entities explicitly

CREATE TABLE IF NOT EXISTS sample_mutation (
    sample_id INTEGER REFERENCES samples(sample_id),
    mutation_id INTEGER REFERENCES mutations(mutation_id),
    PRIMARY KEY (sample_id, mutation_id)
);

CREATE TABLE IF NOT EXISTS protein_mutation (
    protein_id INTEGER REFERENCES proteins(protein_id),
    mutation_id INTEGER REFERENCES mutations(mutation_id),
    PRIMARY KEY (protein_id, mutation_id)
);

CREATE TABLE IF NOT EXISTS protein_annotation (
    protein_id INTEGER REFERENCES proteins(protein_id),
    annotation_id INTEGER,
    annotation_table TEXT,
    PRIMARY KEY (protein_id, annotation_id, annotation_table)
);

CREATE TABLE IF NOT EXISTS mutation_phenotype_evidence (
    phenotype_id INTEGER REFERENCES mutation_phenotypes(phenotype_id),
    source_id INTEGER REFERENCES evidence_sources(source_id),
    PRIMARY KEY (phenotype_id, source_id)
);

CREATE TABLE IF NOT EXISTS sample_outbreak (
    sample_id INTEGER REFERENCES samples(sample_id),
    outbreak_id INTEGER REFERENCES outbreaks(outbreak_id),
    PRIMARY KEY (sample_id, outbreak_id)
);

CREATE TABLE IF NOT EXISTS sample_geo_location (
    sample_id INTEGER REFERENCES samples(sample_id),
    location_id INTEGER REFERENCES geographic_locations(location_id),
    PRIMARY KEY (sample_id, location_id)
);

-- Indexes for query-driven public health searches

CREATE INDEX IF NOT EXISTS idx_samples_run_id ON samples(run_id);
CREATE INDEX IF NOT EXISTS idx_samples_sample_name ON samples(sample_name);
CREATE INDEX IF NOT EXISTS idx_samples_pathogen ON samples(pathogen);
CREATE INDEX IF NOT EXISTS idx_samples_species ON samples(species);
CREATE INDEX IF NOT EXISTS idx_samples_country ON samples(country);
CREATE INDEX IF NOT EXISTS idx_samples_collection_date ON samples(collection_date);
CREATE INDEX IF NOT EXISTS idx_samples_outbreak ON samples(outbreak);
CREATE INDEX IF NOT EXISTS idx_samples_lineage ON samples(lineage);

CREATE INDEX IF NOT EXISTS idx_genomes_sample_id ON genomes(sample_id);
CREATE INDEX IF NOT EXISTS idx_genomes_accession ON genomes(accession);

CREATE INDEX IF NOT EXISTS idx_proteins_uniprot ON proteins(uniprot_accession);
CREATE INDEX IF NOT EXISTS idx_proteins_gene_id ON proteins(gene_id);

CREATE INDEX IF NOT EXISTS idx_mutations_label ON mutations(mutation_label);
CREATE INDEX IF NOT EXISTS idx_mutations_protein_id ON mutations(protein_id);
CREATE INDEX IF NOT EXISTS idx_mutations_position ON mutations(position);

CREATE INDEX IF NOT EXISTS idx_epi_records_dataset_id ON epidemiological_records(dataset_id);
CREATE INDEX IF NOT EXISTS idx_epi_records_date ON epidemiological_records(record_date);
CREATE INDEX IF NOT EXISTS idx_epi_records_country ON epidemiological_records(country);
CREATE INDEX IF NOT EXISTS idx_epi_records_location_id ON epidemiological_records(location_id);

CREATE INDEX IF NOT EXISTS idx_geo_country ON geographic_locations(country);
CREATE INDEX IF NOT EXISTS idx_geo_admin1 ON geographic_locations(admin1);

CREATE INDEX IF NOT EXISTS idx_outbreaks_name ON outbreaks(name);
CREATE INDEX IF NOT EXISTS idx_outbreaks_pathogen ON outbreaks(pathogen);

CREATE INDEX IF NOT EXISTS idx_entity_evidence_table_entity ON entity_evidence(table_name, entity_id);
CREATE INDEX IF NOT EXISTS idx_evidence_sources_name ON evidence_sources(source_name);

CREATE INDEX IF NOT EXISTS idx_samples_admin1 ON samples(admin1);
CREATE INDEX IF NOT EXISTS idx_samples_locality ON samples(locality);
CREATE INDEX IF NOT EXISTS idx_samples_ppx ON samples(ppx_accession);
CREATE INDEX IF NOT EXISTS idx_samples_is_query ON samples(is_query);
CREATE INDEX IF NOT EXISTS idx_tree_tips_sample_id ON tree_tips(sample_id);
CREATE INDEX IF NOT EXISTS idx_tree_tips_label ON tree_tips(label);
CREATE INDEX IF NOT EXISTS idx_tree_tips_country_date ON tree_tips(country, tip_date);
CREATE INDEX IF NOT EXISTS idx_tree_tips_outbreak ON tree_tips(outbreak);
CREATE INDEX IF NOT EXISTS idx_geo_location_code ON geographic_locations(location_code, location_code_type);
CREATE INDEX IF NOT EXISTS idx_epi_records_dataset_id ON epidemiological_records(dataset_id);
CREATE INDEX IF NOT EXISTS idx_epi_records_country_date_measure ON epidemiological_records(country, record_date, measure);
CREATE INDEX IF NOT EXISTS idx_epi_records_measure_classification ON epidemiological_records(measure, case_classification);
CREATE INDEX IF NOT EXISTS idx_epi_records_location_id ON epidemiological_records(location_id);
CREATE INDEX IF NOT EXISTS idx_pipeline_outputs_run_id ON pipeline_outputs(run_id);
CREATE INDEX IF NOT EXISTS idx_pipeline_outputs_stage ON pipeline_outputs(stage);
CREATE INDEX IF NOT EXISTS idx_pipeline_outputs_file_name ON pipeline_outputs(file_name);
