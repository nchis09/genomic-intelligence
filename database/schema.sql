--
-- PostgreSQL database dump
--

\restrict pzG6puWa0aiNZRlVp6flk0RhHcWllMbYxAEhdVTRGB4mwZ8eRah6ZG8DfHcljhR

-- Dumped from database version 16.14 (Homebrew)
-- Dumped by pg_dump version 16.14 (Homebrew)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: public; Type: SCHEMA; Schema: -; Owner: -
--

-- *not* creating schema, since initdb creates it


--
-- Name: SCHEMA public; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON SCHEMA public IS '';


--
-- Name: pg_stat_statements; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pg_stat_statements WITH SCHEMA public;


--
-- Name: EXTENSION pg_stat_statements; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON EXTENSION pg_stat_statements IS 'track planning and execution statistics of all SQL statements executed';


--
-- Name: evidence_strength_enum; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.evidence_strength_enum AS ENUM (
    'strong',
    'moderate',
    'weak',
    'preliminary',
    'insufficient'
);


--
-- Name: frequency_enum; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.frequency_enum AS ENUM (
    'Common',
    'Uncommon',
    'Rare',
    'Novel'
);


--
-- Name: genome_quality_enum; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.genome_quality_enum AS ENUM (
    'HIGH',
    'MODERATE',
    'LOW'
);


--
-- Name: verification_status_enum; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.verification_status_enum AS ENUM (
    'unverified',
    'verified',
    'rejected',
    'needs_review'
);


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: gene_function; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene_function (
    gene_function_id integer NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    species_id character varying(50) NOT NULL,
    gene character varying(20) NOT NULL,
    protein_name character varying(100),
    protein_function text,
    genome_start integer,
    genome_end integer,
    protein_length_aa integer,
    key_domains jsonb,
    functional_sites jsonb,
    known_hotspots jsonb,
    conserved_regions jsonb,
    pdb_ids text[],
    literature_refs text[],
    last_curated date,
    curator character varying(100),
    protein_function_tsv tsvector GENERATED ALWAYS AS (to_tsvector('english'::regconfig, COALESCE(protein_function, ''::text))) STORED
);


--
-- Name: gene_function_gene_function_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.gene_function_gene_function_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: gene_function_gene_function_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.gene_function_gene_function_id_seq OWNED BY public.gene_function.gene_function_id;


--
-- Name: genome_metadata; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genome_metadata (
    genome_accession character varying(30) NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    species_id character varying(50) NOT NULL,
    reference_accession character varying(30) NOT NULL,
    ncbi_taxonomy_id integer,
    strain character varying(255),
    isolate character varying(255),
    collection_date date,
    collection_year smallint,
    collection_country character varying(255),
    collection_country_code character varying(10),
    collection_region character varying(255),
    host character varying(255),
    isolation_source character varying(255),
    genome_length integer,
    completeness character varying(30),
    genome_quality public.genome_quality_enum DEFAULT 'MODERATE'::public.genome_quality_enum,
    source_db character varying(50) DEFAULT 'NCBI'::character varying,
    release_date date,
    last_updated date,
    lineage_id character varying(50)
);


--
-- Name: genotype_phenotype; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genotype_phenotype (
    association_id character varying(20) NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    species_id character varying(50) NOT NULL,
    lineage_id character varying(50),
    mutation_id character varying(20),
    protein character varying(20),
    "position" integer,
    ref_aa character varying(30),
    alt_aa character varying(30),
    genotype_description text,
    phenotype_category text,
    phenotype_specific text,
    effect_size text,
    study_type text,
    model_system text,
    evidence_strength public.evidence_strength_enum NOT NULL,
    first_reported_year smallint,
    confirmed_by text,
    conflicted_by text,
    literature_refs text[],
    record_flagged boolean DEFAULT false,
    flag_reason text,
    last_updated date,
    verification_status public.verification_status_enum DEFAULT 'unverified'::public.verification_status_enum
);


--
-- Name: lineages; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.lineages (
    lineage_id character varying(50) NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    species_id character varying(50) NOT NULL,
    lineage_name character varying(100) NOT NULL,
    clade character varying(100),
    event_classification character varying(100),
    first_country_detected character varying(200),
    first_region_detected character varying(200),
    countries_reported text[],
    regions_reported text[],
    current_distribution text,
    endemic_regions text[],
    first_detected smallint,
    last_detected smallint,
    primary_host character varying(150),
    reservoir character varying(150),
    human_to_human boolean,
    animal_to_human boolean,
    vector_borne boolean DEFAULT false,
    nosocomial_transmission boolean,
    number_genomes_available integer,
    known_recombination boolean DEFAULT false,
    known_reassortment boolean DEFAULT false,
    evolutionary_rate character varying(50),
    who_priority_pathogen boolean DEFAULT false,
    last_updated date,
    verification_status public.verification_status_enum DEFAULT 'unverified'::public.verification_status_enum,
    known_aliases text[]
);


--
-- Name: pathogens; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.pathogens (
    pathogen_id character varying(20) NOT NULL,
    family character varying(100) NOT NULL,
    genus character varying(100) NOT NULL,
    ncbi_taxonomy_id integer,
    notes text,
    last_curated date,
    curator character varying(100)
);


--
-- Name: protein_variants; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.protein_variants (
    variant_id bigint NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    species_id character varying(50) NOT NULL,
    gene character varying(30) NOT NULL,
    reference_accession character varying(30) NOT NULL,
    "position" integer NOT NULL,
    ref_aa character varying(5) NOT NULL,
    alt_aa character varying(5) NOT NULL,
    variant_type character varying(30) NOT NULL,
    hgvs_p character varying(50),
    is_stop boolean DEFAULT false,
    last_updated date,
    genome_count integer DEFAULT 1 NOT NULL,
    first_seen_date date,
    last_seen_date date,
    first_seen_year smallint,
    last_seen_year smallint,
    countries_seen text[] DEFAULT '{}'::text[],
    country_codes text[] DEFAULT '{}'::text[],
    lineage_ids text[] DEFAULT '{}'::text[],
    species_total_genomes integer
);


--
-- Name: protein_variants_variant_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.protein_variants_variant_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: protein_variants_variant_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.protein_variants_variant_id_seq OWNED BY public.protein_variants.variant_id;


--
-- Name: reference_genomes; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.reference_genomes (
    accession character varying(30) NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    species_id character varying(50) NOT NULL,
    genome_role character varying(50),
    genome_length integer,
    segmented boolean DEFAULT false,
    collection_year smallint,
    collection_country character varying(100),
    outbreak_id character varying(50),
    source_database character varying(50) DEFAULT 'NCBI'::character varying,
    gene_coordinates jsonb,
    notes text,
    last_curated date,
    curator character varying(100),
    verification_status public.verification_status_enum DEFAULT 'unverified'::public.verification_status_enum
);


--
-- Name: reference_proteomes; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.reference_proteomes (
    reference_accession character varying(30) NOT NULL,
    species_id character varying(50) NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    gene character varying(30) NOT NULL,
    protein_name character varying(100),
    protein_sequence text NOT NULL,
    genome_start integer NOT NULL,
    genome_end integer NOT NULL,
    strand smallint DEFAULT 1,
    protein_length integer GENERATED ALWAYS AS (length(protein_sequence)) STORED,
    notes text,
    last_curated date
);


--
-- Name: species; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.species (
    species_id character varying(50) NOT NULL,
    pathogen_id character varying(20) NOT NULL,
    species_name character varying(150) NOT NULL,
    common_name character varying(150),
    abbreviation character varying(20),
    ncbi_refseq_accession character varying(30),
    biosafety_level smallint,
    human_pathogen boolean DEFAULT true,
    last_curated date,
    curator character varying(100),
    CONSTRAINT species_biosafety_level_check CHECK (((biosafety_level >= 1) AND (biosafety_level <= 4)))
);


--
-- Name: v_phenotype_surveillance; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.v_phenotype_surveillance AS
 SELECT pathogen_id,
    species_id,
    phenotype_category,
    evidence_strength,
    verification_status,
        CASE
            WHEN ("position" IS NOT NULL) THEN COALESCE((((((protein)::text || ':'::text) || (ref_aa)::text) || "position") || (alt_aa)::text), genotype_description)
            ELSE genotype_description
        END AS genotype_label,
    protein,
    "position",
    ref_aa,
    alt_aa,
    genotype_description,
    count(*) AS candidate_count
   FROM public.genotype_phenotype
  GROUP BY pathogen_id, species_id, phenotype_category, evidence_strength, verification_status, protein, "position", ref_aa, alt_aa, genotype_description;


--
-- Name: v_species_summary; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.v_species_summary AS
 SELECT pathogen_id,
    species_id,
    collection_year,
    collection_country,
    collection_country_code,
    count(*) AS total_genomes,
    min(collection_date) AS earliest_date,
    max(collection_date) AS latest_date
   FROM public.genome_metadata
  GROUP BY pathogen_id, species_id, collection_year, collection_country, collection_country_code;


--
-- Name: v_variant_summary; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.v_variant_summary AS
 SELECT pv.pathogen_id,
    pv.species_id,
    pv.gene,
    pv."position",
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    pv.is_stop,
    pv.genome_count,
    pv.species_total_genomes,
    round(((pv.genome_count)::numeric / (NULLIF(pv.species_total_genomes, 0))::numeric), 4) AS global_frequency,
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
   FROM (public.protein_variants pv
     LEFT JOIN public.gene_function gf ON ((((pv.species_id)::text = (gf.species_id)::text) AND ((pv.gene)::text = (gf.gene)::text))));


--
-- Name: v_variant_with_phenotype; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.v_variant_with_phenotype AS
 SELECT pv.pathogen_id,
    pv.species_id,
    pv.gene,
    pv."position",
    pv.ref_aa,
    pv.alt_aa,
    pv.variant_type,
    pv.hgvs_p,
    pv.genome_count,
    pv.species_total_genomes,
    round(((pv.genome_count)::numeric / (NULLIF(pv.species_total_genomes, 0))::numeric), 4) AS global_frequency,
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
   FROM (public.protein_variants pv
     LEFT JOIN public.genotype_phenotype gp ON ((((pv.pathogen_id)::text = (gp.pathogen_id)::text) AND ((pv.species_id)::text = (gp.species_id)::text) AND ((((pv.gene)::text = (gp.protein)::text) AND (pv."position" = gp."position") AND ((pv.alt_aa)::text = (gp.alt_aa)::text)) OR ((gp.lineage_id IS NOT NULL) AND ((gp.lineage_id)::text = ANY (pv.lineage_ids))) OR ((gp."position" IS NULL) AND (gp.protein IS NOT NULL) AND ((gp.protein)::text = (pv.gene)::text))))));


--
-- Name: gene_function gene_function_id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_function ALTER COLUMN gene_function_id SET DEFAULT nextval('public.gene_function_gene_function_id_seq'::regclass);


--
-- Name: protein_variants variant_id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_variants ALTER COLUMN variant_id SET DEFAULT nextval('public.protein_variants_variant_id_seq'::regclass);


--
-- Name: gene_function gene_function_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_function
    ADD CONSTRAINT gene_function_pkey PRIMARY KEY (gene_function_id);


--
-- Name: gene_function gene_function_species_id_gene_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_function
    ADD CONSTRAINT gene_function_species_id_gene_key UNIQUE (species_id, gene);


--
-- Name: genome_metadata genome_metadata_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_metadata
    ADD CONSTRAINT genome_metadata_pkey PRIMARY KEY (genome_accession);


--
-- Name: genotype_phenotype genotype_phenotype_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genotype_phenotype
    ADD CONSTRAINT genotype_phenotype_pkey PRIMARY KEY (association_id);


--
-- Name: lineages lineages_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.lineages
    ADD CONSTRAINT lineages_pkey PRIMARY KEY (lineage_id);


--
-- Name: pathogens pathogens_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.pathogens
    ADD CONSTRAINT pathogens_pkey PRIMARY KEY (pathogen_id);


--
-- Name: protein_variants protein_variants_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.protein_variants
    ADD CONSTRAINT protein_variants_pkey PRIMARY KEY (variant_id);


--
-- Name: reference_genomes reference_genomes_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.reference_genomes
    ADD CONSTRAINT reference_genomes_pkey PRIMARY KEY (accession);


--
-- Name: reference_proteomes reference_proteomes_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.reference_proteomes
    ADD CONSTRAINT reference_proteomes_pkey PRIMARY KEY (reference_accession, gene);


--
-- Name: species species_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.species
    ADD CONSTRAINT species_pkey PRIMARY KEY (species_id);


--
-- Name: idx_genefunc_conserved_gin; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_conserved_gin ON public.gene_function USING gin (conserved_regions);


--
-- Name: idx_genefunc_domains_gin; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_domains_gin ON public.gene_function USING gin (key_domains);


--
-- Name: idx_genefunc_fts; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_fts ON public.gene_function USING gin (protein_function_tsv);


--
-- Name: idx_genefunc_full_scope; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_full_scope ON public.gene_function USING btree (pathogen_id, species_id, gene);


--
-- Name: idx_genefunc_gene; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_gene ON public.gene_function USING btree (gene);


--
-- Name: idx_genefunc_hotspots_gin; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_hotspots_gin ON public.gene_function USING gin (known_hotspots);


--
-- Name: idx_genefunc_literature_gin; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_literature_gin ON public.gene_function USING gin (literature_refs);


--
-- Name: idx_genefunc_pathogen; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_pathogen ON public.gene_function USING btree (pathogen_id);


--
-- Name: idx_genefunc_pdb_ids_gin; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_pdb_ids_gin ON public.gene_function USING gin (pdb_ids);


--
-- Name: idx_genefunc_sites_gin; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_sites_gin ON public.gene_function USING gin (functional_sites);


--
-- Name: idx_genefunc_species; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genefunc_species ON public.gene_function USING btree (species_id);


--
-- Name: idx_genomemeta_country; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_country ON public.genome_metadata USING btree (collection_country);


--
-- Name: idx_genomemeta_date; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_date ON public.genome_metadata USING btree (collection_date);


--
-- Name: idx_genomemeta_isolate; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_isolate ON public.genome_metadata USING btree (isolate text_pattern_ops);


--
-- Name: idx_genomemeta_lineage; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_lineage ON public.genome_metadata USING btree (pathogen_id, species_id, lineage_id) WHERE (lineage_id IS NOT NULL);


--
-- Name: idx_genomemeta_lineage_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_lineage_id ON public.genome_metadata USING btree (lineage_id);


--
-- Name: idx_genomemeta_pathogen; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_pathogen ON public.genome_metadata USING btree (pathogen_id);


--
-- Name: idx_genomemeta_pathogen_country; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_pathogen_country ON public.genome_metadata USING btree (pathogen_id, collection_country);


--
-- Name: idx_genomemeta_pathogen_date; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_pathogen_date ON public.genome_metadata USING btree (pathogen_id, collection_date);


--
-- Name: idx_genomemeta_scope; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_scope ON public.genome_metadata USING btree (pathogen_id, species_id);


--
-- Name: idx_genomemeta_species; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_species ON public.genome_metadata USING btree (species_id);


--
-- Name: idx_genomemeta_strain; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_strain ON public.genome_metadata USING btree (strain text_pattern_ops);


--
-- Name: idx_genomemeta_temporal_geo; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_genomemeta_temporal_geo ON public.genome_metadata USING btree (species_id, collection_year, collection_country);


--
-- Name: idx_gp_category; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_category ON public.genotype_phenotype USING btree (phenotype_category);


--
-- Name: idx_gp_flagged; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_flagged ON public.genotype_phenotype USING btree (species_id, protein, "position") WHERE (record_flagged = true);


--
-- Name: idx_gp_full_scope; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_full_scope ON public.genotype_phenotype USING btree (pathogen_id, species_id, protein, "position");


--
-- Name: idx_gp_hgvs; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_hgvs ON public.genotype_phenotype USING btree (protein, "position", alt_aa);


--
-- Name: idx_gp_high_evidence; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_high_evidence ON public.genotype_phenotype USING btree (species_id, protein, "position") WHERE (evidence_strength = ANY (ARRAY['strong'::public.evidence_strength_enum, 'moderate'::public.evidence_strength_enum]));


--
-- Name: idx_gp_lineage_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_lineage_id ON public.genotype_phenotype USING btree (lineage_id);


--
-- Name: idx_gp_lineage_match; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_lineage_match ON public.genotype_phenotype USING btree (pathogen_id, species_id, lineage_id) WHERE (lineage_id IS NOT NULL);


--
-- Name: idx_gp_literature_gin; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_literature_gin ON public.genotype_phenotype USING gin (literature_refs);


--
-- Name: idx_gp_mutation_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_mutation_id ON public.genotype_phenotype USING btree (mutation_id);


--
-- Name: idx_gp_mutation_match; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_mutation_match ON public.genotype_phenotype USING btree (pathogen_id, species_id, protein, "position", alt_aa);


--
-- Name: idx_gp_pathogen; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_pathogen ON public.genotype_phenotype USING btree (pathogen_id);


--
-- Name: idx_gp_pathogen_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_pathogen_id ON public.genotype_phenotype USING btree (pathogen_id);


--
-- Name: idx_gp_scope_position; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_scope_position ON public.genotype_phenotype USING btree (species_id, protein, "position");


--
-- Name: idx_gp_species; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_species ON public.genotype_phenotype USING btree (species_id);


--
-- Name: idx_gp_species_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_species_id ON public.genotype_phenotype USING btree (species_id);


--
-- Name: idx_gp_strength; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_strength ON public.genotype_phenotype USING btree (evidence_strength);


--
-- Name: idx_gp_verif; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gp_verif ON public.genotype_phenotype USING btree (verification_status);


--
-- Name: idx_lineages_aliases; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_lineages_aliases ON public.lineages USING gin (known_aliases);


--
-- Name: idx_lineages_countries; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_lineages_countries ON public.lineages USING gin (countries_reported);


--
-- Name: idx_lineages_pathogen; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_lineages_pathogen ON public.lineages USING btree (pathogen_id);


--
-- Name: idx_lineages_pathogen_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_lineages_pathogen_id ON public.lineages USING btree (pathogen_id);


--
-- Name: idx_lineages_species; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_lineages_species ON public.lineages USING btree (species_id);


--
-- Name: idx_lineages_species_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_lineages_species_id ON public.lineages USING btree (species_id);


--
-- Name: idx_lineages_verif; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_lineages_verif ON public.lineages USING btree (verification_status);


--
-- Name: idx_proteinvars_frequency; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_proteinvars_frequency ON public.protein_variants USING btree (pathogen_id, species_id, genome_count DESC);


--
-- Name: idx_proteinvars_hgvs; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_proteinvars_hgvs ON public.protein_variants USING btree (hgvs_p);


--
-- Name: idx_proteinvars_scope; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_proteinvars_scope ON public.protein_variants USING btree (pathogen_id, species_id, gene);


--
-- Name: idx_proteinvars_variant_lookup; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_proteinvars_variant_lookup ON public.protein_variants USING btree (species_id, gene, "position", alt_aa);


--
-- Name: idx_refgenomes_pathogen; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_refgenomes_pathogen ON public.reference_genomes USING btree (pathogen_id);


--
-- Name: idx_refgenomes_species; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_refgenomes_species ON public.reference_genomes USING btree (species_id);


--
-- Name: idx_refgenomes_verif; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_refgenomes_verif ON public.reference_genomes USING btree (verification_status);


--
-- Name: idx_refproteomes_pathogen; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_refproteomes_pathogen ON public.reference_proteomes USING btree (pathogen_id);


--
-- Name: idx_refproteomes_scope; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_refproteomes_scope ON public.reference_proteomes USING btree (pathogen_id, species_id, gene);


--
-- Name: idx_refproteomes_species; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_refproteomes_species ON public.reference_proteomes USING btree (species_id);


--
-- Name: idx_species_pathogen; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_species_pathogen ON public.species USING btree (pathogen_id);


--
-- Name: uq_proteinvars_variant; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX uq_proteinvars_variant ON public.protein_variants USING btree (pathogen_id, species_id, gene, "position", ref_aa, alt_aa, variant_type);


--
-- Name: gene_function gene_function_pathogen_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_function
    ADD CONSTRAINT gene_function_pathogen_id_fkey FOREIGN KEY (pathogen_id) REFERENCES public.pathogens(pathogen_id);


--
-- Name: gene_function gene_function_species_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_function
    ADD CONSTRAINT gene_function_species_id_fkey FOREIGN KEY (species_id) REFERENCES public.species(species_id);


--
-- Name: genome_metadata genome_metadata_lineage_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_metadata
    ADD CONSTRAINT genome_metadata_lineage_id_fkey FOREIGN KEY (lineage_id) REFERENCES public.lineages(lineage_id) ON DELETE SET NULL;


--
-- Name: genome_metadata genome_metadata_pathogen_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_metadata
    ADD CONSTRAINT genome_metadata_pathogen_id_fkey FOREIGN KEY (pathogen_id) REFERENCES public.pathogens(pathogen_id);


--
-- Name: genome_metadata genome_metadata_reference_accession_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_metadata
    ADD CONSTRAINT genome_metadata_reference_accession_fkey FOREIGN KEY (reference_accession) REFERENCES public.reference_genomes(accession);


--
-- Name: genome_metadata genome_metadata_species_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genome_metadata
    ADD CONSTRAINT genome_metadata_species_id_fkey FOREIGN KEY (species_id) REFERENCES public.species(species_id);


--
-- Name: genotype_phenotype genotype_phenotype_lineage_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genotype_phenotype
    ADD CONSTRAINT genotype_phenotype_lineage_id_fkey FOREIGN KEY (lineage_id) REFERENCES public.lineages(lineage_id);


--
-- Name: genotype_phenotype genotype_phenotype_pathogen_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genotype_phenotype
    ADD CONSTRAINT genotype_phenotype_pathogen_id_fkey FOREIGN KEY (pathogen_id) REFERENCES public.pathogens(pathogen_id);


--
-- Name: genotype_phenotype genotype_phenotype_species_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genotype_phenotype
    ADD CONSTRAINT genotype_phenotype_species_id_fkey FOREIGN KEY (species_id) REFERENCES public.species(species_id);


--
-- Name: lineages lineages_pathogen_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.lineages
    ADD CONSTRAINT lineages_pathogen_id_fkey FOREIGN KEY (pathogen_id) REFERENCES public.pathogens(pathogen_id);


--
-- Name: lineages lineages_species_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.lineages
    ADD CONSTRAINT lineages_species_id_fkey FOREIGN KEY (species_id) REFERENCES public.species(species_id);


--
-- Name: reference_genomes reference_genomes_pathogen_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.reference_genomes
    ADD CONSTRAINT reference_genomes_pathogen_id_fkey FOREIGN KEY (pathogen_id) REFERENCES public.pathogens(pathogen_id);


--
-- Name: reference_genomes reference_genomes_species_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.reference_genomes
    ADD CONSTRAINT reference_genomes_species_id_fkey FOREIGN KEY (species_id) REFERENCES public.species(species_id);


--
-- Name: reference_proteomes reference_proteomes_pathogen_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.reference_proteomes
    ADD CONSTRAINT reference_proteomes_pathogen_id_fkey FOREIGN KEY (pathogen_id) REFERENCES public.pathogens(pathogen_id);


--
-- Name: reference_proteomes reference_proteomes_reference_accession_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.reference_proteomes
    ADD CONSTRAINT reference_proteomes_reference_accession_fkey FOREIGN KEY (reference_accession) REFERENCES public.reference_genomes(accession);


--
-- Name: reference_proteomes reference_proteomes_species_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.reference_proteomes
    ADD CONSTRAINT reference_proteomes_species_id_fkey FOREIGN KEY (species_id) REFERENCES public.species(species_id);


--
-- Name: species species_pathogen_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.species
    ADD CONSTRAINT species_pathogen_id_fkey FOREIGN KEY (pathogen_id) REFERENCES public.pathogens(pathogen_id);


--
-- Name: SCHEMA public; Type: ACL; Schema: -; Owner: -
--

REVOKE USAGE ON SCHEMA public FROM PUBLIC;


--
-- PostgreSQL database dump complete
--

\unrestrict pzG6puWa0aiNZRlVp6flk0RhHcWllMbYxAEhdVTRGB4mwZ8eRah6ZG8DfHcljhR

