# Genomic Epidemic Intelligence System

A knowledge-based system that translates pathogen genomic data into epidemic intelligence and early warning alerts by answering epidemiological questions through curated evidence databases.

## Overview

Unlike existing bioinformatics pipelines (Nextflow, Genome Detective, etc.) that focus on genomic analysis, this system focuses on **epidemic intelligence** - converting genomic findings into actionable public health insights by connecting them to historical data, published evidence, and epidemiological context.

### The Problem

Traditional genomic tools tell you:
- Lineage: A
- Mutation: GP:A82V
- Location: Kasese District
- Date: July 2027

**But they don't answer:**
- Should we respond?
- Is this concerning?
- What does this mean for public health?

### Our Solution

**Evidence-based intelligence, not biological interpretation.**

Instead of asking "What does this mutation do biologically?", we ask:
- "Has this mutation been reported before?"
- "Is there published evidence of public health relevance?"
- "What have epidemiologists observed about this variant?"

## Target Pathogens

- **Ebola virus** (EBOV)
- **Rift Valley fever virus** (RVFV)
- **Influenza virus**
- **Dengue virus** (DENV)
- **Mpox virus** (MPXV)

## Architecture: Library + Engine

The system has two halves:

1. **Reference Library (PGIRL)** — curated *facts* and curated *rules*. The knowledge core.
2. **Intelligence Engine** — code that subjects a new genome to the library, links external
   databases, and assembles the Brief. It uses deterministic lookup by default — **no AI/ML required**. An optional local Ollama LLM can be enabled for narrative synthesis, but the deterministic path works without it.

```
NEW GENOME ──► Intelligence Engine ──► queries ──► Reference Library (PGIRL)
                     │                                     │
                     ├── links external DBs (NCBI/WHO/…)   │
                     └────────────► INTELLIGENCE BRIEF ◄───┘
```

### Reference Library tiers
- **Tier 1 — Curated Facts:** taxonomy, reference genomes, gene function, molecular
  epidemiology, mutations, genotype-phenotype associations, outbreaks.
- **Tier 2 — Controlled Vocabularies:** shared `_vocabularies/` (phenotypes, evidence
  strength, risk tiers). Interpretation rules are engine logic in `genomic_intelligence_engine/`,
  not curated data. *No AI — deterministic lookup.*

## Project Structure

```
genomic_intelligence_system/
├── README.md                           # This file
├── workflow.readme.readme              # Detailed workflow documentation
├── reference_library/                  # PGIRL — curated knowledge core
│   └── <pathogen>/                     # ebola | dengue | influenza | rvf | mpox
│       ├── taxonomy/                   # Tier 1
│       ├── reference_genomes/          # Tier 1
│       ├── gene_function/              # Tier 1
│       ├── molecular_epidemiology/     # Tier 1  (lineages/clades/serotypes)
│       ├── mutations/                  # Tier 1  (mutation catalogue)
│       ├── genotype_phenotype/         # Tier 1  (mutation → phenotype associations)
│       └── outbreaks/                  # Tier 1
├── database/                           # DB schema + shared controlled vocabularies
│   ├── schema.sql                      # Single source of truth for current schema
│   ├── _vocabularies/                  # Tier 2 — shared controlled terms
│   │   ├── phenotypes.yaml
│   │   ├── evidence_strength.yaml
│   │   └── risk_tiers.yaml
│   ├── _schemas/                       # YAML schema templates (one per library type)
│   │   ├── taxonomy.schema.yaml
│   │   ├── reference_genome.schema.yaml
│   │   ├── gene_function.schema.yaml
│   │   ├── molecular_epidemiology.schema.yaml
│   │   ├── mutation.schema.yaml
│   │   └── outbreak.schema.yaml
│   ├── views/                          # Surveillance views for the intelligence engine
│   │   └── surveillance_views.sql
│   └── ebola/                          # Pathogen-specific DB pipelines
│       └── protein_variants/
│           ├── fetch_reference_proteomes.py
│           ├── call_variants.py
│           ├── assign_lineages.py        # lineage assignment for DB genomes
│           └── variant_catalogue.sql
├── intelligence_engine/                # Tier 3 — code (no AI/ML)
│   ├── README.md
│   ├── bioinformatics/                 # new genome → species/lineage/mutations
│   │   ├── qc/
│   │   │   └── stage0_quality_control.py
│   │   ├── taxonomic_classification/
│   │   │   └── stage1_classification.py
│   │   ├── nextclade_runner.py         # shared Nextclade wrapper
│   │   ├── templates/                  # staged bioinformatics output templates
│   │   └── validation/                 # metadata schema definitions
│   ├── data_engine/                    # query reference_library + external sources
│   │   ├── sql_querying/               # PostgreSQL curated queries
│   │   │   └── bioinformatics_query.py
│   │   ├── online_querying/            # WHO, CDC, ProMED, HealthMap, NCBI, etc.
│   │   ├── llm_querying/               # structured LLM-based summarisation
│   │   └── epi_questions.json
│   ├── evidence_integration/           # harmonize evidence + cross-evidence statistics
│   │   ├── engine.py                   # core assessment logic
│   │   ├── harmonization.py            # unified EvidenceObjects
│   │   ├── cross_evidence.py           # statistical cross-evidence analyzers
│   │   ├── visualization.py            # evidence-network / geo / temporal charts
│   │   ├── analyzers/
│   │   ├── pipeline/
│   │   ├── tree/
│   │   ├── figures/                    # R figure scripts for the report
│   │   └── examples/
│   └── genomic_intelligence/           # LLM evidence synthesis (no new stats, no recommendations)
├── scripts/                            # Entry points
└── config/
    └── parameters.yaml
```

## Inputs

A submission is a **multi-FASTA + a CSV/TSV metadata sheet**, linked by `sample_id`.
Required per sample: `sample_id`, `collection_date`, `country`, `admin1`, `host`
(location granularity is critical for the "seen here before" logic). Full spec, field
tiers, and a template are in [`input/README.md`](./input/README.md) and
[`intelligence_engine/bioinformatics/modules/schemas/metadata_schema.yaml`](./intelligence_engine/bioinformatics/modules/schemas/metadata_schema.yaml). Anything derivable from the
FASTA (completeness, lineage, mutations, closest reference) is **computed by the engine**,
not submitted.

## Quick start

### Requirements

- [Nextflow](https://www.nextflow.io/) (`>=23.04`)
- Python 3.10+ with the packages listed in `environment.yml`
  - The wrapper scripts auto-detect and prefer Anaconda/Miniconda if installed.
  - Alternatively, use `-profile conda` to let Nextflow build the `pgirl` Conda environment.
- PostgreSQL 16
- [Ollama](https://ollama.com/) (optional; only needed for `--use_llm true`)

### 1. Set up the environment and the reference database

This installs required Python packages, checks external tools, creates the `pgirl` PostgreSQL database, loads the schema, and fetches the canonical reference proteomes.

> **Recommended first run / CI setup:** the curated-data sync downloads data from NCBI, UniProt and PubTator and can take **10-30 minutes**. To load only the schema and canonical reference proteomes (enough to run the deterministic pipeline), use:
>
> ```bash
> PGIRL_SYNC_SOURCES=none ./setup.sh
> ```

Full setup with curated data:

```bash
./setup.sh
```

If you prefer a managed Conda environment, use:

```bash
./setup.sh -profile conda
./run_pipeline.sh -profile conda ...
```

### 2. Run the end-to-end analysis pipeline

#### Fully deterministic (recommended first run; no LLM required)

```bash
./run_pipeline.sh \
  --input_fasta input/input_FASTA.fasta \
  --input_metadata input/input_metadata.csv \
  --outdir output \
  --use_llm false
```

#### With local LLM narrative synthesis

Pull the default local model first (7B — fast on most laptops):

```bash
ollama pull qwen2.5:7b
```

For higher quality narrative at the cost of speed, use the 14B model:

```bash
OLLAMA_MODEL=qwen2.5:14b ./run_pipeline.sh ... --use_llm true
```

Then run with the LLM enabled:

```bash
./run_pipeline.sh \
  --input_fasta input/input_FASTA.fasta \
  --input_metadata input/input_metadata.csv \
  --outdir output \
  --use_llm true
```

`run_pipeline.sh` starts Ollama if it is not already running and stops the server it started once the pipeline finishes. If Ollama is unavailable or an LLM call times out, the pipeline automatically falls back to deterministic output.

All defaults use `db_url=postgresql://localhost:5432/pgirl`, so `--db_url` can be omitted if your database is at that URL.

### Output organisation

After a successful run, the key results are consolidated under `${outdir}/reports/<sample_id>/`:

- `final_report.txt` — combined genomic intelligence brief + report
- `figures/` — copied figures referenced by the report
- `context_used.json` — grounding context passed to the LLM (for traceability)

A machine-readable summary of the whole run is written to:

- `${outdir}/run_summary.json`

Intermediate stage outputs remain under `${outdir}/bioinformatics/`, `${outdir}/data_query/`, `${outdir}/evidence_integration/`, and `${outdir}/genomic_intelligence/`.

### Logs and work directories

`run_pipeline.sh` keeps the project root clean:

- Nextflow log: `${outdir}/.nextflow.log`
- Execution work directory: `${outdir}/work/`

### Common parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_fasta` | `input/input_FASTA.fasta` | Multi-sample input FASTA |
| `--input_metadata` | `input/input_metadata.csv` | Sample metadata CSV |
| `--outdir` | `output` | Root output directory |
| `--db_url` | `postgresql://localhost:5432/pgirl` | PostgreSQL connection URL |
| `--pathogen` | `ebola` | Target pathogen for setup and analysis |
| `--use_llm` | `false` | Enable local Ollama LLM narrative synthesis |

## Ebola Variant Data Pipeline

Per-species, reference-based protein variant calling for all 6 Ebola species:

```bash
# 1. Fetch reference genomes + proteomes (one-time, or when references change)
/Users/christianndekezi/anaconda3/bin/python3 database/ebola/protein_variants/fetch_reference_proteomes.py

# 2. Run the full pipeline (all 6 species)
/Users/christianndekezi/anaconda3/bin/python3 -u database/ebola/protein_variants/call_variants.py --batch-size 50

# 3. Incremental update — only new NCBI genomes
/Users/christianndekezi/anaconda3/bin/python3 -u database/ebola/protein_variants/call_variants.py --batch-size 50 --skip-existing
```

The pipeline stores only **genome metadata** and **variant calls** in PostgreSQL — full
sequences are never persisted.

## Pathogen-Specific Intelligence Modules

### Ebola Intelligence Module

**Epidemiological Questions:**
1. Has this Ebola species/lineage been detected in this region before?
2. Is this a known outbreak strain or novel introduction?
3. What is the closest historical outbreak? When? Where?
4. Are there mutations in GP, VP24, VP35, NP with known epidemiological significance?
5. Has this specific mutation profile been associated with:
   - Increased transmission?
   - Immune escape?
   - Vaccine efficacy changes?
   - Diagnostic escape?
6. Is there evidence of spillover vs sustained human-to-human transmission?
7. What is the confidence level in the intelligence assessment?

**Evidence Categories:**
- Mutation epidemiology (published associations)
- Lineage geographic distribution
- Historical outbreak patterns
- Transmission cluster signatures
- Phenotypic associations (from literature)

### Dengue Intelligence Module

**Epidemiological Questions:**
1. Has this DENV serotype/genotype been detected in this region before?
2. Is this consistent with local endemic circulation or new introduction?
3. Are there mutations in E, NS1, NS3, NS5 with known epidemiological significance?
4. Has this mutation been associated with:
   - Increased transmission?
   - Severe disease (DSS/DHF)?
   - Vaccine escape?
   - Diagnostic failure?
5. Is there evidence of serotype shift?
6. What is the outbreak risk based on historical patterns?

### Influenza Intelligence Module

**Epidemiological Questions:**
1. HA/NA subtype - has this combination been seen in this region?
2. Segment typing - evidence of reassortment?
3. Host adaptation markers - zoonotic potential?
4. Antigenic mutations - vaccine match assessment?
5. Antiviral resistance markers - treatment implications?
6. Pandemic lineage assessment - historical pandemic genotypes?

### Rift Valley Fever Intelligence Module

**Epidemiological Questions:**
1. L/M/S segment characterization - known pathogenicity profiles?
2. Reassortment evidence - novel combinations?
3. Livestock vs human linkage - source attribution?
4. Vector association - known vector competence?
5. Geographic lineage - expansion patterns?

### Mpox Intelligence Module

**Epidemiological Questions:**
1. Clade assignment - known differences in severity/transmission?
2. APOBEC signatures - human adaptation evidence?
3. Gene gain/loss - phenotypic implications?
4. Human adaptation markers - transmission efficiency?
5. Geographic expansion - new territories?

## The Knowledge Base

### Mutation Evidence Database

Each mutation entry includes:
```yaml
mutation: GP:A82V
pathogen: Ebola virus
gene: GP
position: 82
reference: A
alternate: V

epidemiological_associations:
  - phenotype: increased_infectivity
    evidence: Science (2016)
    confidence: moderate
    citation: "Diehl et al. Science 2016;351(6277):1178-1183"
  
  - phenotype: immune_escape
    evidence: limited
    confidence: low
    citation: "No direct evidence"

public_health_relevance: yes
last_curated: 2027-01-15
curator: [expert name]
notes: "Associated with 2014-2016 West Africa outbreak"
```

### Lineage Intelligence Database

Each lineage entry includes:
```yaml
lineage: Sudan virus
pathogen: Ebola virus
geographic_distribution: [South Sudan, Uganda, DRC]
first_detected: 1976
outbreak_history:
  - year: 2000
    location: Uganda
    cases: 425
    cfr: 53%
  - year: 2012
    location: Uganda
    cases: 17
    cfr: 41%
transmission_patterns: sporadic_outbreaks
known_epidemiological_features: []
```

### Outbreak Database

Each outbreak entry includes:
```yaml
outbreak_id: EBOV-UGA-2022
pathogen: Ebola virus
species: Sudan virus
lineage: A
location: Uganda
start_date: 2022-09-20
end_date: 2023-01-11
total_cases: 142
deaths: 55
cfr: 39%
source: probable_spillover
genomic_data_available: yes
key_mutations: [GP:A82V, VP35:T230A]
publications: []
```

## Intelligence Output Format

The system produces a two-part output per sample (or batch):

1. **Part 1 — Genomic Intelligence Brief**: A concise one-page executive snapshot
   that synthesises the key findings into a format readable in under two minutes.
   Designed for public health specialists and emergency operations centres that
   need rapid situational awareness to make immediate decisions.

2. **Part 2 — Genomic Intelligence Report**: The comprehensive evidence document
   with full detail across eight sections. Designed for analysts and experts who
   need the complete evidence chain to conduct deeper assessment.

Both parts are generated together; the brief links to the full report by ID.

### Part 1: Genomic Intelligence Brief

```
================================================================================
                         GENOMIC INTELLIGENCE BRIEF
================================================================================
Brief ID:
Report ID:
Generated:     
Sample ID:     
================================================================================

PATHOGEN IDENTIFICATION
--------------------------------------------------------------------------------
Species:     
Lineage:    
Confidence:  
Closest match:

SAMPLE CONTEXT
--------------------------------------------------------------------------------
Location:        
Collection date: 
Source:          
Genome quality:  

Executive Assessment
--------------------------------------------------------------------------------




Genome Integrity and Analytical Confidence
--------------------------------------------------------------------------------
  



Identity of the Virus
--------------------------------------------------------------------------------
  



Genomic Characteristics
--------------------------------------------------------------------------------
  



Integrated Evidence Assessment
--------------------------------------------------------------------------------



Evolutionary Context
--------------------------------------------------------------------------------



Historical and Epidemiological Context
--------------------------------------------------------------------------------


Integrated Evidence Assessment
--------------------------------------------------------------------------------



Sources of Uncertainty
--------------------------------------------------------------------------------



Overall Genomic Intelligence Assessment
--------------------------------------------------------------------------------


================================================================================
  This brief is auto-generated. It MUST be reviewed and approved by a
  qualified public health genomicist before dissemination or action.
  Full evidence chain: see Genomic Intelligence Report number; ...
================================================================================
```

### Part 2: Genomic Intelligence Report

The comprehensive report provides the full evidence chain behind every
finding in the brief. It has eight sections, each with explicit confidence
levels and evidence citations.

```
================================================================================
                      GENOMIC INTELLIGENCE REPORT
================================================================================
Report ID:     
Brief ID:      
Generated:     
Sample ID:     
Review status: PENDING_EXPERT_REVIEW
================================================================================

SECTION 1: SAMPLE & SURVEILLANCE CONTEXT
--------------------------------------------------------------------------------
Sample ID:              
Sampling source:        
Collection date:        
Country:                
Admin Level 1:         
Admin Level 2:          
Host:                   
Host status:            
Sequencing platform:    
Sequencing protocol:    
Submitting lab:         
Genome completeness:    
Mean coverage depth:  
Genome quality flag:    

SECTION 2: GENOMIC CHARACTERIZATION
--------------------------------------------------------------------------------
Family:                 
Genus:                  
Species:                
Lineage:                
Genotype / serotype:    
Clade:                  
Closest reference:      
  - Identity:           
  - SNPs from reference: 
Closest outbreak genome: 
  - SNPs from closest:   
  - Evolutionary distance:
Phylogenetic placement: 
  - Support:            
  - Tree file:          
Genome length:          
Segmented:              
Annotation source:      
Genes covered:          
Missing regions:       

Confidence: 

SECTION 3: FUNCTIONAL GENOMICS — MUTATIONS & GENOMIC FEATURES
--------------------------------------------------------------------------------
Total SNPs vs closest reference: 
Total amino acid substitutions: 
Insertions / deletions:          

3a. MUTATIONS WITH CURATED EVIDENCE
...................................
  Mutation 1: 
  Type:     
  Effect:   
  Evidence (Ref publications):
    -
    - 
    - 
  Confidence: STRONG
  Public health relevance:
    - Transmissibility:  
    - Virulence:         
    - Diagnostics:       
    - Therapeutics:      
    - Vaccine:           
  Note: 

   Mutation 2: 
  Type:     
  Effect:   
  Evidence (Ref publications):
    -
    - 
    - 
  Confidence: STRONG
  Public health relevance:
    - Transmissibility:  
    - Virulence:         
    - Diagnostics:       
    - Therapeutics:      
    - Vaccine:           
  Note: 

   Mutation n: 
  Type:     
  Effect:   
  Evidence (Ref publications):
    -
    - 
    - 
  Confidence: STRONG
  Public health relevance:
    - Transmissibility:  
    - Virulence:         
    - Diagnostics:       
    - Therapeutics:      
    - Vaccine:           
  Note: 


3b. NOVEL / UNCHARACTERISED MUTATIONS
......................................
  Mutation: 
  Type:     synonymous / non-coding
  Evidence: 
  Confidence: 
  Public health relevance: 



SECTION 4: MOLECULAR EPIDEMIOLOGY
--------------------------------------------------------------------------------
4a. HISTORICAL OCCURRENCE
..........................
  Species xxx (the one detected in this sample ) first described:      ... (country)
  Total confirmed outbreaks:         xx (1976, 1979, 2000-2001, 2004, 2011,2012, 2014, 2022 these years should also be repace and in each year we put the actual outbreak number )
  Cumulative cases (all outbreaks):  ~
  Cumulative CFR (all outbreaks):    ~%

4b. LINEAGE DISTRIBUTION
.........................
  Lineage A:
    - First detected:    .., country and area if available =
    - affected counties so far Countries:   list them 
    - Last detected:     year , countiry
    - Associated outbreaks:
        * SUDV-SSD-1976  (Nzara, 284 cases, 53% CFR)
        * SUDV-UGA-2000  (Gulu, 425 cases, 53% CFR)
        * SUDV-UGA-2022  (Mubende, 142 cases, 39% CFR)
    - Cross-border spread documented: Yes (South Sudan ↔ Uganda, multiple events)
    (those infomation should be removed as well they are just helping with the reporting style)

4c. GEOGRAPHIC CONTEXT
.......................
  Endemic in:               name counties if any
  The current region or district:
    - Previous  detection:    None (last filovirus: EBOV 2012 Kibaale)
    - Proximity to previous outbreak: ~ km from name county or region =
    - Cross-border proximity:    country (≤ xxx km)

4d. RESERVOIR & HOST RANGE
............................
  Known reservoirs:          Fruit bats (Hypsignathus monstrosus, Epomops franqueti,
                              Myonycteris torquata) — suspected, not definitive
  Accidental hosts:          Humans, non-human primates (chimpanzees, gorillas)
  Vectors:                   None (direct contact transmission)

4e. TEMPORAL TRENDS
....................
  Inter-outbreak interval:   2-8 years (species-level)
  Current detection gap:     5 years since last SUDV in Uganda (2022)
  Assessment:                Detection within expected inter-outbreak interval;
                              not an anomalous re-emergence interval

4f. INTRODUCTION / SPREAD ASSESSMENT
.....................................
  Novel introduction to district:    YES/no
    - SUDV not previously detected in Kasese District
    - Closest to previous outbreak of this kind: Mubende (2022), 280 km
  Cross-border spread risk:          MODERATE
    - Kasese borders DRC; cross-border population movement is common
  Evidence of epidemiological link:  POSSIBLE
    - 3 SNPs from 2022 Mubende strain
    - Temporal gap (5 years) is plausible for persistent transmission chain or
      re-emergence from reservoir
  Confidence: MODERATE

SECTION 5: EVIDENCE SUMMARY & KNOWLEDGE GAPS
--------------------------------------------------------------------------------
5a. EVIDENCE INVENTORY
.......................
  Peer-reviewed literature cited:    3
  International guidance cited:      1 (WHO SUDV risk assessment)
  Reference databases consulted:     NCBI RefSeq, OutbreakDB (curated)
  Curated mutation records matched:  1 (GP:A82V)
  Novel mutations requiring review:  2 (VP35:T230A, intergenic C→T)

5b. CONFIDENCE SUMMARY
.......................
  Species identification:            HIGH
  Lineage assignment:                 HIGH
  Phylogenetic placement:             MODERATE
  Mutation effect interpretation:     STRONG (GP:A82V) / INSUFFICIENT (VP35:T230A)
  Molecular epidemiology context:     MODERATE
  Introduction / spread assessment:   MODERATE
  Overall report confidence:          MODERATE

5c. KNOWLEDGE GAPS & UNCERTAINTIES
...................................
  - VP35:T230A has no curated functional evidence; wet-lab characterisation needed
  - Single genome available; cannot assess within-outbreak diversity or transmission
    clusters
  - 3' terminal region incomplete (285 bp); cannot rule out 3' UTR variants
  - Reservoir source of this introduction unknown (no bat sampling data)
  - No contemporaneous genomes from neighbouring districts for comparison
  - SUDV vaccine efficacy data limited (no licensed SUDV vaccine as of report date)

SECTION 6: GENOMIC INTELLIGENCE ASSESSMENT
--------------------------------------------------------------------------------
6a. WHAT IS KNOWN
..................
  - The sample is Sudan ebolavirus, lineage A, clade A1
  - It is 3 SNPs from the 2022 Mubende outbreak strain
  - It carries GP:A82V, a well-characterised mutation associated with increased
    transmissibility (strong evidence)
  - It carries VP35:T230A, a novel substitution in a functional domain with no
    curated evidence
  - SUDV has caused 8 outbreaks since 1976 with ~54% cumulative CFR
  - This is the first SUDV detection in Kasese District

6b. WHAT HAS CHANGED vs PREVIOUS SURVEILLANCE
..............................................
  - SUDV not detected in Uganda since 2022 (5-year gap)
  - First detection in Kasese District (geographic expansion)
  - GP:A82V present (also present in 2022 strain — conserved, not new)
  - VP35:T230A is novel — not seen in any prior SUDV genome
  - No evidence of recombination or major genomic rearrangement

6c. PUBLIC HEALTH SIGNIFICANCE
...............................
  - SUDV is a Biosafety Level 4 pathogen with high CFR (~54% historically)
  - No licensed vaccine exists for SUDV (unlike EBOV)
  - GP:A82V may increase transmission efficiency
  - Novel VP35 mutation may alter innate immune evasion — significance unknown
  - Detection in a new district near the DRC border raises cross-border risk
  - Overall: MODERATE public health significance
    (known lineage, known mutation, but new location + novel VP35 change)



SECTION 7: EVIDENCE-BASED CONSIDERATIONS FOR DECISION-MAKING
--------------------------------------------------------------------------------
  This section presents the evidence and factors relevant to public health
  decision-making. It does NOT prescribe specific actions — the analyst
  reviewing this report is responsible for determining appropriate responses
  based on local context, available resources, and applicable regulations.

7a. SURVEILLANCE CONSIDERATIONS
................................
  Evidence:
    - Only 1 genome available from this event — cannot assess diversity or
      transmission clusters
    - VP35:T230A is novel; its frequency in subsequent cases will determine
      whether it is an isolated event or an emerging variant
    - Neighbouring districts (Bundibugyo, Kabarole, Ntoroko) also border DRC
      but have no contemporaneous genomic data
    - No genomic data available from DRC side of the border for comparison
    - GP:A82V was present in the 2022 Mubende strain — tracking its persistence
      informs understanding of SUDV evolution between outbreaks
  Factors for the analyst to weigh:
    - How many additional cases can be sequenced, and how quickly?
    - Is there genomic surveillance infrastructure in neighbouring districts?
    - What is the current sequencing turnaround time in-country?
    - Are there data-sharing agreements with DRC for cross-border comparison?

7b. SOURCE INVESTIGATION CONSIDERATIONS
........................................
  Evidence:
    - 3 SNPs from 2022 Mubende strain — close genetic relationship
    - 5-year temporal gap between 2022 and this detection
    - Two plausible scenarios cannot be distinguished with current data:
      a) Re-emergence from a reservoir (new spillover from bats)
      b) Continued cryptic circulation from 2022 (persistent transmission chain)
    - Kasese has no prior SUDV history; nearest filovirus event was EBOV in
      Kibaale (2012, ~100 km away)
    - No bat sampling data from Kasese District
    - Ebola virus is known to persist in immune-privileged sites (eye, CNS,
      semen) in survivors — relevant if scenario (b) is considered
  Factors for the analyst to weigh:
    - Is there epidemiological evidence of travel links to Mubende?
    - Are there known survivors from the 2022 outbreak in or near Kasese?
    - What is the local bushmeat / bat contact exposure profile?
    - Can animal reservoir sampling be initiated in the area?

7c. LABORATORY CONSIDERATIONS
..............................
  Evidence:
    - Species identification confidence is HIGH (97.8% identity to reference,
      7/7 genes covered) — but this is the first SUDV detection in a new
      district, so independent confirmation strengthens the assessment
    - VP35:T230A is a novel call in a functional domain — Sanger confirmation
      would rule out sequencing artifact
    - 3' terminal 285 bp has insufficient coverage — UTR variants cannot be
      excluded in this region
    - Genome was sequenced on MinION (higher error rate than Illumina);
      homopolymer regions may need validation
  Factors for the analyst to weigh:
    - Is independent confirmation feasible and at what turnaround?
    - Can additional sequencing resolve the 3' terminal region?
    - Should raw reads be deposited for independent verification?
    - If VP35:T230A is confirmed, is functional characterisation warranted?

7d. EXPERT REVIEW REQUIREMENTS
................................
  The following findings require expert judgement beyond what the curated
  knowledge base can provide:

  - VP35:T230A:
      Requires: Virologist / immunologist
      Question:  Does this substitution in the IFN antagonist domain alter
                 innate immune evasion capacity?
      Why:       VP35 is a known IFN antagonist; changes in this domain could
                 affect virulence, but no functional data exists for T230A

  - Phylogenetic placement:
      Requires: Phylogenetics expert
      Question:  Is the clade A1 assignment robust? Can recombination be
                 ruled out?
      Why:       Single genome with moderate bootstrap support; recombination
                 analysis requires multiple sequences and specialized methods

  - Report sign-off:
      Requires: Public health genomicist
      Question:  Are the findings and confidence levels appropriate for
                 dissemination?
      Why:       This is the first SUDV detection in a new district with a
                 novel mutation — the assessment must be validated before
                 it informs public health decisions

SECTION 8: VISUALISATIONS
--------------------------------------------------------------------------------
  [Figure 1] Phylogenetic tree — SUDV genomes with this sample highlighted
             File: gir-2027-0715-ebov-001.fig1.tree.png
             Type: Maximum likelihood, UFboot support, time-scaled

  [Figure 2] Mutation map — genome diagram with all SNPs annotated
             File: gir-2027-0715-ebov-001.fig2.mutation_map.png
             Type: Linear genome with gene tracks + variant positions

  [Figure 3] Genome annotation — gene structure and coverage plot
             File: gir-2027-0715-ebov-001.fig3.annotation.png
             Type: Gene tracks + per-position coverage depth

  [Figure 4] Geographic distribution — map of SUDV outbreaks + this detection
             File: gir-2027-0715-ebov-001.fig4.geo_map.png
             Type: Choropleth + point markers, 1976-present

  [Figure 5] Temporal trend — SUDV detections over time with lineage colours
             File: gir-2027-0715-ebov-001.fig5.temporal.png
             Type: Timeline / epidemic curve with lineage annotation

================================================================================
                         END OF REPORT
================================================================================
```

### Report Design Principles

1. **Every claim has a confidence level** — HIGH / MODERATE / LOW / INSUFFICIENT
2. **Every mutation cites its evidence** — peer-reviewed literature, guidance, or "no evidence available"
3. **Negative findings are reported** — absence of known concerning variants is itself actionable intelligence
4. **Knowledge gaps are explicit** — the report states what is NOT known, not just what is known
5. **Evidence is presented, not prescriptions** — the system provides the evidence base; the analyst decides what actions to take
6. **Escalation triggers are defined** — quantitative criteria for upgrading risk classification based on subsequent data
7. **Expert review is mandatory** — automated report is always flagged for human sign-off before dissemination
8. **Visualisations are referenced** — each figure is linked to a generated file for inclusion in briefings
9. **Brief precedes report** — the executive snapshot enables rapid decision-making; the full report provides the evidence chain

## Risk Classification Framework

- **ROUTINE** - Known patterns, no concerning features
- **MONITOR** - Some unusual features, enhanced surveillance warranted
- **INVESTIGATE** - Novel or concerning features requiring field investigation
- **HIGH PRIORITY** - Significant public health concern, urgent response
- **EMERGENCY** - Immediate public health action required

## Development Approach

1. **Literature Curation** - Systematically review literature for each pathogen
2. **Knowledge Base Development** - Build structured evidence databases
3. **Question Framework Design** - Define epidemiological questions per pathogen
4. **Intelligence Engine** - Build evidence lookup and assessment logic
5. **Validation** - Test against historical outbreaks
6. **Iterative Refinement** - Continuous knowledge base updates

## Key Principles

1. **Evidence over interpretation** - Report what is known, not what we think
2. **Confidence levels** - Always indicate certainty of assessments
3. **Expert-in-the-loop** - Flag for expert review when evidence is limited
4. **Transparency** - Show the evidence chain for every assessment
5. **Curated knowledge** - Quality-controlled evidence, not automated mining

## Differentiation from Existing Tools

| Feature | Traditional Tools | This System |
|---------|------------------|-------------|
| Focus | Genomic analysis | Epidemiological intelligence |
| Output | Mutations, phylogeny | Public health recommendations |
| Mutation interpretation | Biological prediction | Evidence-based associations |
| Context | Limited | Historical outbreak data |
| Questions answered | What is it? | What should we do? |

## Future Enhancements

- Integration with real-time surveillance systems
- Automated literature updates
- Geographic information system integration
- Multi-pathogen outbreak correlation
- API integration with public health databases
- monitor mutations and see how they are evolving
- we can think of showing the mutation on the protein structure if possible
- also put a like a seq logo for that particular mutaiton to help view its prefrency vs other protein frequencies
- also think of having a part which can allow local integration to the local public health epidemiological infomation



- This should help make the risk assessment more accurate and reliable
- Read the guideline for making a risk assessment 
- Views I would add to support this:
    - I can create SQL views that the engine (and analysts) can query directly:

    - v_mutation_surveillance — per-mutation counts by species, country, year
    - v_mutation_geography — first seen, last seen, countries where observed
    - v_mutation_trends — yearly counts to detect emergence/increase
    - v_genome_variant_profile — all variants per genome, ready for comparison
    - v_mutation_with_phenotype — join variants to known phenotype associations

### Missing pieces we should add
To make this fully useful, we still need:

Lineage/clade assignment for each genome 


