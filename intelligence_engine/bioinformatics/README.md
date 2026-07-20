# PGIRL Bioinformatics Pipeline

This directory contains the modular, pathogen-aware bioinformatics workflow for
the Genomic Epidemic Intelligence System. The pipeline accepts assembled
consensus genomes (FASTA) and produces a normalised `bio_output.json` for each
sample, which is consumed by the downstream data engine, evidence integration
and genomic intelligence stages.

## Architecture

```
intelligence_engine/bioinformatics/
├── pipeline.py                     # Top-level orchestrator
├── nextclade_runner.py             # Nextclade CLI wrapper
└── modules/
    ├── quality_control.py          # Stage 0: input validation & QC
    ├── taxonomic_classification.py # Stage 1: pathogen identification
    ├── schemas/                    # Metadata & taxonomy schemas
    │   ├── metadata_schema.yaml
    │   └── ebola_taxonomy.yaml
    ├── reference_selection.py      # Stage 2: curated reference + context
    ├── alignment.py                # Stage 3: alignment (planned enhancement)
    ├── variant_calling.py          # Stage 4: variant calling (Nextclade reuse)
    ├── lineage_assignment.py       # Stage 5: clade, lineage, QC & mutations
    ├── phylogenetics.py            # Stage 6: ML tree (Nextclade stub)
    ├── recombination.py            # Stage 7: recombination / reassortment
    ├── comparative_genomics.py     # Stage 8: gene content / dN/dS
    └── normalisation.py            # Stage 9: bio_output.json assembler
└── pathogen_workflows/
    └── ebola.py                    # Ebola-specific stage wiring
```

## Quick start

Run the full Ebola pipeline on the bundled example data:

```bash
/Users/christianndekezi/anaconda3/bin/python3 -m \
  intelligence_engine.bioinformatics.pipeline \
  --fasta input/input_FASTA.fasta \
  --metadata input/input_metadata.csv \
  --pathogen ebola \
  --species-id EBOV \
  --output-dir output/bioinformatics
```

Supported Ebola species identifiers: `EBOV`, `SUDV`, `BDBV`, `RESTV`, `TAFV`, `BOMV`.

## Pipeline stages

| Stage | Module | Status | Notes |
|-------|--------|--------|-------|
| 0 | `modules/quality_control.py` | Implemented | metadata validation, seqkit stats, genome QC flags |
| 1 | `modules/taxonomic_classification.py` | Implemented | Kraken2 + NCBI BLAST; Ebola taxonomy bundled; optional in Ebola workflow |
| 2 | `modules/reference_selection.py` | Implemented | queries PGIRL DB, fetches NCBI reference, gathers context |
| 3 | `modules/alignment.py` | Stub | will add MAFFT standalone alignment |
| 4 | `modules/variant_calling.py` | Implemented (Nextclade-backed) | amino-acid variants primary; nucleotide changes as supporting evidence |
| 5 | `modules/lineage_assignment.py` | Implemented | Nextclade clade/QC/mutation parser |
| 6 | `modules/phylogenetics.py` | Stub | currently copies Nextclade placement tree |
| 7 | `modules/recombination.py` | Stub | RDP5 / reassortment planned |
| 8 | `modules/comparative_genomics.py` | Partial | GC content; dN/dS planned |
| 9 | `modules/normalisation.py` | Implemented | assembles `bio_output.json` |

## Downstream integration

The `bio_output.json` files feed directly into:

```bash
# Data engine / DB queries
python3 -m intelligence_engine.data_engine.sql_querying.bioinformatics_query \
  --bioinformatics-dir output/bioinformatics \
  --output-dir output/data_query

# Epidemiological queries (requires a running Ollama/LLM backend)
python3 -m intelligence_engine.data_engine.online_querying.epi_query_engine \
  --bio-output output/bioinformatics/EBOV-UGA-2027-001/bio_output.json \
  --db-query-results output/data_query/EBOV-UGA-2027-001/db_query_results.json \
  --output output/data_query/EBOV-UGA-2027-001/epi_output.json
```

## Implementation notes

- Nextclade is used as the validated tool for Ebola lineage assignment,
  mutation calling and QC. Module stubs are in place so each stage can later
  be swapped with a standalone implementation (MAFFT + IQ-TREE2 + custom
  variant caller) for unsupported pathogens.
- The reference selection module queries the PGIRL PostgreSQL database
  (`pgirl`) for curated reference genomes and proteomes; make sure the DB
  is running on `postgresql://localhost:5432/pgirl` or set `--db-url`.
- Pathogen-aware dispatch happens at `pipeline.py`; additional pathogens can
  be added by implementing a module under `pathogen_workflows/` and wiring it
  in the orchestrator.
