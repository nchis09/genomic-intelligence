# Input Specification

A submission = **two files**, linked by `sample_id`:

1. **Sequences** — a multi-FASTA (`*.fasta` / `*.fa`) containing one or many consensus genomes.
2. **Metadata** — a CSV or TSV sheet, one row per genome.

We assume the user has already done everything up to producing the consensus FASTA
(basecalling, assembly, consensus). Anything derivable from the FASTA is **computed by the
engine**, not submitted.

## 1. FASTA rules
- One record per genome. The record header (text after `>`, up to the first whitespace)
  **must equal the `sample_id`** in the metadata sheet.
- Segmented pathogens (RVF = L/M/S; influenza = 8 segments): use headers
  `<sample_id>|<segment>`, e.g. `EBOV-UGA-2027-001` (single) or `RVFV-KEN-2027-004|L`.
- IUPAC nucleotides allowed; `N` for ambiguous. The engine reports % ambiguous itself.

```
>EBOV-UGA-2027-001
ATGGGT...consensus sequence...
>EBOV-UGA-2027-002
ATGGGT...consensus sequence...
```

## 2. Metadata rules
- CSV (comma) or TSV (tab). UTF-8. First row = header with the field names below.
- One row per `sample_id`. Every FASTA record must have a matching row and vice-versa.
- **Empty cell = "not provided / unknown".** Controlled fields must use exact allowed values.
- Full field definitions + allowed values: [`metadata_schema.yaml`](../intelligence_engine/bioinformatics/modules/schemas/metadata_schema.yaml)
- Copy [`metadata_template.csv`](./metadata_template.csv) to start (delete the EXAMPLE rows).

### Field tiers (summary)
| Tier | Fields | Behaviour if missing |
|------|--------|----------------------|
| **Required** | `sample_id`, `collection_date`, `country`, `admin1`, `host` | Sample rejected |
| **Recommended** | `admin2`, `host_species`*, `travel_history`, `travel_locations`, `sample_type`, `symptom_onset_date`, `vaccination_status` | Runs, but flags reduced confidence / unanswered questions |
| **Optional** | `vaccine_details`, `latitude`, `longitude`, `outcome`, `suspected_exposure`, `epi_link_id`, `existing_accession`, `notes` | Used if present |

\* `host_species` is **required when `host = animal`**.

### Why these fields
The genome answers *what it is*; the metadata answers the epidemiology:
- `country` + `admin1`/`admin2` + `collection_date` → **"has this been seen here before / novel introduction?"**
- `host` (+ `host_species`) → **spillover vs human-to-human**
- `travel_history` → **importation vs local transmission**
- `vaccination_status` → interpret **vaccine-escape** mutations meaningfully

**Geographic granularity is make-or-break.** `country` alone is too coarse to answer
"seen here before" — always provide at least `admin1`, ideally `admin2` (district).

## 3. Computed by the engine (do NOT submit)
`genome_length`, `percent_ambiguous_bases`, `genome_completeness`, `species`,
`lineage_or_clade`, `called_mutations`, `closest_reference_genome`.

## 4. Validation (planned)
On ingest the engine will check: FASTA↔metadata `sample_id` match, required fields present,
date formats, controlled-vocabulary values, and `travel_locations` present when
`travel_history = yes`. Invalid rows are reported, not silently dropped.
