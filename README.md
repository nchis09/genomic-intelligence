# pgirl/genomic-intelligence

> [!WARNING]
> This pipeline is under active development. Its workflow, parameters, and outputs may change before a stable release.

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/pgirl/genomic-intelligence)
[![GitHub Actions CI Status](https://github.com/pgirl/genomic-intelligence/actions/workflows/nf-test.yml/badge.svg)](https://github.com/pgirl/genomic-intelligence/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/pgirl/genomic-intelligence/actions/workflows/linting.yml/badge.svg)](https://github.com/pgirl/genomic-intelligence/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.3-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.3)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/pgirl/genomic-intelligence)

````text
  ____  ___  _____
 / ___||_ _||  ___|
| |  _  | | | |_
| |_| | | | |  _|
 \____||___||_|


  o===o       o===o       o===o       o===o
 /     \     /     \     /     \     /     \
o       o---o       o---o       o---o       o
 \     /     \     /     \     /     \     /
  o===o       o===o       o===o       o===o
````

## Introduction

**pgirl/genomic-intelligence** is a Nextflow DSL2 pipeline for multi-pathogen genomic epidemic intelligence. It takes consensus FASTA sequences and sample metadata as input, identifies the pathogen and species with Nextclade, and then assembles genomic, phenotypic, literature, and epidemiological evidence for those samples into a single PostgreSQL **knowledge warehouse** that can be queried for risk assessment.

The pipeline is organised as a **pathogen router**: samples are grouped by the species Nextclade assigns, and each group is dispatched to a pathogen-specific workflow. Ebola (`orthoebolavirus`) is currently the only registered pathogen; species groups without a registered workflow are reported in an `unsupported` summary and skipped with a warning.

### Pipeline steps

1. **Classification** — download every configured Nextclade dataset, run each sample against all of them, and assign the pathogen/species with the best QC score ([`Nextclade`](https://github.com/nextstrain/nextclade))
2. **Pathogen routing** — split samples into per-species groups and dispatch each group to its pathogen workflow
3. **Bioinformatics** — Nextstrain/Augur build per species ([`nextstrain/ebola`](https://github.com/nextstrain/ebola)), plus a model-aware maximum-likelihood tree from the subsampled sequences ([`MAFFT`](https://mafft.cbrc.jp/alignment/software/) + [`IQ-TREE 2`](http://www.iqtree.org/))
4. **Epidemiological data** — search and download matching disease datasets from the Humanitarian Data Exchange (`rhdx`)
5. **Literature retrieval** — Europe PMC search per species and evidence domain, PubMed metadata fetch, deduplication, [`ASReview`](https://asreview.nl/) title/abstract screening, open-access PDF download, PDF-to-text conversion, rule-based structured evidence extraction, and evidence QC
6. **Phenotype annotation** — discover UniProt accessions for the query samples' proteins and annotate them with `UniProtExtractR`, [`rbioapi`](https://cran.r-project.org/package=rbioapi), and Pfam HMM scans ([`HMMER`](http://hmmer.org/))
7. **Knowledge warehouse** — start a shared PostgreSQL instance, ingest every species' outputs into the schema defined by `database/knowledge_schema.sql`, then stop the server

Most stages can be turned off individually (for example `--skip_literature_search`, `--skip_phenotype_annotation`, `--skip_hmm_annotation`, `--skip_iqtree`, `--skip_epi_data`, `--skip_knowledge_warehouse`); see `nextflow.config` for the full parameter list.

## Prerequisites

- **Conda** or **Mamba** (for environment management)
- **Java 11–24** (required by Nextflow; Java 25 is not yet supported)
- **Git** (to clone nextstrain/ebola during setup)
- **[nf-metro](https://github.com/seqeralabs/nf-metro)** (optional) — auto-generates a metro-map diagram of each run's task graph. If `nf-metro` isn't on `PATH`, the pipeline lazily creates a dedicated conda env from `envs/pgirl_nf_metro.yml` on first use (requires `conda`); if that also fails, the diagram is skipped with a warning and the pipeline continues normally.

### 1. Clone the repository

```bash
git clone https://github.com/nchis09/genomic-intelligence.git
cd genomic-intelligence
```

### 2. Set up the Nextstrain Ebola data

The phylogenetics stage relies on background sequences from the [nextstrain/ebola](https://github.com/nextstrain/ebola) repository. These files are **not** included in the repository and must be cloned locally:

```bash
git clone https://github.com/nextstrain/ebola.git data/nextstrain_ebola
```

The main pipeline automatically runs the Nextstrain Ebola ingest workflow for each detected species, so the required background `sequences.fasta` and `metadata.tsv` files are generated on demand. You can therefore proceed straight to **Usage**.

If you prefer to pre-generate the background data (e.g. to save time during repeated runs, to run offline, or to inspect the ingest outputs), run the ingest Snakemake manually:

```bash
cd data/nextstrain_ebola/ingest
snakemake --snakefile Snakefile \
  --cores 4 \
  --config species=[bdbv,sudv,ebov] \
  --rerun-incomplete \
  --nolock \
  data/bdbv/sequences.fasta data/bdbv/metadata.tsv \
  data/sudv/sequences.fasta data/sudv/metadata.tsv \
  data/ebov/sequences.fasta data/ebov/metadata.tsv
cd ../../..
```

To skip the auto-ingest step and use pre-generated files, pass `--skip_nextstrain_ingest true` to the main pipeline.

## Usage

> [!NOTE]
> If you are new to Nextflow, refer to the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html) to set up the runtime.

Provide a consensus FASTA file and a metadata TSV file:

- `sequences.fasta` — one or more consensus sequences (multiple pathogen species can be mixed; the pipeline assigns species using Nextclade).
- `metadata.tsv` — sample metadata with at least `strain`, `date`, and `country` columns.

### Segmented genomes (influenza)

Influenza isolates are submitted as **one FASTA record per segment**, all in the same `sequences.fasta` (segments from multiple isolates can be mixed freely). Because the pipeline screens each record independently against every dataset, the only requirement is a consistent naming convention:

```text
>{LABID}_{segment}
```

- `LABID` — your lab/sample identifier; everything before the last `_<segment>` token. All segments of one isolate share the same `LABID`.
- `segment` — any label: `seg4`, `seg6`, `HA`, `NA`, `pb2`, … This is **advisory only** — the pipeline verifies which segment each record actually is during Nextclade screening, so a mislabelled segment is still routed correctly.

```text
>SAMPLE001_seg4    ← HA segment of isolate SAMPLE001
>SAMPLE001_seg6    ← NA segment of isolate SAMPLE001
>SAMPLE002_seg4    ← HA of a different isolate
```

`metadata.tsv` uses the same format for every pathogen — **one row per sample/isolate** (keyed by `strain`, `sample_id`, `accession` or `name`), which is automatically fanned out to every `LABID_*` segment record. For non-segmented pathogens (Ebola, …) this is simply one row per genome:

```text
strain      date         country
SAMPLE001   2026-03-05   Uganda
SAMPLE002   2026-03-07   Uganda
UG_01       2026-01-01   DRC
```

If your records have real GenBank accessions, either use the accession as the record header/metadata id itself, or keep it in an `accession` column — when the pipeline needs to rewrite `accession` to the record id (segment fan-out), the original value is preserved as `genbank_accession`.

Header rules that will break tools if violated:

- No spaces in the id (text after the first space is ignored).
- No `()`, `,`, `:`, `;`, `[`, `]` or `/` — these break tree/Newick handling downstream.
- Every record id must be unique — never use the same header twice.

### How lineage assignment works

Screening is two-pass:

1. **Per record** — every record is scored against every dataset; `qc.overallScore` (lower is better) picks the best hit. Records whose best score exceeds `--nextclade_max_score` (default `100`; ~0-30 means a real same-lineage match) are treated as **untyped** — they cannot decide anything on their own, so contaminant or junk sequences cannot create phantom pathogen groups.
2. **Per isolate** — records sharing a `LABID` prefix form an isolate, and the isolate's lineage is decided by its **HA record first, then NA**. All other segments inherit that call (`lineage_source=inherited` in `species_assignments.tsv`).

Consequences:

- An isolate only enters a seasonal-flu build when it has a typed **HA or NA** record of a supported lineage (`h1n1pdm`, `h3n2`, `vic` — generic Flu B and Yamagata records are placed in the `vic` build). Submitting only internal segments (PB2, PB1, PA, NP, MP, NS) yields a screened, reported group marked internal-only — **no tree build is run** for it.
- `species_assignments.tsv` reports both levels: `species` = the isolate's final call, `record_species` = what the record itself typed as, plus `segment`, `lineage_source` (`ha`/`na`/`record`/`inherited`/`none`) and `reassortment_suspected` (flagged when a typed record disagrees with its isolate's consensus — except avian NA records, which can only hit seasonal NA datasets anyway).

The literature stages query NCBI Entrez, so a contact email is required via `--pubmed_email` (use `--pdf_email` to set a different contact for the Unpaywall PDF resolver).

Literature tuning flags (defaults live in `nextflow.config`; override with `--<flag>`):

| Flag | Default | What it does |
| --- | --- | --- |
| `--literature_min_year` | current − 5 | Earliest publication year searched. |
| `--literature_max_year` | current | Latest publication year searched. |
| `--literature_search_max_results` | `200` | Max Europe PMC hits per species × domain. |
| `--asreview_top_n` | `25` | Maximum papers per species × domain carried into PDF + LLM stages — the main volume knob. `asreview_n_stop` is auto-raised to ≥ `top_n` so the ranking can fill the requested slots. |
| `--asreview_min_year` | `literature_min_year` | Publication-year floor applied during screening. |
| `--asreview_n_prior_included` / `--asreview_n_prior_excluded` | `5` / `5` | Keyword-seeded ASReview priors. |
| `--skip_literature_search` / `--skip_pubmed_metadata` / `--skip_literature_deduplication` / `--skip_literature_screening` / `--skip_literature_pdf` / `--skip_literature_text` / `--skip_literature_evidence` / `--skip_evidence_qc` | `false` | Skip individual literature stages. |

Epidemiological data flags (WHO FluNet for influenza; HDX for Ebola):

|| Flag | Default | What it does |
|| --- | --- | --- |
|| `--epi_min_year` | current − 5 | Earliest ISO year for WHO FluNet surveillance fetch. |
|| `--epi_max_year` | current | Latest ISO year for WHO FluNet surveillance fetch. |
|| `--skip_epi_data` | `false` | Skip epidemiological data fetch. |

Now, you can run the pipeline using:

```bash
nextflow run main.nf \
  --fasta input/input_FASTA.fasta \
  --metadata input/metadata.tsv \
  --outdir results \
  --pubmed_email you@example.org \
  -profile conda \
  -resume
```

Run this command from the root directory of a local clone of the repository. The `conda` profile provides the required software environment, and `-resume` reuses successfully completed tasks from a previous run when possible.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Outputs

Results are published under `--outdir` (default `results/`), mostly one subdirectory per stage and per detected species:

| Directory | Contents |
| --- | --- |
| `results/nextclade/` | Per-sample Nextclade output for every screened dataset. |
| `results/classification/` | `species_assignments.tsv`, `species_groups.json`, and the per-species FASTA/metadata splits. |
| `results/route/` | Router summary, including any species groups with no registered pathogen workflow. |
| `results/nextstrain_ebola/{species}/` | Nextstrain/Augur build: `auspice/*.json` (interactive dataset) and `results/` (including `tree.nwk`). |
| `results/bioinformatics/{pathogen}_{species}/` | MAFFT alignment and IQ-TREE 2 maximum-likelihood tree. |
| `results/epidemiological_data/{species}/` | HDX search summary and downloaded epidemiological records. |
| `results/literature_retrieval/` | One subdirectory per stage: `literature_search`, `literature_metadata`, `literature_deduplicated`, `literature_screened`, `literature_pdfs`, `literature_text`, `literature_evidence`. |
| `results/evidence_qc/{species}/` | QC report plus `clean/` and `failed/` evidence JSON sets. |
| `results/phenotype_annotation/{pathogen}_{species}/` | Accession discovery tables, query protein FASTA/mutations, and UniProtExtractR / rbioapi / HMM annotation results. |
| `results/knowledge_warehouse/` | Per-species ingestion logs, the shared PostgreSQL data directory, and a SQL dump of the run's database. |

Additionally, `results/pipeline_info/pipeline_metro_map_*.html` — an auto-generated [nf-metro](https://github.com/seqeralabs/nf-metro) metro-map diagram of the run's actual Nextflow task graph (skipped with a warning if `nf-metro` is unavailable; see Prerequisites).

For a high-level conceptual overview of the pipeline's architecture (rather than the literal per-run task graph), see [`docs/architecture_overview.html`](docs/architecture_overview.html).

### Schema visualization

After the pipeline finishes, generate an interactive SchemaSpy report of the knowledge-warehouse database:

```bash
conda run -n pgirl_schemaspy python bin/run_schemaspy.py --outdir results
```

This writes `results/pipeline_info/schemaspy/index.html` and requires the SchemaSpy JAR and PostgreSQL JDBC driver in `assets/schemaspy/` (or set `SCHEMASPY_JAR` and `PGJDBC_JAR` environment variables).

### Knowledge warehouse & downstream analysis

The pipeline builds a PostgreSQL **knowledge warehouse** that links sample metadata, genomic features, phylogenetic trees, literature evidence, and epidemiological records. The warehouse is populated by `bin/build_knowledge_db.py` and is defined by `database/knowledge_schema.sql`.

Post-run visualisation:

- `bin/run_schemaspy.py` generates an interactive SchemaSpy ER report.
- DBeaver can be connected to the live database for an interactive ER diagram.

**Ongoing work:** We are extending this warehouse to support epidemiological queries for risk assessment, such as outbreak detection, transmission mapping, and mutation-phenotype associations.

### Key helper scripts in `bin/`

| Script | Purpose |
| --- | --- |
| `build_knowledge_db.py` | Build and populate the PostgreSQL knowledge warehouse. |
| `start_shared_db.py` / `stop_shared_db.py` | Start and stop the shared PostgreSQL server. |
| `run_schemaspy.py` | Generate a SchemaSpy HTML report of the warehouse schema. |
| `extract_query_proteins.py` | Discover UniProt accessions and extract query proteins/mutations for phenotype annotation. |
| `annotate_uniprotextractr.R` / `annotate_rbioapi.R` | Annotate the discovered proteins with function, GO, pathway, and interaction data. |
| `parse_hmmscan.py` | Parse `hmmscan` output into Pfam domain, sequence, and summary tables. |
| `literature_search.py` / `fetch_pubmed_metadata.py` | Search Europe PMC and fetch PubMed metadata per species and evidence domain. |
| `deduplicate_literature.R` / `screen_literature.R` | Deduplicate records and run ASReview title/abstract screening. |
| `fetch_literature_pdfs.py` / `extract_literature_pdf_text.py` | Resolve and download open-access PDFs, then extract their text. |
| `extract_literature_evidence.py` / `extract_evidence_medspacy.py` | Extract structured evidence from full text using rules and clinical NLP. |
| `run_evidence_qc.py` | Quality-control extracted evidence and emit clean/failed sets. |
| `fetch_rhdx.R` | Search and download epidemiological datasets from the Humanitarian Data Exchange. |

## Credits

pgirl/genomic-intelligence was originally written in collaboration with the Uganda Virus Research Institute (UVRI), Robert Koch Institute (RKI), Global Outbreak Alert and Response Network (GOARN), and the WHO Hub for Pandemic and Epidemic Intelligence.


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- If you use pgirl/genomic-intelligence for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->


An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

