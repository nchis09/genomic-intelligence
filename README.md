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
6. **Phenotype annotation** — discover UniProt accessions for the query samples' proteins and annotate them with [`UniprotR`](https://github.com/Proteomicslab57357/UniprotR), `UniProtExtractR`, [`rbioapi`](https://cran.r-project.org/package=rbioapi), and Pfam HMM scans ([`HMMER`](http://hmmer.org/))
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

The phylogenetics stage needs local background sequences and metadata from the Nextstrain ingest workflow. These are gitignored in the upstream repo and must be generated locally before the main pipeline runs. Run the ingest Snakemake for each species you plan to analyse (replace `<species>` with `bdbv`, `sudv`, `ebov`, or a list like `[bdbv,sudv]`):

```bash
cd data/nextstrain_ebola/ingest
snakemake --snakefile Snakefile \
  --cores 4 \
  --config species=[<species>] \
  --rerun-incomplete \
  --nolock \
  data/<species>/sequences.fasta data/<species>/metadata.tsv
cd ../../..
```

For example, to prepare the BDBV and SUDV background data:

```bash
cd data/nextstrain_ebola/ingest
snakemake --snakefile Snakefile \
  --cores 4 \
  --config species=[bdbv,sudv] \
  --rerun-incomplete \
  --nolock \
  data/bdbv/sequences.fasta data/bdbv/metadata.tsv \
  data/sudv/sequences.fasta data/sudv/metadata.tsv
cd ../../..
```

## Usage

> [!NOTE]
> If you are new to Nextflow, refer to the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html) to set up the runtime.

Provide a consensus FASTA file and a metadata TSV file:

- `sequences.fasta` — one or more consensus sequences (multiple pathogen species can be mixed; the pipeline assigns species using Nextclade).
- `metadata.tsv` — sample metadata with at least `strain`, `date`, and `country` columns.

The literature stages query NCBI Entrez, so a contact email is required via `--pubmed_email` (use `--pdf_email` to set a different contact for the Unpaywall PDF resolver).

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
| `results/phenotype_annotation/{pathogen}_{species}/` | Accession discovery tables, query protein FASTA/mutations, and UniprotR / UniProtExtractR / rbioapi / HMM annotation results. |
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
| `annotate_uniprotr.R` / `annotate_uniprotextractr.R` / `annotate_rbioapi.R` | Annotate the discovered proteins with function, GO, pathway, and interaction data. |
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

