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

## Introduction

**pgirl/genomic-intelligence** is a Nextflow DSL2 pipeline for multi-pathogen genomic epidemic intelligence. It takes consensus FASTA sequences and sample metadata as input, performs pathogen classification using Nextclade, builds phylogenetic trees using the Nextstrain framework, and produces annotated tree visualisations with amino-acid mutation heatmaps.

### Pipeline steps

1. **Classification** — Nextclade dataset download and sequence classification ([`Nextclade`](https://github.com/nextstrain/nextclade))
2. **Phylogenetics** — Nextstrain/augur phylogenetic tree building per species
3. **Annotation** — Extract tip metadata and top amino-acid mutation matrices from Nextstrain JSONs
4. **Visualisation** — Annotated rectangular phylogeny with branches colored by country and mutation heatmap ([`ggtree`](https://bioconductor.org/packages/ggtree/))
5. **QC report** — Aggregate quality metrics ([`MultiQC`](http://multiqc.info/))

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

Then apply the augur compatibility patch (renames `_resolve_filepath` → `resolve_filepath`):

```bash
PATCH_FILE="data/nextstrain_ebola/shared/vendored/snakemake/config.smk"
if grep -q "_resolve_filepath" "$PATCH_FILE" 2>/dev/null; then
  sed -i'' -e 's/_resolve_filepath/resolve_filepath/g' "$PATCH_FILE"
  echo "Patch applied ✓"
fi
```

## Usage

> [!NOTE]
> If you are new to Nextflow, refer to the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html) to set up the runtime.

Provide a consensus FASTA file and a metadata TSV file:

- `sequences.fasta` — one or more consensus sequences (multiple pathogen species can be mixed; the pipeline assigns species using Nextclade).
- `metadata.tsv` — sample metadata with at least `strain`, `date`, and `country` columns.

Now, you can run the pipeline using:

```bash
nextflow run main.nf \
  --fasta input/input_FASTA.fasta \
  --metadata input/metadata.tsv \
  --outdir results \
  -profile conda \
  -resume
```

Run this command from the root directory of a local clone of the repository. The `conda` profile provides the required software environment, and `-resume` reuses successfully completed tasks from a previous run when possible.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Outputs

For each detected species the pipeline produces:

- `results/nextstrain_ebola/{species}/auspice/*.json` — interactive Nextstrain/Auspice dataset.
- `results/nextstrain_ebola/{species}/results/` — Nextstrain Augur results including `tree.nwk`.
- `results/nextstrain_ebola/annotations/{species}*_tip_metadata.tsv` and `*_mutation_matrix.tsv` — extracted tip metadata and mutation matrix.
- `results/figures/{species}_tree_heatmap.png` — static annotated phylogeny with per-genome mutation heatmap.
- `results/figures/{species}_geo_map.png` — static world map showing the geographic distribution of samples.

Additionally, `results/pipeline_info/pipeline_metro_map_*.html` — an auto-generated [nf-metro](https://github.com/seqeralabs/nf-metro) metro-map diagram of the run's actual Nextflow task graph (skipped with a warning if `nf-metro` is unavailable; see Prerequisites).

## Credits

pgirl/genomic-intelligence was originally written in collaboration with the Uganda Virus Research Institute (UVRI), Robert Koch Institute (RKI), Global Outbreak Alert and Response Network (GOARN), and the WHO Hub for Pandemic and Epidemic Intelligence.


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- If you use pgirl/genomic-intelligence for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->


An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

