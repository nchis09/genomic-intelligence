# pgirl/genomic-intelligence

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

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Provide a consensus FASTA file and a metadata TSV file:

- `sequences.fasta` — one or more consensus sequences (multiple pathogen species can be mixed; the pipeline assigns species using Nextclade).
- `metadata.tsv` — sample metadata with at least `strain`, `date`, and `country` columns.

Now, you can run the pipeline using:

```bash
nextflow run nchis09/genomic-intelligence \
   -profile <docker/singularity/.../institute> \
   --fasta sequences.fasta \
   --metadata metadata.tsv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

pgirl/genomic-intelligence was originally written by GOARN fellowship.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use pgirl/genomic-intelligence for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
