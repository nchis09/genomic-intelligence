/*
 * Local module: NEXTSTRAIN_EBOLA
 *
 * Run the Nextstrain Ebola Snakemake workflow which:
 *   1. Downloads/uses curated background sequences from ingest/data/{species}/
 *   2. Injects user query sequences as additional_inputs
 *   3. Runs the full augur pipeline (filter, align, tree, refine, ancestral, traits, export)
 *
 * Conda env is managed by Nextflow (no manual activation needed).
 */

process NEXTSTRAIN_EBOLA {
    tag "$meta.id"
    label 'process_high'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"

    input:
    tuple val(meta), path(fasta), path(metadata)

    output:
    tuple val(meta), path("${meta.species}/auspice/*.json")   , emit: auspice
    tuple val(meta), path("${meta.species}/results/")         , emit: results_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def species        = meta.species ?: meta.id.replaceAll(/^.*_/, '')
    def builds_list    = task.ext.builds ?: "${species}/all-outbreaks"
    def nextstrain_dir = "${projectDir}/data/nextstrain_ebola"
    def configfile     = task.ext.configfile ?: ''
    """
    NEXTSTRAIN_DIR="${nextstrain_dir}"
    INGEST_DIR="\${NEXTSTRAIN_DIR}/ingest/data/${species}"

    # Create a config override that:
    # 1. Uses LOCAL ingest data (not S3) as background sequences
    # 2. Injects user query sequences as additional_inputs
    cat > pgirl_override.yaml <<EOF
inputs:
  - name: local_background
    species: ${species}
    sequences: \${INGEST_DIR}/sequences.fasta
    metadata: \${INGEST_DIR}/metadata.tsv

additional_inputs:
  - name: pgirl_query
    species: ${species}
    sequences: \$(realpath ${fasta})
    metadata: \$(realpath ${metadata})

builds:
  - ${builds_list}
EOF

    PHYLO_DIR="\${NEXTSTRAIN_DIR}/phylogenetic"

    mkdir -p ${species}

    snakemake \\
        --snakefile \${PHYLO_DIR}/Snakefile \\
        --cores ${task.cpus} \\
        --directory ${species} \\
        --configfile pgirl_override.yaml \\
        ${configfile ? "--configfile ${configfile}" : ''} \\
        --rerun-incomplete \\
        --nolock
    """
}
