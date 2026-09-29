/*
 * Local module: MUTATION_SUMMARY
 *
 * LLM-assisted summaries for the mutation profile at pipeline time:
 *   - mutation_summary.json: plain-language "mutation landscape" note
 *     (query-vs-background separation, elevated proteins, phenotype-linked
 *     mutations, query-enriched positions)
 *   - protein_summaries.tsv: one digest per protein condensed from the
 *     UniProt function/domain text attached to mutated positions.
 * The dashboard renders these files — no live LLM call at view time.
 */

process MUTATION_SUMMARY {
    tag "$meta.id - $meta.species"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_evidence.yml"

    input:
    tuple val(meta), path(mutation_profile_dir, stageAs: "mutation_profile")

    output:
    tuple val(meta), path("mutation_summary.json"),  emit: summary
    tuple val(meta), path("protein_summaries.tsv"),  emit: protein_summaries

    when:
    task.ext.when == null || task.ext.when

    script:
    def ollama_host  = task.ext.ollama_host  ?: params.ollama_host  ?: 'http://localhost:11434'
    def ollama_model = task.ext.ollama_model ?: params.ollama_model ?: ''
    def ollama_temp  = task.ext.ollama_temperature ?: params.ollama_temperature ?: '0.1'
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\$PATH"

    python3 ${projectDir}/bin/summarize_mutation_profile.py \
        --profile-dir mutation_profile \
        --species "${meta.species}" \
        --outdir . \
        --ollama-host "${ollama_host}" \
        --ollama-model "${ollama_model}" \
        --temperature ${ollama_temp}
    """
}
