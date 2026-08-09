/*
 * Local module: LITERATURE_EVIDENCE
 *
 * Extracts structured, quotable, confidence-scored evidence from full-text
 * .txt files using MedSpaCy (spaCy + clinical NLP) with rule-based matching,
 * ConText negation detection, and section-aware confidence scoring.
 * Merges PubMed metadata (title, year, authors, journal, DOI) into output JSON.
 */

process LITERATURE_EVIDENCE {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_evidence.yml"

    input:
    tuple val(meta), path("*.txt"), path("*.json")
    path rules_yml

    output:
    tuple val(meta), path("evidence/*.json"), optional: true, emit: evidence

    when:
    !params.skip_literature_evidence && !params.skip_literature_text && (task.ext.when == null || task.ext.when)

    script:
    def metadata_arg = params.skip_pubmed_metadata ? "" : "--metadata-dir metadata"
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    # --- Stage metadata JSONs into metadata/ subdirectory ---
    mkdir -p metadata
    mv *.json metadata/ 2>/dev/null || true

    # --- Run MedSpaCy evidence extraction ---
    # Output `*.json` files are written to evidence/ and published by
    # the publishDir block in conf/modules.config.
    set -e
    mkdir -p evidence
    python ${projectDir}/bin/extract_evidence_medspacy.py \
        --input-dir . \
        --outdir evidence/ \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --rules-yml "${rules_yml}" \
        ${metadata_arg}
    """
}
