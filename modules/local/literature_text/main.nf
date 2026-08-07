/*
 * Local module: LITERATURE_TEXT
 *
 * Extracts plain text from downloaded PDFs using pymupdf or pdfplumber.
 */

process LITERATURE_TEXT {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_text.yml"

    input:
    tuple val(meta), path("*.pdf")

    output:
    tuple val(meta), path("*.txt"), optional: true, emit: text
    path "pdf_text_summary.json", emit: summary

    when:
    !params.skip_literature_grobid && (task.ext.when == null || task.ext.when)

    script:
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    python ${projectDir}/bin/extract_literature_pdf_text.py \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --format txt
    """
}
