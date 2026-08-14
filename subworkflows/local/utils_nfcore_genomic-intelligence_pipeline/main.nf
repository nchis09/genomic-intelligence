//
// Subworkflow with functionality specific to the pgirl/genomic-intelligence pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    // OR from direct --fasta / --metadata CLI arguments
    //

    if (params.fasta) {
        // Direct CLI mode: --fasta <file> --metadata <file>
        // Species will be auto-detected by Nextclade downstream
        def fasta_file = file(params.fasta, checkIfExists: true)
        def meta_file  = params.metadata ? file(params.metadata, checkIfExists: true) : file('NO_FILE')
        def sample_id  = fasta_file.baseName.replaceAll(/\.[^.]+$/, '')
        def meta_map   = [ id: sample_id ]

        channel
            .of([ meta_map, fasta_file, meta_file ])
            .set { ch_samplesheet }
    } else {
        channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fasta, metadata ->
                return [ meta, fasta, metadata ]
            }
            .set { ch_samplesheet }
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

        //
        // Print a consolidated summary of samples/pathogens that were skipped
        // because no workflow exists yet for their detected pathogen family
        // (see subworkflows/local/pathogen_router/main.nf).
        //
        def unsupported_file = file("${outdir}/pipeline_info/unsupported_pathogens.tsv")
        if (unsupported_file.exists()) {
            def lines = unsupported_file.readLines().drop(1)  // skip header
            if (lines) {
                log.warn ""
                log.warn "=============================================================="
                log.warn "  WARNING: ${lines.size()} sample(s) skipped — no workflow yet"
                log.warn "=============================================================="
                lines.each { line ->
                    def (sample_id, pathogen, species) = line.split('\t')
                    log.warn "  - ${sample_id} (pathogen: '${pathogen}', species: '${species}')"
                }
                log.warn "See ${unsupported_file} for full details."
                log.warn "=============================================================="
            }
        }

        //
        // Auto-generate a metro-map style diagram of this run's actual task
        // graph from the Mermaid DAG Nextflow just wrote (dag.file in
        // nextflow.config, *.mmd). Best-effort only: any failure here only
        // logs a warning and never fails the pipeline.
        //
        renderNfMetroDiagram(outdir)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    if (!params.input && !params.fasta) {
        error("Please provide either --input <samplesheet.csv> or --fasta <file.fasta> --metadata <file.tsv>")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    // No-op for FASTA+metadata input; schema validation handles required fields
    return input
}

//
// Best-effort: render the Mermaid DAG that Nextflow wrote for this run
// (nextflow.config dag.file, *.mmd) into an nf-metro metro-map diagram.
//
// Nextflow only finalizes the .mmd file AFTER every workflow.onComplete
// handler has returned (it's part of session shutdown, not concurrent with
// onComplete), so the file can never be read synchronously from in here —
// no amount of in-process polling helps. Instead, launch a fully detached
// background shell script (`&` + `disown`) that survives independently of
// this JVM: it polls for the file, lazily resolves nf-metro (PATH -> a
// dedicated conda env built from envs/pgirl_nf_metro.yml -> give up), and
// renders the diagram. All output/errors go to a log file next to the
// diagram; nothing here can fail or delay the pipeline itself.
//
def renderNfMetroDiagram(outdir) {
    def dag_mmd    = "${outdir}/pipeline_info/pipeline_dag_${params.trace_report_suffix}.mmd"
    def metro_html = "${outdir}/pipeline_info/pipeline_metro_map_${params.trace_report_suffix}.html"
    def metro_log  = "${outdir}/pipeline_info/.nf_metro.log"
    def env_prefix = "${workflow.projectDir}/.nf-metro-env"
    def env_yml    = "${workflow.projectDir}/envs/pgirl_nf_metro.yml"

    def assignments_tsv = "${outdir}/classification/species_assignments.tsv"
    def metro_mmd = "${dag_mmd}.metro.mmd"

    def script = """
        (
          for ((i=0; i<120; i++)); do
            [ -f '${dag_mmd}' ] && break
            sleep 1
          done
          if [ ! -f '${dag_mmd}' ]; then
            echo "nf-metro: DAG file not found after waiting: ${dag_mmd}"
            exit 0
          fi
          # Wait for species assignments (best-effort; not required for rendering)
          for ((i=0; i<60; i++)); do
            [ -f '${assignments_tsv}' ] && break
            sleep 1
          done
          # Rewrite the Mermaid DAG for hierarchical pathogen layout
          python3 ${workflow.projectDir}/bin/metro_layout_postprocess.py \\
            --input '${dag_mmd}' \\
            --assignments '${assignments_tsv}' \\
            --output '${metro_mmd}' || echo "nf-metro: post-processing failed, falling back to raw DAG"
          [ -f '${metro_mmd}' ] || cp '${dag_mmd}' '${metro_mmd}'

          NF_METRO_BIN=""
          if command -v nf-metro >/dev/null 2>&1; then
            NF_METRO_BIN="nf-metro"
          elif [ -x '${env_prefix}/bin/nf-metro' ]; then
            NF_METRO_BIN='${env_prefix}/bin/nf-metro'
          else
            CONDA_BIN=""
            for c in "\$CONDA_EXE" "\$CONDA_PREFIX/bin/conda" "\$HOME/anaconda3/bin/conda" "\$HOME/miniconda3/bin/conda" "\$HOME/miniforge3/bin/conda"; do
              if [ -n "\$c" ] && [ -x "\$c" ]; then CONDA_BIN="\$c"; break; fi
            done
            if [ -z "\$CONDA_BIN" ]; then
              echo "nf-metro: conda not found -- skipping (install with 'conda install bioconda::nf-metro' or 'pip install nf-metro')."
              exit 0
            fi
            echo "nf-metro: setting up dedicated conda env at ${env_prefix} (one-time, may take a minute)..."
            "\$CONDA_BIN" env create -f '${env_yml}' -p '${env_prefix}'
            if [ -x '${env_prefix}/bin/nf-metro' ]; then
              NF_METRO_BIN='${env_prefix}/bin/nf-metro'
            else
              echo "nf-metro: conda env creation failed -- skipping."
              exit 0
            fi
          fi
          if "\$NF_METRO_BIN" render '${metro_mmd}' --from-nextflow -o '${metro_html}' --title 'pgirl/genomic-intelligence' --section-x-gap 80 --section-y-gap 120 --track-gap 3 --y-spacing 80; then
            echo "nf-metro: pipeline metro-map diagram written to ${metro_html}"
          else
            echo "nf-metro: render failed"
          fi
        ) > '${metro_log}' 2>&1 &
        disown
    """.toString()

    try {
        ['bash', '-c', script].execute()
        log.info "nf-metro: scheduled background render (log: ${metro_log})"
    } catch (Exception e) {
        log.warn "nf-metro: failed to launch background render task (${e.message})."
    }
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
