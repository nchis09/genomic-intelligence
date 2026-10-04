/*
 * Subworkflow: AVIAN_INFLUENZA_WORKFLOW
 *
 * Pathogen-specific workflow for avian influenza species groups
 * (h5nx via the community iav-h5 Nextclade dataset, narrowed to h5n1 when
 * GenoFLU resolves a clean genotype):
 *   1. GENOFLU_ASSIGN          — vendored GenoFLU-multi per-isolate
 *      constellation genotyping; resolves h5nx -> h5n1 when every
 *      isolate assigns a clean genotype
 *   2. NEXTSTRAIN_AVIAN_INGEST — NCBI GenBank background for the resolved
 *      serotype(s) (all H5 while the species stays h5nx)
 *   3. NEXTSTRAIN_AVIAN        — segment-focused augur build
 *      (subtype/segment/all-time) with per-segment query sequences
 *      injected as additional_inputs
 *
 * Stage 1 wires only the bioinformatics path. Epi data, phenotype
 * annotation and literature retrieval are later stages; their emit
 * channels stay empty (safe to .mix()) while kw_input carries NO_FILE
 * placeholders so the shared knowledge warehouse can still load
 * avian samples + the screening/build trees.
 *
 * Input:  ch_species_data         - channel of [ meta, fasta, metadata ]
 *         ch_nextclade_aligned    - channel of [ meta(dataset), fasta ] (collected)
 *         ch_species_assignments  - path: species_assignments.tsv (broadcast)
 */

include { GENOFLU_ASSIGN        } from '../../../modules/local/genoflu_assign/main'
include { NEXTSTRAIN_AVIAN_INGEST } from '../../../modules/local/nextstrain_avian_ingest/main'
include { NEXTSTRAIN_AVIAN        } from '../../../modules/local/nextstrain_avian/main'
include { LITERATURE_RETRIEVAL    } from '../literature_retrieval/main'
include { EPIDEMIOLOGICAL_DATA    } from '../epidemiological_data/main'

// Resolve the build subtype from a genoflu_results.tsv: GenoFLU's
// genotype panel is entirely H5N1, so a clean call on EVERY isolate
// narrows h5nx -> h5n1. Any "Not assigned" or parse failure returns null
// (stays broad). Top-level function because strict DSL2 cannot invoke
// closures defined inside `main:` from nested operators.
def genoflu_species(tsv) {
    try {
        if (!tsv || !tsv.exists() || tsv.size() == 0) { return null }
        def lines = tsv.readLines()
        if (lines.size() < 2) { return null }
        def hdr = lines[0].split('\t') as List
        def gi = hdr.indexOf('Genotype')
        if (gi < 0) { return null }
        def genos = lines.drop(1)
            .collect { it.split('\t')[gi]?.trim() }
            .findAll { it }
        if (genos && genos.every { !it.startsWith('Not assigned') }) { return 'h5n1' }
    } catch (Throwable ignored) {}
    return null
}

workflow AVIAN_INFLUENZA_WORKFLOW {
    take:
    ch_species_data        // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_nextclade_aligned   // channel: [ val(meta), path(fasta) ] (value channel; each meta has .dataset)
    ch_species_assignments // path: species_assignments.tsv (broadcast/value channel)

    main:
    //
    // GenoFLU genotyping — runs ONLY on groups the pathogen router already
    // resolved as avian_influenza (this subworkflow never sees other
    // pathogen groups). A clean genotype call on every isolate narrows
    // h5nx -> h5n1; a reassortant or unresolved group stays broad so the
    // wide H5 ingest/build still applies.
    //
    GENOFLU_ASSIGN(ch_species_data.map { meta, fasta, metadata -> [ meta, fasta, metadata ] })

    ch_species_data = ch_species_data
        .join(GENOFLU_ASSIGN.out.results, by: 0)
        .map { meta, fasta, metadata, results ->
            def resolved = genoflu_species(results) ?: meta.species
            [ meta + [species: resolved, id: "${meta.pathogen}_${resolved}"], fasta, metadata ]
        }

    //
    // Background data: either freshly ingested NCBI files or the staged
    // copy left in data/avian_flu/ingest/ncbi/results by a previous run.
    //
    ch_ingest_bg = params.skip_nextstrain_ingest
        ? ch_species_data.map { meta, fasta, metadata ->
            def bg = new File("${projectDir}/data/avian_flu/ingest/ncbi/results")
            // path() inputs need nio Paths — listFiles() returns io.Files
            def seqs = (bg.listFiles()?.findAll { it.name ==~ /^sequences_.*\.fasta/ } ?: [])
                .collect { file(it.toString()) }
            [ meta, file("${bg}/metadata.tsv"), seqs ]
        }
        : NEXTSTRAIN_AVIAN_INGEST(ch_species_data.map { it[0] }).background

    //
    // Pick the aligned query FASTA from the Nextclade run for the group's
    // winning dataset (meta.dataset on the group == meta.dataset on the
    // nextclade run). Used to extract per-segment query sequences in the
    // build so whole-genome and per-segment inputs both work.
    //
    // NOTE: the collected aligned channel can arrive flat ([meta,path,meta,…])
    // because collect() flattens tuple emissions; re-pair with collate(2) —
    // or use as-is if it arrives already nested.
    ch_aligned_for_group = ch_species_data
        .combine(ch_nextclade_aligned)
        .map { it ->
            def meta  = it[0]
            def tail  = it.drop(3)
            def alist = (tail && tail[0] instanceof List) ? tail : tail.collate(2)
            def hit = alist.find { pair -> pair[0].dataset == meta.dataset }
            [ meta, hit ? hit[1] : file('NO_FILE_aligned') ]
        }

    ch_avian_input = ch_species_data
        .join(ch_aligned_for_group)
        .join(ch_ingest_bg)
        .combine(ch_species_assignments)
        .map { meta, fasta, metadata, aligned, bg_metadata, bg_fastas, assignments ->
            [ meta, fasta, metadata, aligned, bg_metadata, bg_fastas, assignments ]
        }

    NEXTSTRAIN_AVIAN(ch_avian_input)

    //
    // Literature retrieval — Europe PMC per-domain search, PubMed metadata,
    // dedup, ASReview screening, OA PDF download, PDF→text, evidence
    // extraction + per-domain summaries for the resolved avian species.
    //
    LITERATURE_RETRIEVAL(ch_species_data)
    ch_lit_results   = LITERATURE_RETRIEVAL.out.lit_results
    ch_lit_evidence  = LITERATURE_RETRIEVAL.out.lit_evidence
    ch_lit_summaries = LITERATURE_RETRIEVAL.out.lit_summaries

    //
    // Epidemiological data — WHO FluMart H5 surveillance (human sentinel detections)
    //
    EPIDEMIOLOGICAL_DATA(ch_species_data)
    ch_epi_raw         = EPIDEMIOLOGICAL_DATA.out.epi_raw
    ch_epi_search_summary = EPIDEMIOLOGICAL_DATA.out.search_summary

    //
    // Knowledge-warehouse bundle: same tuple shape as EBOLA_WORKFLOW's,
    // with NO_FILE placeholders for the stages not wired for avian yet
    // (epi data, phenotype annotation, HMM, query-protein discovery).
    //
    if (!params.skip_knowledge_warehouse) {
        ch_species_assignments_kw = ch_species_assignments.map { it }

        def no_file_meta     = file('NO_FILE_metadata')
        def no_file_epi      = file('NO_FILE_epi')
        def no_file_summary  = file('NO_FILE_epi_summary')
        def no_file_uniprotr = file('NO_FILE_uniprotr')
        def no_file_extractr = file('NO_FILE_extractr')
        def no_file_rbioapi  = file('NO_FILE_rbioapi')
        def no_file_tree     = file('NO_FILE_tree')
        def no_file_hmm      = file('NO_FILE_hmm')
        def no_file_query    = file('NO_FILE_query_data')

        ch_metadata_for_kw = ch_species_data
            .map { meta, fasta, metadata ->
                [ meta, (metadata && !metadata.name.startsWith('NO_FILE')) ? metadata : no_file_meta ]
            }

        ch_bioinfo_for_kw = NEXTSTRAIN_AVIAN.out.results_dir
            .map { meta, dir -> [ meta, dir ?: file('NO_FILE_bioinfo') ] }

        // The build emits a glob of Auspice JSONs (one per segment) plus
        // tip-frequencies sidecars — the kw bundle slot takes a single
        // file, so prefer the ha tree and drop the sidecars.
        ch_auspice_for_kw = NEXTSTRAIN_AVIAN.out.auspice
            .map { meta, json ->
                def jsons = json instanceof List ? json : [ json ]
                def pick = jsons.find { it && it.name ==~ /.*_ha_.*\.json/ && !it.name.contains('tip-frequencies') }
                    ?: jsons.find { it && it.name ==~ /.*\.json/ && !it.name.contains('tip-frequencies') }
                [ meta, pick ?: file('NO_FILE_auspice'), no_file_tree ]
            }

        ch_epi_for_kw = ch_epi_raw
            .map { meta, epi_dir ->
                [ meta, (epi_dir && !epi_dir.name.startsWith('NO_FILE')) ? epi_dir : no_file_epi ]
            }
        ch_epi_summary_for_kw = ch_epi_search_summary
            .map { meta, summary ->
                [ meta, (summary && !summary.name.startsWith('NO_FILE')) ? summary : no_file_summary ]
            }

        ch_kw_input = ch_metadata_for_kw
            .join(ch_bioinfo_for_kw, by: 0)
            .join(ch_auspice_for_kw, by: 0)
            .join(ch_epi_for_kw, by: 0)
            .join(ch_epi_summary_for_kw, by: 0)
            .combine(ch_species_assignments_kw)
            .map { meta, metadata, bioinfo_dir, auspice, tree, epi_dir, epi_summary, assignments ->
                [ meta, assignments, metadata, epi_dir, epi_summary,
                  bioinfo_dir, no_file_uniprotr, no_file_extractr, no_file_rbioapi,
                  auspice, tree, [ no_file_query ], no_file_hmm ]
            }

        // Dataflow wait: join the bundle on this species' literature
        // summaries — the bundle can't emit until EVIDENCE_SUMMARY publishes
        // (the lit chain's last stage). remainder:true releases it with
        // ready=false when the species produced no lit outputs or the chain
        // was skipped, so the warehouse still builds.
        ch_kw_input = ch_kw_input
            .map { it -> [ it[0].species ] + it }
            .join(
                ch_lit_summaries.map { meta, sums -> [ meta.species, sums ] },
                by: 0, remainder: true
            )
            .map { it ->
                [ it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[8], it[9], it[10], it[11], it[12], it[13], it[14] != null ]
            }
    } else {
        ch_kw_input = channel.empty()
    }

    emit:
    kw_input           = ch_kw_input                        // channel: kw bundle (see EBOLA_WORKFLOW for shape)
    auspice            = NEXTSTRAIN_AVIAN.out.auspice       // channel: [ meta, json ]
    results            = NEXTSTRAIN_AVIAN.out.results_dir   // channel: [ meta, dir ]
    mutations          = channel.empty()                    // later stage: mutation profile
    query_summary      = channel.empty()
    uniprotr_results   = channel.empty()
    extractr_results   = channel.empty()
    rbioapi_results    = channel.empty()
    epi_raw            = ch_epi_raw
    epi_search_summary = ch_epi_search_summary
    lit_results        = ch_lit_results                       // channel: [ meta, [ literature result files ] ]
    lit_evidence       = ch_lit_evidence                      // channel: [ meta, evidence_extracted.tsv ]
}
