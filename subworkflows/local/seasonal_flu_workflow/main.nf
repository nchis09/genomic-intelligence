/*
 * Subworkflow: SEASONAL_FLU_WORKFLOW
 *
 * Pathogen-specific workflow for seasonal influenza species groups
 * (h1n1pdm, h3n2, vic — vic also carries yam/generic-B query groups, which
 * have no dedicated ingest source or upstream lineage):
 *   1. NEXTSTRAIN_FLU_INGEST — GenSpectrum/LAPIS public background, all 8
 *      genome segments
 *   2. NEXTSTRAIN_FLU — vendored nextstrain/seasonal-flu build
 *      ({lineage}_pgirl, ha+na, 6y) with the group's query segments injected
 *      as additional_inputs
 *
 * Stage 1 wires only the bioinformatics path. Epi data, phenotype
 * annotation and literature retrieval are later stages; their emit
 * channels stay empty (safe to .mix()) while kw_input carries NO_FILE
 * placeholders so the shared knowledge warehouse can still load
 * influenza samples + the screening/build trees.
 *
 * Input:  ch_species_data         - channel of [ meta, fasta, metadata ]
 *         ch_nextclade_aligned    - channel: collected [ meta(dataset), fasta ]
 *         ch_species_assignments  - path: species_assignments.tsv (broadcast)
 */

include { NEXTSTRAIN_FLU_INGEST } from '../../../modules/local/nextstrain_flu_ingest/main'
include { NEXTSTRAIN_FLU        } from '../../../modules/local/nextstrain_flu/main'
include { LITERATURE_RETRIEVAL  } from '../literature_retrieval/main'

workflow SEASONAL_FLU_WORKFLOW {
    take:
    ch_species_data        // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_nextclade_aligned   // channel: collected [ val(meta), path(fasta) ]; each meta has .dataset
    ch_species_assignments // path: species_assignments.tsv (broadcast/value channel)

    main:
    //
    // Background data: freshly ingested GenSpectrum files, or the staged
    // copy in data/seasonal_flu/ingest-open/results/{lineage}.
    //
    ch_ingest_bg = params.skip_nextstrain_ingest
        ? ch_species_data.map { meta, fasta, metadata ->
            def bg = new File("${projectDir}/data/seasonal_flu/ingest-open/results/${meta.species}")
            // path() inputs need nio Paths — listFiles() returns io.Files
            def seqs = (bg.listFiles()?.findAll { it.name ==~ /.*\.fasta/ } ?: [])
                .collect { file(it.toString()) }
            [ meta, file("${bg}/metadata.tsv"), seqs ]
        }
        : NEXTSTRAIN_FLU_INGEST(ch_species_data.map { it[0] }).background

    //
    // Collect every aligned output for this lineage's segment datasets —
    // the build can now emit a tree per segment, not just ha/na. The
    // group's winning dataset may be yam/b (mapped onto vic), so match on
    // the lineage prefix rather than meta.dataset directly. Two parallel
    // lists are emitted: the aligned FASTAs and the segment each was typed
    // against (token 3 of nextstrain/flu/{lineage}/{segment}/...).
    //
    // NOTE: the collected aligned channel can arrive flat ([meta,path,meta,…])
    // because collect() flattens tuple emissions; re-pair with collate(2) —
    // or use as-is if it arrives already nested.
    ch_aligned_segments = ch_species_data
        .combine(ch_nextclade_aligned)
        .map { it ->
            def meta  = it[0]
            def tail  = it.drop(3)
            def alist = (tail && tail[0] instanceof List) ? tail : tail.collate(2)
            def hits = alist.findAll { pair ->
                pair[0].dataset?.startsWith("nextstrain/flu/${meta.species}/")
            }
            def segs = hits.collect { pair -> pair[0].dataset.tokenize('/')[3] }
            [ meta, hits.collect { it[1] }, segs ]
        }

    ch_flu_input = ch_species_data
        .join(ch_aligned_segments)
        .join(ch_ingest_bg)
        .combine(ch_species_assignments)
        .map { meta, fasta, metadata, aligned_fastas, aligned_segs, bg_metadata, bg_fastas, assignments ->
            [ meta, fasta, metadata, aligned_fastas, aligned_segs, bg_metadata, bg_fastas, assignments ]
        }

    NEXTSTRAIN_FLU(ch_flu_input)

    //
    // Literature retrieval — Europe PMC per-domain search, PubMed metadata,
    // dedup, ASReview screening, OA PDF download, PDF→text, evidence
    // extraction + per-domain summaries for the resolved seasonal species.
    //
    LITERATURE_RETRIEVAL(ch_species_data)
    ch_lit_results   = LITERATURE_RETRIEVAL.out.lit_results
    ch_lit_evidence  = LITERATURE_RETRIEVAL.out.lit_evidence
    ch_lit_summaries = LITERATURE_RETRIEVAL.out.lit_summaries

    //
    // Knowledge-warehouse bundle: same tuple shape as the other pathogen
    // workflows, with NO_FILE placeholders for the stages not wired for
    // influenza yet (epi data, phenotype annotation, HMM, query-protein
    // discovery).
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

        ch_bioinfo_for_kw = NEXTSTRAIN_FLU.out.results_dir
            .map { meta, dir -> [ meta, dir ?: file('NO_FILE_bioinfo') ] }

        // The build emits one Auspice JSON per segment (ha + na) plus
        // tip-frequencies sidecars — the kw bundle slot takes a single file,
        // so prefer the ha tree (primary surveillance layer).
        ch_auspice_for_kw = NEXTSTRAIN_FLU.out.auspice
            .map { meta, json ->
                def jsons = json instanceof List ? json : [ json ]
                def pick = jsons.find { it && it.name ==~ /.*_ha\.json/ }
                    ?: jsons.find { it && it.name ==~ /.*\.json/ }
                [ meta, pick ?: file('NO_FILE_auspice'), no_file_tree ]
            }

        ch_kw_input = ch_metadata_for_kw
            .join(ch_bioinfo_for_kw, by: 0)
            .join(ch_auspice_for_kw, by: 0)
            .combine(ch_species_assignments_kw)
            .map { meta, metadata, bioinfo_dir, auspice, tree, assignments ->
                [ meta, assignments, metadata, no_file_epi, no_file_summary,
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
    kw_input           = ch_kw_input                        // channel: kw bundle (see AVIAN_INFLUENZA_WORKFLOW)
    auspice            = NEXTSTRAIN_FLU.out.auspice         // channel: [ meta, json ]
    results            = NEXTSTRAIN_FLU.out.results_dir     // channel: [ meta, dir ]
    mutations          = channel.empty()                    // later stage: mutation profile
    query_summary      = channel.empty()
    uniprotr_results   = channel.empty()
    extractr_results   = channel.empty()
    rbioapi_results    = channel.empty()
    epi_raw            = channel.empty()
    epi_search_summary = channel.empty()
    lit_results        = ch_lit_results                       // channel: [ meta, [ literature result files ] ]
    lit_evidence       = ch_lit_evidence                      // channel: [ meta, evidence_extracted.tsv ]
}
