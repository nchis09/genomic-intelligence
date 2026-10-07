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
 *   3. SEGMENT_MUTATIONS + PHENOTYPE_ANNOTATION — per-segment mutation
 *      TSVs (all built segments), UniProt/UniProtExtractR/rbioapi and
 *      hmmscan coverage for every protein the query isolates carry.
 *
 * Literature retrieval, epi data and the knowledge-warehouse bundle are
 * wired; skipped stages contribute NO_FILE placeholders so the shared
 * knowledge warehouse still loads influenza samples + the build trees.
 *
 * Input:  ch_species_data         - channel of [ meta, fasta, metadata ]
 *         ch_nextclade_aligned    - channel: collected [ meta(dataset), fasta ]
 *         ch_species_assignments  - path: species_assignments.tsv (broadcast)
 */

include { NEXTSTRAIN_FLU_INGEST } from '../../../modules/local/nextstrain_flu_ingest/main'
include { NEXTSTRAIN_FLU        } from '../../../modules/local/nextstrain_flu/main'
include { SEGMENT_MUTATIONS     } from '../../../modules/local/segment_mutations/main'
include { SEGMENT_SIGNATURES    } from '../../../modules/local/segment_signatures/main'
include { PHENOTYPE_ANNOTATION  } from '../phenotype_annotation/main'
include { LITERATURE_RETRIEVAL  } from '../literature_retrieval/main'
include { EPIDEMIOLOGICAL_DATA  } from '../epidemiological_data/main'

workflow SEASONAL_FLU_WORKFLOW {
    take:
    ch_species_data        // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_nextclade_aligned   // channel: collected [ val(meta), path(fasta) ]; each meta has .dataset
    ch_species_assignments // path: species_assignments.tsv (broadcast/value channel)
    ch_nextclade_tsvs      // path: all Nextclade TSVs (collected, broadcast)

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
    // Epidemiological data — WHO FluMart strain-specific weekly surveillance
    //
    EPIDEMIOLOGICAL_DATA(ch_species_data)
    ch_epi_raw         = EPIDEMIOLOGICAL_DATA.out.epi_raw
    ch_epi_search_summary = EPIDEMIOLOGICAL_DATA.out.search_summary

    //
    // Per-segment mutation tables + phenotype annotation.
    // The build emits one Auspice JSON per segment; SEGMENT_MUTATIONS
    // produces mutations_{seg}.tsv per built segment (nextclade on the
    // build's aligned.fasta, or Auspice-derived), and EXTRACT_QUERY_PROTEINS
    // merges them so UniProt/HMM coverage spans every segment — not just
    // the ha tree used for phylo-neighbor discovery.
    //
    ch_auspice_results = NEXTSTRAIN_FLU.out.auspice
        .join(NEXTSTRAIN_FLU.out.results_dir, by: 0)
        .map { meta, json, dir ->
            def jsons = json instanceof List ? json : [ json ]
            def pick = jsons.find { it && it.name ==~ /.*_ha\.json/ }
                ?: jsons.find { it && it.name ==~ /.*\.json/ }
            [ meta, pick ?: file('NO_FILE_auspice'), jsons, dir ]
        }

    ch_segmut_input = NEXTSTRAIN_FLU.out.results_dir
        .join(NEXTSTRAIN_FLU.out.auspice, by: 0, remainder: true)
        .map { meta, dir, json -> [ meta, dir, json ?: file('NO_FILE_auspice') ] }

    // SEGMENT_MUTATIONS is cheap local parsing and feeds both the phenotype
    // annotation chain AND the segment-signature tables, so it runs outside
    // the skip_phenotype_annotation guard.
    SEGMENT_MUTATIONS(ch_segmut_input)

    if (!params.skip_phenotype_annotation) {
        PHENOTYPE_ANNOTATION(
            ch_auspice_results,
            SEGMENT_MUTATIONS.out.mutations_tsvs,
            ch_species_assignments
        )
    }

    //
    // Influenza signature tables — constellation/reassortment, DI screen,
    // curated-marker hits and group_id diversity. Needs the group's mutation
    // TSVs + build dir + metadata; GenoFLU is avian-only so a NO_FILE
    // placeholder fills that input slot here.
    //
    if (!params.skip_segment_signatures) {
        def no_file_genoflu = file('NO_FILE_genoflu')
        ch_sig_input = SEGMENT_MUTATIONS.out.mutations_tsvs
            .join(NEXTSTRAIN_FLU.out.results_dir, by: 0)
            .join(ch_species_data.map { meta, fasta, metadata -> [ meta, metadata ] }, by: 0)
            .map { meta, tsvs, dir, metadata ->
                [ meta, tsvs, dir, metadata, no_file_genoflu ]
            }
        SEGMENT_SIGNATURES(
            ch_sig_input,
            ch_species_assignments,
            ch_nextclade_tsvs,
            Channel.fromPath("${projectDir}/database/flu_markers.yml").first()
        )
    }

    // Per-meta NO_FILE fallbacks so kw joins and emits stay populated when
    // phenotype annotation is skipped (results_dir is always emitted, so it
    // is a safe base channel unlike the optional auspice emit).
    ch_mutations_out = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.mutations

    ch_pheno_summary = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.query_summary

    ch_uniprotr_results = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.uniprotr_results

    ch_extractr_results = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.extractr_results

    ch_rbioapi_results = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.rbioapi_results

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

        ch_epi_for_kw = ch_epi_raw
            .map { meta, epi_dir ->
                [ meta, (epi_dir && !epi_dir.name.startsWith('NO_FILE')) ? epi_dir : no_file_epi ]
            }
        ch_epi_summary_for_kw = ch_epi_search_summary
            .map { meta, summary ->
                [ meta, (summary && !summary.name.startsWith('NO_FILE')) ? summary : no_file_summary ]
            }

        // Phenotype-annotation outputs for the kw bundle — per-meta
        // NO_FILE-safe channels keyed on results_dir (always emitted).
        ch_uniprotr_for_kw = params.skip_phenotype_annotation
            ? NEXTSTRAIN_FLU.out.results_dir.map { meta, _dir -> [ meta, no_file_uniprotr ] }
            : PHENOTYPE_ANNOTATION.out.uniprotr_results
                .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_uniprotr ] }
        ch_extractr_for_kw = params.skip_phenotype_annotation
            ? NEXTSTRAIN_FLU.out.results_dir.map { meta, _dir -> [ meta, no_file_extractr ] }
            : PHENOTYPE_ANNOTATION.out.extractr_results
                .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_extractr ] }
        ch_rbioapi_for_kw = params.skip_phenotype_annotation
            ? NEXTSTRAIN_FLU.out.results_dir.map { meta, _dir -> [ meta, no_file_rbioapi ] }
            : PHENOTYPE_ANNOTATION.out.rbioapi_results
                .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_rbioapi ] }
        ch_hmm_for_kw = params.skip_phenotype_annotation
            ? NEXTSTRAIN_FLU.out.results_dir.map { meta, _dir -> [ meta, no_file_hmm ] }
            : PHENOTYPE_ANNOTATION.out.hmm_results
                .map { meta, files -> [ meta, files ?: no_file_hmm ] }

        // EXTRACT_QUERY_PROTEINS outputs bundled as one per-meta list —
        // plus the per-segment mutation TSVs from SEGMENT_MUTATIONS.
        ch_query_data_for_kw = params.skip_phenotype_annotation
            ? NEXTSTRAIN_FLU.out.results_dir.map { meta, _dir -> [ meta, [ no_file_query ] ] }
            : PHENOTYPE_ANNOTATION.out.discovery
                .join(PHENOTYPE_ANNOTATION.out.accessions, by: 0)
                .join(PHENOTYPE_ANNOTATION.out.uniprot_tsv, by: 0)
                .join(PHENOTYPE_ANNOTATION.out.mutations, by: 0)
                .join(PHENOTYPE_ANNOTATION.out.query_proteins, by: 0)
                .join(PHENOTYPE_ANNOTATION.out.query_summary, by: 0)
                .join(SEGMENT_MUTATIONS.out.mutations_tsvs, by: 0, remainder: true)
                .map { meta, discovery, accessions, uniprot_tsv, mutations, proteins, summary, seg_tsvs ->
                    [ meta, [ discovery, accessions, uniprot_tsv, mutations, proteins, summary ] +
                        (seg_tsvs instanceof List ? seg_tsvs : [ seg_tsvs ]).findAll { it != null } ]
                }

        ch_kw_input = ch_metadata_for_kw
            .join(ch_bioinfo_for_kw, by: 0)
            .join(ch_auspice_for_kw, by: 0)
            .join(ch_epi_for_kw, by: 0)
            .join(ch_epi_summary_for_kw, by: 0)
            .join(ch_uniprotr_for_kw, by: 0, remainder: true)
            .join(ch_extractr_for_kw, by: 0, remainder: true)
            .join(ch_rbioapi_for_kw, by: 0, remainder: true)
            .join(ch_query_data_for_kw, by: 0, remainder: true)
            .join(ch_hmm_for_kw, by: 0, remainder: true)
            .combine(ch_species_assignments_kw)
            .map { meta, metadata, bioinfo_dir, auspice, tree, epi_dir, epi_summary,
                   uniprotr_dir, extractr_dir, rbioapi_dir, query_data, hmm, assignments ->
                [ meta, assignments, metadata, epi_dir, epi_summary,
                  bioinfo_dir, uniprotr_dir ?: no_file_uniprotr,
                  extractr_dir ?: no_file_extractr, rbioapi_dir ?: no_file_rbioapi,
                  auspice, tree, query_data ?: [ no_file_query ], hmm ?: no_file_hmm ]
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
    mutations          = ch_mutations_out                   // channel: [ meta, tsv ] (empty when skipped)
    query_summary      = ch_pheno_summary                   // channel: [ meta, json ]
    uniprotr_results   = ch_uniprotr_results                // channel: [ meta, dir ]
    extractr_results   = ch_extractr_results                // channel: [ meta, dir ]
    rbioapi_results    = ch_rbioapi_results                 // channel: [ meta, dir ]
    epi_raw            = ch_epi_raw
    epi_search_summary = ch_epi_search_summary
    lit_results        = ch_lit_results                       // channel: [ meta, [ literature result files ] ]
    lit_evidence       = ch_lit_evidence                      // channel: [ meta, evidence_extracted.tsv ]
}
