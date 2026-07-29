/*
 * Local module: PHYLO_ANNOTATE
 *
 * Extract tip/node metadata and mutation matrices from Nextstrain Auspice
 * JSON for phylogenetic tree plotting.
 *
 * AA mutations for QUERY tips are now taken directly from each sample's own
 * Nextclade run (aaSubstitutions/aaDeletions/aaInsertions, relative to the
 * reference genome) instead of muts.json's branch-specific ancestral
 * reconstruction (which only reflects mutations on that tip's terminal
 * branch relative to its parent, not the full path from the reference).
 * BACKGROUND tips (curated Nextstrain ingest data) have no Nextclade run of
 * their own, so they still fall back to muts.json's per-node aa_muts.
 *
 * The mutation-matrix panel now shows the top-5 GLOBALLY most frequent AA
 * mutations (not top-5 per gene), and mutations are cross-referenced
 * against rbioapi_results (mutagenesis.tsv/variation.tsv query_mutation_match)
 * to flag ones with documented phenotypic evidence. A separate
 * *_protein_burden.tsv summarises total AA mutation counts per viral
 * protein across the whole dataset (for the 7-bar summary panel).
 */

process PHYLO_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container null

    input:
    tuple val(meta), path(auspice_json), path(results_dir), path(rbioapi_dir)
    path nextclade_jsons
    path species_assignments

    output:
    tuple val(meta), path("*_tip_metadata.tsv")     , emit: tip_metadata
    tuple val(meta), path("*_mutation_matrix.tsv")  , emit: mutation_matrix, optional: true
    tuple val(meta), path("*_mutation_legend.tsv")  , emit: mutation_legend, optional: true
    tuple val(meta), path("*_protein_burden.tsv")   , emit: protein_burden, optional: true
    tuple val(meta), path("*_node_metadata.tsv")    , emit: node_metadata, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def query  = params.query_samples ?: (meta.query_samples ?: '')
    // The 7 canonical Ebola proteins for the burden summary panel.
    // GP's alternate reading frame (Nextclade cdsName "GP_003") is folded
    // into the "GP" bucket rather than shown as its own bar.
    def proteins = task.ext.proteins ?: 'NP,VP35,VP40,GP,VP30,VP24,L'
    """
    #!/usr/bin/env python3
    import json, sys, glob, csv
    from pathlib import Path
    from collections import defaultdict

    auspice_file = "${auspice_json}"
    results_dir  = "${results_dir}"
    rbioapi_dir  = "${rbioapi_dir}"
    prefix       = "${prefix}"
    query_list   = [s.strip() for s in "${query}".split(",") if s.strip()]
    protein_list = [p.strip() for p in "${proteins}".split(",") if p.strip()]
    protein_alias = {"GP_003": "GP"}  # fold alternate ORF into its parent protein

    # Find muts.json and branch_lengths.json inside results_dir (background-tip fallback)
    muts_files = glob.glob(f"{results_dir}/**/muts.json", recursive=True)
    bl_files   = glob.glob(f"{results_dir}/**/branch_lengths.json", recursive=True)
    muts_file  = muts_files[0] if muts_files else "NO_FILE"
    bl_file    = bl_files[0] if bl_files else "NO_FILE"

    # --- Helper ---
    def get_value(node_attr, key):
        val = node_attr.get(key)
        if isinstance(val, dict):
            return val.get("value")
        return val

    def collect_tip_attrs(tree):
        tips = []
        children = tree.get("children")
        if not children:
            row = {"label": tree["name"]}
            node_attrs = tree.get("node_attrs", {})
            for k in node_attrs:
                row[k] = get_value(node_attrs, k)
            tips.append(row)
        else:
            for child in children:
                tips.extend(collect_tip_attrs(child))
        return tips

    # --- Tip metadata from Auspice JSON ---
    with open(auspice_file) as f:
        auspice = json.load(f)
    tips = collect_tip_attrs(auspice["tree"])
    for row in tips:
        row["is_query"] = row["label"] in query_list
    tip_labels = {row["label"] for row in tips}

    # -------------------------------------------------------------------
    # Query-tip AA mutations: direct from each sample's own Nextclade run
    # -------------------------------------------------------------------
    # species_assignments.tsv columns: sample, pathogen, species, qc_score, best_dataset_file
    sample_to_tsv = {}
    assignments_file = "${species_assignments}"
    if Path(assignments_file).exists() and Path(assignments_file).name != "NO_FILE":
        with open(assignments_file) as f:
            reader = csv.DictReader(f, delimiter="\\t")
            for row in reader:
                sample = row.get("sample", "").strip()
                best_tsv = row.get("best_dataset_file", "").strip()
                if sample and best_tsv:
                    sample_to_tsv[sample] = best_tsv

    # Map each Nextclade TSV basename to its matching JSON path (same
    # basename, .json extension) among the collected nextclade_jsons.
    nextclade_json_paths = "${nextclade_jsons}".split()
    tsv_stem_to_json = {}
    for jp in nextclade_json_paths:
        stem = Path(jp).stem
        if stem and Path(jp).name != "NO_FILE":
            tsv_stem_to_json[stem] = jp

    def load_nextclade_result(sample):
        best_tsv = sample_to_tsv.get(sample)
        if not best_tsv:
            return None
        json_path = tsv_stem_to_json.get(Path(best_tsv).stem)
        if not json_path or not Path(json_path).exists():
            return None
        with open(json_path) as f:
            data = json.load(f)
        for r in data.get("results", []):
            if r.get("seqName") == sample:
                return r
        return None

    # tip_muts_data: { tip_label: { protein: set(mutation_labels like "Y151F") } }
    # mutation_source: { "PROTEIN:MUT": "nextclade" | "ancestral" }
    tip_muts_data = {}
    mutation_source = {}

    for sample in tip_labels & set(sample_to_tsv.keys()):
        result = load_nextclade_result(sample)
        if not result:
            continue
        by_protein = defaultdict(set)
        total_aa = 0
        for sub in result.get("aaSubstitutions", []):
            gene = sub.get("cdsName", "")
            mut = f"{sub.get('refAa', '')}{int(sub.get('pos', 0)) + 1}{sub.get('qryAa', '')}"
            by_protein[gene].add(mut)
            mutation_source[f"{gene}:{mut}"] = "nextclade"
            total_aa += 1
        tip_muts_data[sample] = dict(by_protein)
        for row in tips:
            if row["label"] == sample:
                row["nuc_mutation_count"] = result.get("totalSubstitutions", "")
                row["aa_mutation_count"] = total_aa

    # -------------------------------------------------------------------
    # Background-tip AA mutations: fall back to muts.json ancestral aa_muts
    # (no Nextclade run exists for curated ingest/background sequences)
    # -------------------------------------------------------------------
    muts = None
    if Path(muts_file).exists():
        with open(muts_file) as f:
            muts = json.load(f)
        for node, data in muts["nodes"].items():
            if node.startswith("NODE_") or node in tip_muts_data:
                continue
            aa = data.get("aa_muts", {})
            by_protein = {}
            for gene, mut_list in aa.items():
                present = set(mut_list)
                if present:
                    by_protein[gene] = present
                    for m in present:
                        mutation_source.setdefault(f"{gene}:{m}", "ancestral")
            if by_protein:
                tip_muts_data[node] = by_protein
            for row in tips:
                if row["label"] == node:
                    row["nuc_mutation_count"] = len(data.get("muts", []))
                    row["aa_mutation_count"] = sum(len(v) for v in aa.values())

    # -------------------------------------------------------------------
    # Global mutation frequency (across ALL tips with data) -> top 5 overall
    # -------------------------------------------------------------------
    mutation_counts = defaultdict(int)  # "PROTEIN:MUT" -> count across tips
    for lab, by_protein in tip_muts_data.items():
        for gene, muts_set in by_protein.items():
            for m in muts_set:
                mutation_counts[f"{gene}:{m}"] += 1

    top5 = sorted(mutation_counts.items(), key=lambda x: (-x[1], x[0]))[:5]
    selected_columns = [col for col, _ in top5]

    # -------------------------------------------------------------------
    # Protein-level mutation burden: total AA mutation occurrences per
    # protein, summed across the whole dataset (raw counts, 7 bars)
    # -------------------------------------------------------------------
    protein_burden = defaultdict(int)
    for col, cnt in mutation_counts.items():
        gene = col.split(":", 1)[0]
        gene = protein_alias.get(gene, gene)
        if gene in protein_list:
            protein_burden[gene] += cnt

    with open(f"{prefix}_protein_burden.tsv", "w") as f:
        f.write("protein\\ttotal_aa_mutations\\n")
        for protein in protein_list:
            f.write(f"{protein}\\t{protein_burden.get(protein, 0)}\\n")

    # -------------------------------------------------------------------
    # Phenotype relevance: cross-reference query_mutation_match from
    # rbioapi_results (mutagenesis.tsv / variation.tsv), added by
    # annotate_rbioapi.R. A mutation_label there is formatted the same way
    # as our columns here (e.g. "GP:Y151F").
    # -------------------------------------------------------------------
    relevant_mutations = set()
    if Path(rbioapi_dir).exists() and Path(rbioapi_dir).name != "NO_FILE":
        for tsv_path in list(Path(rbioapi_dir).glob("*_mutagenesis.tsv")) + list(Path(rbioapi_dir).glob("*_variation.tsv")):
            try:
                with open(tsv_path) as f:
                    reader = csv.DictReader(f, delimiter="\\t")
                    for row in reader:
                        matches = row.get("query_mutation_match", "")
                        for m in matches.split(";"):
                            m = m.strip()
                            if m:
                                relevant_mutations.add(m)
            except Exception as e:
                print(f"  WARNING: could not parse {tsv_path}: {e}", file=sys.stderr)

    # --- Write mutation matrix + legend (top 5 overall) ---
    if selected_columns:
        with open(f"{prefix}_mutation_matrix.tsv", "w") as f:
            f.write("label\\t" + "\\t".join(selected_columns) + "\\n")
            for row in tips:
                lab = row["label"]
                present = tip_muts_data.get(lab, {})
                vals = []
                for col in selected_columns:
                    gene, mut = col.split(":", 1)
                    vals.append("1" if mut in present.get(gene, set()) else "0")
                f.write(lab + "\\t" + "\\t".join(vals) + "\\n")

        with open(f"{prefix}_mutation_legend.tsv", "w") as f:
            f.write("mutation\\tcount\\tsource\\thas_phenotype_evidence\\n")
            for col, cnt in top5:
                has_evidence = col in relevant_mutations
                source = mutation_source.get(col, "")
                f.write(f"{col}\\t{cnt}\\t{source}\\t{has_evidence}\\n")

    # --- Write tip metadata ---
    if tips:
        all_keys = set()
        for row in tips:
            all_keys.update(row.keys())
        header = ["label", "is_query"] + sorted(k for k in all_keys if k not in {"label", "is_query"})
        with open(f"{prefix}_tip_metadata.tsv", "w") as f:
            f.write("\\t".join(header) + "\\n")
            for row in tips:
                f.write("\\t".join(str(row.get(h, "")) for h in header) + "\\n")

    # --- Node metadata from branch_lengths (internal nodes; unaffected by
    # the Nextclade change since internal nodes have no direct sequence) ---
    nodes = {}
    if Path(bl_file).exists():
        with open(bl_file) as f:
            bl = json.load(f)
        for node, data in bl["nodes"].items():
            nodes.setdefault(node, {"label": node})
            nodes[node]["branch_length"] = data.get("branch_length")

    if muts is not None:
        for node, data in muts["nodes"].items():
            nodes.setdefault(node, {"label": node})
            nodes[node]["nuc_mutation_count"] = len(data.get("muts", []))
            aa = data.get("aa_muts", {})
            nodes[node]["aa_mutation_count"] = sum(len(v) for v in aa.values())

    if nodes:
        all_keys = set()
        for row in nodes.values():
            all_keys.update(row.keys())
        header = ["label"] + sorted(k for k in all_keys if k != "label")
        with open(f"{prefix}_node_metadata.tsv", "w") as f:
            f.write("\\t".join(header) + "\\n")
            for node in sorted(nodes.keys()):
                row = nodes[node]
                f.write("\\t".join(str(row.get(h, "")) for h in header) + "\\n")

    print(f"Annotation complete for {prefix}")
    """
}
