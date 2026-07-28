/*
 * Local module: PHYLO_ANNOTATE
 *
 * Extract tip/node metadata and mutation matrices from Nextstrain
 * Auspice JSON and muts JSON for phylogenetic tree plotting.
 */

process PHYLO_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container null

    input:
    tuple val(meta), path(auspice_json), path(results_dir)

    output:
    tuple val(meta), path("*_tip_metadata.tsv")    , emit: tip_metadata
    tuple val(meta), path("*_mutation_matrix.tsv")  , emit: mutation_matrix, optional: true
    tuple val(meta), path("*_mutation_legend.tsv")  , emit: mutation_legend, optional: true
    tuple val(meta), path("*_node_metadata.tsv")    , emit: node_metadata, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def query  = params.query_samples ?: (meta.query_samples ?: '')
    """
    #!/usr/bin/env python3
    import json, sys, glob
    from pathlib import Path

    auspice_file = "${auspice_json}"
    results_dir  = "${results_dir}"
    prefix       = "${prefix}"
    query_list   = [s.strip() for s in "${query}".split(",") if s.strip()]

    # Find muts.json and branch_lengths.json inside results_dir
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

    # --- Mutation counts from muts.json ---
    if Path(muts_file).exists():
        with open(muts_file) as f:
            muts = json.load(f)
        gene_mutation_counts = {}
        tip_muts_data = {}
        for node, data in muts["nodes"].items():
            if node.startswith("NODE_"):
                continue
            aa = data.get("aa_muts", {})
            tip_muts_data[node] = {}
            for gene, mut_list in aa.items():
                if gene not in gene_mutation_counts:
                    gene_mutation_counts[gene] = {}
                present = set(mut_list)
                tip_muts_data[node][gene] = present
                for m in present:
                    gene_mutation_counts[gene][m] = gene_mutation_counts[gene].get(m, 0) + 1
            for row in tips:
                if row["label"] == node:
                    row["nuc_mutation_count"] = len(data.get("muts", []))
                    row["aa_mutation_count"] = sum(len(v) for v in aa.values())

        # Top 5 mutations per gene
        selected_columns = []
        for gene in sorted(gene_mutation_counts.keys()):
            counts = gene_mutation_counts[gene]
            if not counts:
                continue
            top = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[:5]
            for mut, cnt in top:
                selected_columns.append(f"{gene}:{mut}")

        # Write mutation matrix
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

    # --- Node metadata from branch_lengths ---
    nodes = {}
    if Path(bl_file).exists():
        with open(bl_file) as f:
            bl = json.load(f)
        for node, data in bl["nodes"].items():
            nodes.setdefault(node, {"label": node})
            nodes[node]["branch_length"] = data.get("branch_length")

    if Path(muts_file).exists():
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
