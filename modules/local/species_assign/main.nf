/*
 * Local module: SPECIES_ASSIGN
 *
 * Parse Nextclade TSV outputs from multiple dataset runs,
 * determine the best-matching dataset per sample (lowest qc.overallScore),
 * and group samples by pathogen/species.
 *
 * Outputs one FASTA + metadata file per detected species group.
 */

process SPECIES_ASSIGN {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python>=3.12"
    container null

    input:
    tuple val(meta), path(fasta), path(metadata)
    path nextclade_tsvs  // all Nextclade TSV outputs (collected)

    output:
    path "species_assignments.tsv"                        , emit: assignments
    path "*_species_*.fasta"                              , emit: species_fasta
    path "*_species_*.metadata.tsv"                       , emit: species_metadata
    path "species_groups.json"                            , emit: species_groups_json

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    import json, csv, os, sys
    from pathlib import Path
    from collections import defaultdict

    # --- Parse all Nextclade TSV files ---
    # Each TSV is from a different dataset run. Filename encodes the dataset info.
    # We need to figure out which dataset each TSV came from.
    # Nextclade TSVs have columns: seqName, clade, qc.overallScore, qc.overallStatus, ...

    tsv_files = sorted(Path(".").glob("*.tsv"))
    # Exclude the metadata input file from TSV glob
    metadata_file = "${metadata}"
    tsv_files = [f for f in tsv_files if f.name != Path(metadata_file).name
                 and f.name != "species_assignments.tsv"
                 and not f.name.endswith("_species_assignments.tsv")]

    # Parse each TSV and collect scores per sample per dataset-file
    # Key insight: each TSV file comes from running against one specific dataset
    dataset_results = {}  # { tsv_filename: { sample_name: qc_score } }

    for tsv_file in tsv_files:
        try:
            with open(tsv_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\\t")
                scores = {}
                for row in reader:
                    seq_name = row.get("seqName", "").strip()
                    try:
                        score = float(row.get("qc.overallScore", "999999"))
                    except (ValueError, TypeError):
                        score = 999999.0
                    if seq_name:
                        scores[seq_name] = score
                if scores:
                    dataset_results[tsv_file.name] = scores
        except Exception as e:
            print(f"Warning: could not parse {tsv_file}: {e}", file=sys.stderr)

    if not dataset_results:
        print("ERROR: No valid Nextclade TSV results found", file=sys.stderr)
        sys.exit(1)

    # --- Determine dataset name from TSV filename ---
    # Nextclade TSV files are named like: <prefix>.tsv
    # The prefix comes from the NEXTCLADE_RUN module ext.prefix or meta.id
    # We need to map filenames to dataset names
    # Convention: TSV files are named <sample_id>_<dataset_suffix>.tsv
    # e.g., "input_FASTA_bdbv.tsv" or "input_FASTA_sudan.tsv"
    # The dataset suffix is the last part of the nextclade dataset name

    # Extract dataset identifiers from filenames
    # Files are named: {meta.id}.tsv but we run once per dataset
    # Since combine() creates one run per dataset, each TSV has same meta.id
    # We differentiate by the actual content — use the clade/scores to determine

    # Better approach: read the nextclade dataset name from the TSV if available
    # Or infer from directory structure. Since we collect all TSVs, let's use
    # a simpler heuristic: the dataset with the BEST average score wins.

    # Collect all sample names across all datasets
    all_samples = set()
    for scores in dataset_results.values():
        all_samples.update(scores.keys())

    # For each sample, find the dataset (TSV file) with the lowest qc.overallScore
    sample_best = {}  # { sample_name: (best_score, best_tsv_filename) }
    for sample in all_samples:
        best_score = float("inf")
        best_tsv = None
        for tsv_name, scores in dataset_results.items():
            if sample in scores and scores[sample] < best_score:
                best_score = scores[sample]
                best_tsv = tsv_name
        sample_best[sample] = (best_score, best_tsv)

    # --- Map TSV filenames to dataset/species (generic, dataset-agnostic) ---
    # NEXTCLADE_DATASETGET downloads each dataset into a directory named after
    # its LAST path segment (Nextflow `path.name`), e.g.:
    #   nextstrain/orthoebolavirus/bdbv        -> dir "bdbv"
    #   nextstrain/sars-cov-2/wuhan-hu-1/orfs   -> dir "orfs"
    #   nextstrain/flu/h1n1pdm/ha/MW626062      -> dir "MW626062"
    # CLASSIFICATION then names each Nextclade run "<sample>_<dir>", so the
    # TSV filename always contains that exact leaf directory name — this is
    # what we must match on. Pathogen family is inferred from the path
    # segment right after the source prefix (nextstrain/community/...);
    # species is a human-readable label built from everything after that.
    # Pathogen-specific renaming needed by a downstream workflow (e.g.
    # mapping Ebola's nextclade suffix to its Nextstrain ingest dir name)
    # is handled inside that pathogen's own workflow, not here.
    nextclade_datasets_str = "${params.nextclade_datasets}".strip("[]")
    dataset_list = [d.strip().strip("'").strip('"') for d in nextclade_datasets_str.split(",")]

    # Build a map: leaf_dir_name -> (pathogen, species)
    leaf_to_info = {}
    for ds in dataset_list:
        parts = [p for p in ds.strip("/").split("/") if p]
        if len(parts) < 2:
            continue
        # First path segment is the dataset source (nextstrain/community/...);
        # the pathogen family is the next segment.
        if parts[0] in ("nextstrain", "community") and len(parts) > 2:
            pathogen = parts[1]
            remainder = parts[2:]
        else:
            pathogen = parts[0]
            remainder = parts[1:]
        leaf = remainder[-1] if remainder else pathogen
        species = "-".join(remainder) if remainder else pathogen
        leaf_to_info[leaf] = (pathogen, species)

    # Map each TSV file to a species based on filename containing the leaf dir name
    tsv_to_species = {}
    for tsv_name in dataset_results.keys():
        matched = False
        # Try longest leaf names first to avoid partial/ambiguous matches
        for leaf in sorted(leaf_to_info.keys(), key=len, reverse=True):
            pathogen, species = leaf_to_info[leaf]
            if leaf.lower() in tsv_name.lower():
                tsv_to_species[tsv_name] = (pathogen, species)
                matched = True
                break
        if not matched:
            # Fallback: use first dataset's pathogen, species = "unknown"
            if leaf_to_info:
                first_pathogen = list(leaf_to_info.values())[0][0]
                tsv_to_species[tsv_name] = (first_pathogen, "unknown")
            else:
                tsv_to_species[tsv_name] = ("unknown", "unknown")

    # --- Assign each sample to its best species ---
    sample_assignments = {}  # { sample: (pathogen, species, score) }
    for sample, (score, best_tsv) in sample_best.items():
        if best_tsv and best_tsv in tsv_to_species:
            pathogen, _ = tsv_to_species[best_tsv]
            # Use the short leaf name so downstream processes and DB keys agree
            leaf = [l for l in leaf_to_info if best_tsv.lower().find(l.lower()) != -1]
            leaf = leaf[0] if leaf else "unknown"
            species = leaf
            sample_assignments[sample] = (pathogen, species, score)
        else:
            sample_assignments[sample] = ("unknown", "unknown", score)

    # --- Group samples by species ---
    species_groups = defaultdict(list)  # { (pathogen, species): [sample_names] }
    for sample, (pathogen, species, score) in sample_assignments.items():
        species_groups[(pathogen, species)].append(sample)

    # --- Write assignments TSV ---
    with open("species_assignments.tsv", "w") as f:
        f.write("sample\\tpathogen\\tspecies\\tqc_score\\tbest_dataset_file\\n")
        for sample in sorted(sample_assignments.keys()):
            pathogen, species, score = sample_assignments[sample]
            best_tsv = sample_best[sample][1] or ""
            f.write(f"{sample}\\t{pathogen}\\t{species}\\t{score}\\t{best_tsv}\\n")

    # --- Split FASTA and metadata by species ---
    # Simple FASTA parser (no BioPython dependency)
    def parse_fasta(filepath):
        sequences = {}
        current_id = None
        current_seq = []
        with open(filepath) as f:
            for line in f:
                line = line.rstrip()
                if line.startswith(">"):
                    if current_id:
                        sequences[current_id] = "\\n".join([f">{current_id}"] + current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id:
                sequences[current_id] = "\\n".join([f">{current_id}"] + current_seq)
        return sequences

    fasta_file = "${fasta}"
    sequences = parse_fasta(fasta_file)

    # Read metadata
    meta_rows = {}
    meta_header = None
    if os.path.exists(metadata_file) and Path(metadata_file).name != "NO_FILE":
        with open(metadata_file) as f:
            reader = csv.DictReader(f, delimiter="\\t")
            meta_header = reader.fieldnames
            for row in reader:
                # Try common ID columns
                row_id = row.get("strain") or row.get("name") or row.get("sample") or row.get("seqName") or ""
                if row_id:
                    meta_rows[row_id] = row

    # Write per-species files
    species_groups_output = []
    prefix = "${prefix}"

    for (pathogen, species), samples in species_groups.items():
        # Write FASTA
        fasta_out = f"{prefix}_species_{species}.fasta"
        with open(fasta_out, "w") as f:
            for sample in samples:
                if sample in sequences:
                    f.write(sequences[sample] + "\\n")

        # Write metadata
        meta_out = f"{prefix}_species_{species}.metadata.tsv"
        with open(meta_out, "w") as f:
            if meta_header:
                writer = csv.DictWriter(f, fieldnames=meta_header, delimiter="\\t")
                writer.writeheader()
                for sample in samples:
                    if sample in meta_rows:
                        writer.writerow(meta_rows[sample])
            else:
                # Minimal metadata with just sample names
                f.write("strain\\n")
                for sample in samples:
                    f.write(f"{sample}\\n")

        species_groups_output.append({
            "pathogen": pathogen,
            "species": species,
            "samples": samples,
            "fasta": str(Path(fasta_out).resolve()),
            "metadata": str(Path(meta_out).resolve())
        })

    # Write JSON summary for downstream channel construction
    with open("species_groups.json", "w") as f:
        json.dump(species_groups_output, f, indent=2)

    print(f"Species assignment complete: {len(species_groups)} group(s) detected")
    for (pathogen, species), samples in species_groups.items():
        print(f"  {pathogen}/{species}: {len(samples)} sample(s)")
    """
}
