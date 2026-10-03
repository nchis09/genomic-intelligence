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
    def nextclade_species_overrides_str = (params.nextclade_species_overrides ?: [:])
        .collect { ds, ov -> "${ds}|${ov[0]}|${ov[1]}" }.join(';')
    """
    #!/usr/bin/env python3
    import json, csv, os, sys, re
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
    # CLASSIFICATION names each Nextclade run "<sample>_<dataset-suffix>"
    # where the suffix is the dataset path minus its source segment, joined
    # with "_" (e.g. nextstrain/orthoebolavirus/bdbv -> "orthoebolavirus_bdbv",
    # community/moncla-lab/iav-h5/ha/all-clades ->
    # "moncla-lab_iav-h5_ha_all-clades"). We match each TSV filename on that
    # suffix — matching on the leaf dir name alone is ambiguous
    # (nextstrain/mpox/all-clades and community/.../all-clades share it).
    # Pathogen family is inferred from the path segment right after the
    # source prefix for official datasets (nextstrain/<virus>/...), or after
    # the maintainer segment for community datasets
    # (community/<maintainer>/<virus>/...). species is a human-readable label
    # built from everything after that; a params override can pin a clean
    # (pathogen, species) pair per dataset.
    nextclade_datasets_str = "${params.nextclade_datasets}".strip("[]")
    dataset_list = [d.strip().strip("'").strip('"') for d in nextclade_datasets_str.split(",")]

    # params.nextclade_species_overrides: map of dataset -> [pathogen, species],
    # rendered as "ds1|pathogen|species;ds2|pathogen|species;..."
    overrides_str = "${nextclade_species_overrides_str}".strip()
    overrides = {}
    if overrides_str:
        for ent in overrides_str.split(";"):
            if "|" in ent:
                ds_p, ov_p, ov_s = ent.split("|")
                overrides[ds_p.strip()] = (ov_p.strip(), ov_s.strip())

    # Build a map: dataset-suffix -> (pathogen, species, dataset)
    suffix_to_info = {}
    for ds in dataset_list:
        parts = [p for p in ds.strip("/").split("/") if p]
        if len(parts) < 2:
            continue
        suffix = "_".join(parts[1:])
        if ds in overrides:
            pathogen, species = overrides[ds]
        elif parts[0] == "community" and len(parts) > 3:
            # community/<maintainer>/<virus>/<rest...>: the virus family is
            # the segment after the maintainer name.
            pathogen = parts[2]
            remainder = parts[3:]
            species = remainder[-1] if remainder else pathogen
        elif parts[0] == "nextstrain" and len(parts) > 2:
            pathogen = parts[1]
            remainder = parts[2:]
            species = remainder[-1] if remainder else pathogen
        else:
            pathogen = parts[0]
            remainder = parts[1:]
            species = remainder[-1] if remainder else pathogen
        suffix_to_info[suffix] = (pathogen, species, ds)

    # Map each TSV file to a species based on filename containing the suffix
    tsv_to_species = {}
    for tsv_name in dataset_results.keys():
        matched = False
        for suffix in sorted(suffix_to_info.keys(), key=len, reverse=True):
            pathogen, species, ds = suffix_to_info[suffix]
            if suffix.lower() in tsv_name.lower():
                tsv_to_species[tsv_name] = (pathogen, species, ds)
                matched = True
                break
        if not matched:
            if suffix_to_info:
                first_pathogen, first_species, first_ds = list(suffix_to_info.values())[0]
                tsv_to_species[tsv_name] = (first_pathogen, first_species, first_ds)
            else:
                tsv_to_species[tsv_name] = ("unknown", "unknown", "")

    # --- Two-pass assignment: record typing, then isolate consensus ---
    #
    # A record's own best dataset is only EVIDENCE. Routing is decided per
    # isolate (records sharing a segment-suffix prefix): a typed HA record
    # decides the isolate's lineage, an NA record is the fallback, and
    # isolates with no typed HA/NA record are marked internal_only so the
    # router can screen+report them without spawning a tree build.
    #
    # qc.overallScore is a penalty score (0 = perfect placement); records
    # above params.nextclade_max_score are treated as untyped so that junk
    # or cross-lineage hits cannot create phantom groups.

    MAX_SCORE = float("${params.nextclade_max_score}")

    # Influenza record naming: <isolate>_seg<N> etc. — advisory only;
    # Nextclade decides what each record actually is.
    SEG_SUFFIX = re.compile(
        r"[_|-](seg(?:ment)?\\d+|ha|na|pb2|pb1|pa|np|mp|ns|m1|m2|ns1|nep)\$",
        re.IGNORECASE)
    # Segment encoded in a dataset path, e.g. nextstrain/flu/h3n2/ha/… or
    # community/moncla-lab/iav-h5/ha/all-clades. Non-segment datasets
    # (whole-genome pathogens) resolve to "genome".
    SEG_IN_DS = re.compile(r"/(ha|na|pb2|pb1|pa|np|mp|ns)(?:/|\$)")

    def segment_of(ds):
        m = SEG_IN_DS.search(ds)
        return m.group(1) if m else "genome"

    def isolate_of(sample):
        return SEG_SUFFIX.sub("", sample)

    # Pathogen families whose species labels are lineage-typed via HA/NA —
    # for these, a non-HA/NA decider means the isolate had no typed HA/NA.
    FLU_PATHOGENS = {"influenza", "avian_influenza", "flu"}

    record_info = {}  # sample -> typing evidence
    for sample, (score, best_tsv) in sample_best.items():
        if best_tsv and best_tsv in tsv_to_species:
            pathogen, species, ds = tsv_to_species[best_tsv]
        else:
            pathogen, species, ds = "unknown", "unknown", ""
        record_info[sample] = {
            "score": score, "tsv": best_tsv or "",
            "pathogen": pathogen, "species": species, "ds": ds,
            "segment": segment_of(ds),
            "typed": best_tsv is not None and score <= MAX_SCORE,
        }

    isolates = defaultdict(list)  # isolate_key -> [samples]
    for sample in sample_best:
        isolates[isolate_of(sample)].append(sample)

    final = {}  # sample -> final row
    for isolate, members in isolates.items():
        typed = [s for s in members if record_info[s]["typed"]]
        by_score = lambda s: record_info[s]["score"]
        ha = [s for s in typed if record_info[s]["segment"] == "ha"]
        na = [s for s in typed if record_info[s]["segment"] == "na"]
        if ha:
            decider, source = min(ha, key=by_score), "ha"
        elif na:
            decider, source = min(na, key=by_score), "na"
        elif typed:
            decider, source = min(typed, key=by_score), "record"
        else:
            decider = source = None

        if decider:
            d = record_info[decider]
            consensus = (d["pathogen"], d["species"], d["ds"])
            # NA-segment disagreement does not flag reassortment for avian
            # isolates — no avian NA datasets exist, so an avian isolate's
            # NA always best-matches a seasonal NA dataset.
            reassort = any(
                record_info[s]["typed"]
                and record_info[s]["species"] != consensus[1]
                and not (consensus[0] == "avian_influenza"
                         and record_info[s]["segment"] == "na")
                for s in members)
            internal_only = (d["pathogen"] in FLU_PATHOGENS
                             and source not in ("ha", "na"))
            for s in members:
                ri = record_info[s]
                final[s] = {
                    "isolate": isolate, "pathogen": consensus[0],
                    "species": consensus[1], "ds": consensus[2],
                    "record_species": ri["species"], "score": ri["score"],
                    "tsv": ri["tsv"], "segment": ri["segment"],
                    "lineage_source": source if s == decider else "inherited",
                    "reassortment": reassort, "internal_only": internal_only,
                }
        else:
            for s in members:
                ri = record_info[s]
                final[s] = {
                    "isolate": isolate, "pathogen": "unclassified",
                    "species": "unclassified", "ds": "",
                    "record_species": ri["species"], "score": ri["score"],
                    "tsv": ri["tsv"], "segment": ri["segment"],
                    "lineage_source": "none",
                    "reassortment": False, "internal_only": False,
                }

    # --- Group samples by isolate consensus ---
    species_groups = defaultdict(list)   # { (pathogen, species): [samples] }
    group_meta = {}                      # { (pathogen, species): {dataset, internal_only} }
    for sample, a in final.items():
        if a["pathogen"] == "unclassified":
            continue
        key = (a["pathogen"], a["species"])
        species_groups[key].append(sample)
        group_meta.setdefault(key, {"dataset": a["ds"], "internal_only": True})
    # A group is buildable as soon as ONE member isolate was decided by a
    # typed HA/NA record — internal-only isolates riding along is fine and
    # correct (they are screened context). Only all-internal groups skip.
    for key, members in species_groups.items():
        group_meta[key]["internal_only"] = all(
            final[s]["internal_only"] for s in members)

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

    # FASTA records that never appeared in any Nextclade TSV (alignment
    # failed in every dataset) are reported as unclassified too, so the
    # assignments table covers every submitted record.
    for sample in sequences:
        if sample not in final:
            final[sample] = {
                "isolate": isolate_of(sample), "pathogen": "unclassified",
                "species": "unclassified", "ds": "",
                "record_species": "unknown", "score": "",
                "tsv": "", "segment": "",
                "lineage_source": "none",
                "reassortment": False, "internal_only": False,
            }

    # --- Write assignments TSV ---
    with open("species_assignments.tsv", "w") as f:
        f.write("sample\\tisolate\\tpathogen\\tspecies\\trecord_species\\tqc_score"
                "\\tbest_dataset_file\\tsegment\\tlineage_source\\treassortment_suspected\\n")
        for sample in sorted(final.keys()):
            a = final[sample]
            f.write(f"{sample}\\t{a['isolate']}\\t{a['pathogen']}\\t{a['species']}"
                    f"\\t{a['record_species']}\\t{a['score']}\\t{a['tsv']}"
                    f"\\t{a['segment']}\\t{a['lineage_source']}"
                    f"\\t{str(a['reassortment']).lower()}\\n")

    # Read metadata
    meta_rows = {}
    meta_header = None
    if os.path.exists(metadata_file) and Path(metadata_file).name != "NO_FILE":
        with open(metadata_file) as f:
            reader = csv.DictReader(f, delimiter="\\t")
            meta_header = reader.fieldnames
            for row in reader:
                # Try common ID columns (in priority order)
                row_id = (row.get("strain") or row.get("accession") or row.get("name")
                          or row.get("sample") or row.get("sample_id") or row.get("seqName") or "")
                if row_id:
                    meta_rows[row_id] = row

    # Metadata is keyed at isolate level, but segmented genomes (influenza)
    # submit one FASTA record per segment named '<isolate>_<seg>' or
    # '<isolate>|<seg>'. Fall back to the isolate prefix when a record has
    # no exact metadata row.
    def lookup_meta(sample):
        if sample in meta_rows:
            return meta_rows[sample], False
        base = SEG_SUFFIX.sub("", sample)
        if base != sample and base in meta_rows:
            return meta_rows[base], True
        return None, False

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

        # Write metadata. Downstream augur merge keys rows on `accession`, so
        # guarantee that column exists: copy the detected ID column into
        # `accession` when the source metadata used a different name
        # (e.g. sample_id, strain, seqName).
        meta_out = f"{prefix}_species_{species}.metadata.tsv"
        with open(meta_out, "w") as f:
            if meta_header:
                out_header = list(meta_header)
                if "accession" not in out_header:
                    out_header = ["accession"] + out_header
                # All user-supplied sequences are queries; the nextstrain
                # subsample config filters on is_query=="true".
                if "is_query" not in out_header:
                    out_header.append("is_query")
                # genbank_accession holds a user-supplied real accession when
                # the pipeline accession (record id) would overwrite it
                if "genbank_accession" not in out_header:
                    out_header.append("genbank_accession")
                writer = csv.DictWriter(f, fieldnames=out_header, delimiter="\\t",
                                        extrasaction="ignore")
                writer.writeheader()
                for sample in samples:
                    row_src, is_isolate_match = lookup_meta(sample)
                    if row_src is not None:
                        row = dict(row_src)
                        # Isolate-level rows matched by prefix must carry the
                        # record-level accession — augur merge keys on the
                        # FASTA header (segment id), not the isolate id. A
                        # real GenBank accession on the row is preserved as
                        # genbank_accession rather than lost.
                        if is_isolate_match or not row.get("accession"):
                            orig = row.get("accession", "")
                            if orig and orig != sample:
                                row["genbank_accession"] = orig
                            row["accession"] = sample
                        row["is_query"] = "true"
                        writer.writerow(row)
            else:
                # Minimal metadata with just sample names
                f.write("accession\\tstrain\\tis_query\\n")
                for sample in samples:
                    f.write(f"{sample}\\t{sample}\\ttrue\\n")

        species_groups_output.append({
            "pathogen": pathogen,
            "species": species,
            "dataset": group_meta.get((pathogen, species), {}).get("dataset", ""),
            "internal_only": group_meta.get((pathogen, species), {}).get("internal_only", False),
            "samples": samples,
            "fasta": str(Path(fasta_out).resolve()),
            "metadata": str(Path(meta_out).resolve())
        })

    # Write JSON summary for downstream channel construction
    with open("species_groups.json", "w") as f:
        json.dump(species_groups_output, f, indent=2)

    print(f"Species assignment complete: {len(species_groups)} group(s) detected")
    for (pathogen, species), samples in species_groups.items():
        note = " (internal-only — screened/reported, no tree build)" \\
            if group_meta.get((pathogen, species), {}).get("internal_only") else ""
        print(f"  {pathogen}/{species}: {len(samples)} sample(s){note}")
    """
}
