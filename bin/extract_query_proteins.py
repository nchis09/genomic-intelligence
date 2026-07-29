#!/usr/bin/env python3
"""
extract_query_proteins.py

Dynamically discover UniProt accessions from pipeline outputs by:
  1. Finding phylogenetic neighbors of the query sample in the tree
  2. Finding reference tips sharing the same outbreak/strain
  3. Finding reference tips sharing identical AA mutations (per gene)
  4. Mapping GenBank (INSDC) accessions → UniProt via cross-reference search
  5. Falling back to an organism+gene UniProt search for uncovered genes
  6. Searching for reviewed (Swiss-Prot) canonical accessions per gene, with
     a genus-level cross-species fallback — required because curated
     mutagenesis/variation feature data (queried downstream by
     annotate_rbioapi.R) only exists on reviewed entries, never on the
     unreviewed isolate-specific accessions found by strategies 1-5
  7. Downloading UniProtKB TSV for discovered accessions

AA mutations (query + reference tips) are sourced from the Nextstrain
build's own nextclade.tsv (codon-aware, reference-relative substitution
calls), not from walking the Auspice tree. Nucleotide substitution count
is also read from nextclade.tsv, but is informational only — it does not
drive any discovery strategy below (only non-synonymous/AA mutations do).

Inputs:
  --auspice        Auspice JSON file
  --results_dir    Nextstrain results directory (contains nextclade.tsv)
  --query_samples  Comma-separated query sample names
  --species        Species ID (e.g., bdbv, ebov, sudv)
  --prefix         Output file prefix
  --max_neighbors  Max phylogenetic neighbors to collect (default: 5)

Outputs:
  {prefix}_discovery.tsv            - Per-gene discovery table (ref, INSDC, reason, uniprot_acc)
  {prefix}_all_accessions.txt       - Combined unique UniProt accessions (for UniprotR)
  {prefix}_uniprot_download.tsv     - UniProtKB TSV export (for UniProtExtractR)
  {prefix}_query_mutations.tsv      - AA mutations per query sample
  {prefix}_query_proteins.fasta     - Reconstructed protein sequences
  {prefix}_discovery_summary.json   - Discovery stats and sample metadata
"""

import argparse
import csv
import glob
import json
import re
import sys
import time
import urllib.request
import urllib.error
import urllib.parse


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--auspice", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--query_samples", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--prefix", default="query")
    parser.add_argument("--max_neighbors", type=int, default=5)
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------

def get_all_tips(node, tips=None):
    """Collect all tip nodes from the tree."""
    if tips is None:
        tips = []
    if not node.get("children"):
        tips.append(node)
    else:
        for child in node["children"]:
            get_all_tips(child, tips)
    return tips


def find_tip(node, name):
    """Find a tip node by name."""
    if not node.get("children"):
        return node if node["name"] == name else None
    for child in node["children"]:
        result = find_tip(child, name)
        if result:
            return result
    return None


def find_path(node, target, path):
    """Find the path from root to a target tip, accumulating nodes."""
    path.append(node)
    if node["name"] == target:
        return True
    for child in node.get("children", []):
        if find_path(child, target, path):
            return True
    path.pop()
    return False


def get_value(node_attrs, key):
    """Extract value from Auspice node_attrs (handles dict or scalar)."""
    val = node_attrs.get(key)
    if isinstance(val, dict):
        return val.get("value")
    return val


def get_insdc(tip):
    """Extract INSDC (GenBank) accession from a tip node."""
    na = tip.get("node_attrs", {})
    insdc = na.get("INSDC_accession", {})
    if isinstance(insdc, dict):
        val = insdc.get("value", "")
    else:
        val = insdc or ""
    # Strip version suffix (.1, .2) for UniProt xref search
    return val.split(".")[0] if val else ""


def load_nextclade_mutations(results_dir):
    """
    Load per-sample reference-relative AA mutations from the Nextstrain
    build's own nextclade.tsv (codon-aware, reference-relative substitution
    calls — same reference/annotation as the tree).

    Returns:
      gene_muts_by_sample : { seqName: { gene: [mutation, ...] } }
      nuc_counts_by_sample: { seqName: int }  (informational only)
    """
    nc_files = glob.glob(f"{results_dir}/**/nextclade.tsv", recursive=True)
    gene_muts_by_sample = {}
    nuc_counts_by_sample = {}
    if not nc_files:
        print(f"  WARNING: no nextclade.tsv found under {results_dir}", file=sys.stderr)
        return gene_muts_by_sample, nuc_counts_by_sample

    with open(nc_files[0]) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row.get("seqName", "")
            if not name:
                continue
            gene_muts = {}
            for entry in (row.get("aaSubstitutions") or "").split(","):
                entry = entry.strip()
                if not entry or ":" not in entry:
                    continue
                gene, mut = entry.split(":", 1)
                gene_muts.setdefault(gene, []).append(mut)
            gene_muts_by_sample[name] = gene_muts
            try:
                nuc_counts_by_sample[name] = int(row.get("totalSubstitutions") or 0)
            except ValueError:
                nuc_counts_by_sample[name] = 0
    return gene_muts_by_sample, nuc_counts_by_sample


# ---------------------------------------------------------------------------
# Discovery strategies
# ---------------------------------------------------------------------------

def discover_phylo_neighbors(tree, query_name, all_ref_tips, max_n=5):
    """
    Find the N closest reference tips by walking up from the query tip
    through successively larger ancestral subtrees.
    """
    path = []
    find_path(tree, query_name, path)

    seen = set()
    neighbors = []
    # Walk up from parent toward root, collecting tips from sibling subtrees
    for i in range(len(path) - 2, -1, -1):
        ancestor = path[i]
        subtree_tips = get_all_tips(ancestor)
        for t in subtree_tips:
            if t["name"] != query_name and t["name"] not in seen:
                seen.add(t["name"])
                neighbors.append(t)
        if len(neighbors) >= max_n:
            break

    return neighbors[:max_n]


def discover_outbreak_matches(query_tip, all_ref_tips):
    """Find reference tips sharing the same outbreak value."""
    query_outbreak = get_value(query_tip.get("node_attrs", {}), "outbreak")
    if not query_outbreak:
        return []
    matches = []
    for t in all_ref_tips:
        ref_outbreak = get_value(t.get("node_attrs", {}), "outbreak")
        if ref_outbreak == query_outbreak:
            matches.append(t)
    return matches


def discover_mutation_matches(query_gene_muts, all_ref_tips, nc_gene_muts_by_sample):
    """
    For each gene and its mutations in the query, find reference tips
    that share the exact same mutation IN THE SAME GENE. Reference tip
    mutations are looked up from the same nextclade.tsv-derived dict as
    the query, so both sides are compared on the same reference-relative
    basis.
    Returns: list of (tip, gene, mutation_label)
    """
    matches = []
    for gene, mut_list in query_gene_muts.items():
        for mut_str in mut_list:
            for t in all_ref_tips:
                ref_gene_muts = nc_gene_muts_by_sample.get(t["name"], {}).get(gene, [])
                if mut_str in ref_gene_muts:
                    matches.append((t, gene, f"{gene}:{mut_str}"))
    return matches


# ---------------------------------------------------------------------------
# UniProt helpers
# ---------------------------------------------------------------------------

def genbank_to_uniprot(genbank_ids):
    """
    Map GenBank (INSDC) accessions to UniProt accessions using xref search.
    Returns: dict { genbank_id: [ { accession, gene, protein_name } ] }
    """
    results = {}
    for gid in genbank_ids:
        if not gid:
            continue
        # Query UniProt for entries cross-referenced to this EMBL/GenBank ID
        query = f"xref:embl-{gid}"
        fields = "accession,gene_names,protein_name,organism_name"
        url = (
            f"https://rest.uniprot.org/uniprotkb/search?"
            f"query={urllib.parse.quote(query)}"
            f"&format=tsv&fields={fields}&size=50"
        )
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PGIRL-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                content = resp.read().decode("utf-8")
            lines = content.strip().split("\n")
            entries = []
            if len(lines) > 1:
                header = lines[0].split("\t")
                for line in lines[1:]:
                    cols = line.split("\t")
                    row = dict(zip(header, cols))
                    entries.append({
                        "accession": row.get("Entry", ""),
                        "gene": row.get("Gene Names", ""),
                        "protein_name": row.get("Protein names", ""),
                        "organism": row.get("Organism", ""),
                    })
            results[gid] = entries
            if entries:
                print(f"    xref:{gid} → {len(entries)} UniProt entries", file=sys.stderr)
            else:
                print(f"    xref:{gid} → no UniProt entries found", file=sys.stderr)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"    xref:{gid} → ERROR: {e}", file=sys.stderr)
            results[gid] = []
        time.sleep(0.3)  # Rate limiting
    return results


GENE_SEARCH_MAP = {
    "GP": "GP", "GP_003": "GP", "NP": "NP",
    "VP35": "VP35", "VP40": "VP40", "VP30": "VP30",
    "VP24": "VP24", "L": "L",
}

ORGANISM_MAP = {
    "bdbv": "Bundibugyo ebolavirus", "ebov": "Zaire ebolavirus",
    "sudv": "Sudan ebolavirus", "tafv": "Tai Forest ebolavirus",
    "restv": "Reston ebolavirus",
}

# Genus/family-level term used to fall back to a related species' reviewed
# UniProt entry when the exact species has no reviewed Swiss-Prot coverage.
GENUS_FALLBACK_NAME = "ebolavirus"


def search_uniprot_by_organism(organism_name, genome_genes):
    """
    Fallback: search UniProt by organism name to find protein accessions
    for each gene. Used when xref mapping returns poor coverage.
    Returns: dict { gene: [ { accession, gene, protein_name, organism } ] }
    """
    gene_search_map = GENE_SEARCH_MAP
    results = {}
    # Clean organism name for search
    org_query = organism_name.strip()
    if not org_query:
        return results

    for gene in genome_genes:
        search_gene = gene_search_map.get(gene, gene)
        query = f'(organism_name:"{org_query}") AND (gene:{search_gene})'
        fields = "accession,gene_names,protein_name,organism_name"
        url = (
            f"https://rest.uniprot.org/uniprotkb/search?"
            f"query={urllib.parse.quote(query)}"
            f"&format=tsv&fields={fields}&size=10"
        )
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PGIRL-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                content = resp.read().decode("utf-8")
            lines = content.strip().split("\n")
            entries = []
            if len(lines) > 1:
                header = lines[0].split("\t")
                for line in lines[1:]:
                    cols = line.split("\t")
                    row = dict(zip(header, cols))
                    entries.append({
                        "accession": row.get("Entry", ""),
                        "gene": row.get("Gene Names", ""),
                        "protein_name": row.get("Protein names", ""),
                        "organism": row.get("Organism", ""),
                    })
            results[gene] = entries
            if entries:
                print(f"    organism+gene:{gene} → {len(entries)} UniProt entries", file=sys.stderr)
            else:
                print(f"    organism+gene:{gene} → no entries", file=sys.stderr)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"    organism+gene:{gene} → ERROR: {e}", file=sys.stderr)
            results[gene] = []
        time.sleep(0.3)
    return results


def search_uniprot_reviewed_canonical(organism_name, genus_name, genome_genes):
    """
    Find reviewed (Swiss-Prot) canonical UniProt accessions per gene.

    Curated feature annotations used by rbioapi's mutagenesis/variation
    endpoints (EBI Proteins API) only exist on reviewed Swiss-Prot entries,
    never on the unreviewed (TrEMBL) isolate-specific accessions produced
    by the xref/organism discovery strategies above. This strategy runs
    unconditionally for every gene to surface those curated accessions.

    First tries the exact species organism name; if no reviewed entry
    exists for that species (common for e.g. Bundibugyo/Tai Forest/Reston
    ebolavirus, which have sparse Swiss-Prot coverage), falls back to a
    broader genus-level search and picks up a reviewed entry from a
    related species instead (tagged as cross-species).

    NOTE: cross-species entries may use different residue numbering than
    the reference sequence used for phylogenetic mutation calling, so
    query_mutation_match position overlaps for those rows are best-effort
    and not guaranteed to be biologically equivalent positions.

    Returns: dict { gene: [ { accession, gene, protein_name, organism,
                               cross_species } ] }
    """
    def _search(org_query, gene_query):
        query = f'reviewed:true AND (organism_name:"{org_query}") AND (gene:{gene_query})'
        fields = "accession,gene_names,protein_name,organism_name"
        url = (
            f"https://rest.uniprot.org/uniprotkb/search?"
            f"query={urllib.parse.quote(query)}"
            f"&format=tsv&fields={fields}&size=10"
        )
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PGIRL-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                content = resp.read().decode("utf-8")
            lines = content.strip().split("\n")
            entries = []
            if len(lines) > 1:
                header = lines[0].split("\t")
                for line in lines[1:]:
                    cols = line.split("\t")
                    row = dict(zip(header, cols))
                    entries.append({
                        "accession": row.get("Entry", ""),
                        "gene": row.get("Gene Names", ""),
                        "protein_name": row.get("Protein names", ""),
                        "organism": row.get("Organism", ""),
                    })
            return entries
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"    reviewed_canonical ERROR: {e}", file=sys.stderr)
            return []

    results = {}
    org_query = (organism_name or "").strip()
    genus_query = (genus_name or "").strip()

    for gene in genome_genes:
        search_gene = GENE_SEARCH_MAP.get(gene, gene)
        entries = []
        cross_species = False

        if org_query:
            entries = _search(org_query, search_gene)

        if not entries and genus_query:
            entries = _search(genus_query, search_gene)
            cross_species = True

        for e in entries:
            e["cross_species"] = cross_species

        results[gene] = entries
        if entries:
            label = "cross-species fallback" if cross_species else "species-specific"
            print(f"    reviewed_canonical:{gene} → {len(entries)} entries ({label})", file=sys.stderr)
        else:
            print(f"    reviewed_canonical:{gene} → no reviewed entries found (species or genus)", file=sys.stderr)
        time.sleep(0.3)

    return results


def download_uniprot_tsv(accessions):
    """
    Download UniProtKB TSV for the given accessions from the REST API.
    Columns match what UniProtExtractR expects.
    """
    fields = [
        "accession", "reviewed", "id", "protein_name", "gene_names",
        "organism_name", "length",
        "cc_subcellular_location", "ft_transmem", "ft_domain", "ft_signal",
        "ft_motif", "protein_families", "cc_pathway", "ft_dna_bind",
        "cc_disease", "cc_function", "ft_binding", "cc_ptm",
        "ft_carbohyd", "ft_disulfid", "ft_lipid", "ft_region",
        "go_p", "go_f", "go_c", "go_id",
    ]

    unique_acc = sorted(set(accessions))
    if not unique_acc:
        return None

    # UniProt's search API rejects overly long OR-combined queries
    # (observed HTTP 400 once the accession count grew past ~100 after
    # adding the reviewed-canonical discovery strategy). Batch requests
    # to stay well under that limit, then concatenate results (keeping
    # only the first header line).
    batch_size = 90
    batches = [unique_acc[i:i + batch_size] for i in range(0, len(unique_acc), batch_size)]

    combined_lines = None
    for batch in batches:
        acc_query = " OR ".join(batch)
        query = f"accession:({acc_query})"
        url = (
            f"https://rest.uniprot.org/uniprotkb/search?"
            f"query={urllib.parse.quote(query)}"
            f"&format=tsv"
            f"&fields={','.join(fields)}"
            f"&size=500"
        )
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PGIRL-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                content = resp.read().decode("utf-8")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"  WARNING: UniProt TSV download failed for a batch of {len(batch)} accessions: {e}", file=sys.stderr)
            continue

        lines = content.strip().split("\n")
        if combined_lines is None:
            combined_lines = lines
        else:
            combined_lines.extend(lines[1:])  # skip repeated header
        time.sleep(0.3)

    if combined_lines is None:
        return None
    return "\n".join(combined_lines) + "\n"


# ---------------------------------------------------------------------------
# Mutation / sequence helpers
# ---------------------------------------------------------------------------

def parse_mutation(mut_str):
    """Parse mutation string like 'Y387H' into (ref_aa, position, alt_aa)."""
    match = re.match(r"([A-Z*])(\d+)([A-Z*])", mut_str)
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None, None, None


def apply_mutations(seq, mutations):
    """Apply AA mutations to a reference sequence, return mutated sequence."""
    seq_list = list(seq)
    for mut_str in mutations:
        ref_aa, pos, alt_aa = parse_mutation(mut_str)
        if pos is not None and pos <= len(seq_list):
            seq_list[pos - 1] = alt_aa
    return "".join(seq_list)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # --- Load Auspice JSON ---
    print(f"Loading Auspice JSON: {args.auspice}", file=sys.stderr)
    with open(args.auspice) as f:
        auspice = json.load(f)

    tree = auspice["tree"]
    meta = auspice.get("meta", {})
    root_seqs = auspice.get("root_sequence", {})
    genome_genes = [g for g in meta.get("genome_annotations", {}).keys() if g != "nuc"]

    # --- Load reference-relative AA mutations (query + all reference tips) ---
    print(f"Loading nextclade.tsv mutations from: {args.results_dir}", file=sys.stderr)
    nc_gene_muts, nc_nuc_counts = load_nextclade_mutations(args.results_dir)
    print(f"  {len(nc_gene_muts)} samples with reference-relative mutations loaded", file=sys.stderr)

    # --- Collect all reference tips ---
    query_names = [s.strip() for s in args.query_samples.split(",") if s.strip()]
    all_tips = get_all_tips(tree)
    all_ref_tips = [t for t in all_tips if t["name"] not in query_names]
    print(f"Tree: {len(all_tips)} tips, {len(all_ref_tips)} references, {len(query_names)} query", file=sys.stderr)

    # --- Process each query sample ---
    all_mutations = []
    all_discoveries = []  # rows for discovery.tsv
    all_uniprot_accessions = set()
    reconstructed_seqs = {}
    sample_summaries = []

    for sample_name in query_names:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"Processing query sample: {sample_name}", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        tip = find_tip(tree, sample_name)
        if not tip:
            print(f"  WARNING: '{sample_name}' not found in tree", file=sys.stderr)
            continue

        na = tip.get("node_attrs", {})
        sample_meta = {k: get_value(na, k) for k in na}

        # --- Query mutations (reference-relative, from nextclade.tsv) ---
        if sample_name not in nc_gene_muts:
            print(f"  WARNING: '{sample_name}' not found in nextclade.tsv; no reference-relative mutations available", file=sys.stderr)
        gene_mutations = nc_gene_muts.get(sample_name, {})
        nuc_count = nc_nuc_counts.get(sample_name, 0)  # informational only
        print(f"  Mutations: {sum(len(v) for v in gene_mutations.values())} AA across {len(gene_mutations)} genes, {nuc_count} nuc (informational)", file=sys.stderr)

        for gene in sorted(gene_mutations.keys()):
            for mut_str in gene_mutations[gene]:
                ref_aa, pos, alt_aa = parse_mutation(mut_str)
                all_mutations.append({
                    "sample": sample_name, "gene": gene,
                    "position": pos, "ref_aa": ref_aa, "alt_aa": alt_aa,
                    "mutation_label": f"{gene}:{mut_str}",
                })

        # --- Strategy 1: Phylogenetic neighbors ---
        print(f"\n  [Strategy 1] Phylogenetic neighbors (max {args.max_neighbors}):", file=sys.stderr)
        phylo_neighbors = discover_phylo_neighbors(tree, sample_name, all_ref_tips, args.max_neighbors)
        for t in phylo_neighbors:
            insdc = get_insdc(t)
            print(f"    {t['name']} (INSDC={insdc})", file=sys.stderr)
            all_discoveries.append({
                "sample": sample_name, "gene": "*",
                "ref_tip": t["name"], "insdc": insdc,
                "reason": "phylogenetic_neighbor",
            })

        # --- Strategy 2: Outbreak matching ---
        print(f"\n  [Strategy 2] Outbreak matching:", file=sys.stderr)
        outbreak_matches = discover_outbreak_matches(tip, all_ref_tips)
        outbreak_insdc = set()
        for t in outbreak_matches:
            insdc = get_insdc(t)
            if insdc and t["name"] not in [n["name"] for n in phylo_neighbors]:
                outbreak_insdc.add(insdc)
                all_discoveries.append({
                    "sample": sample_name, "gene": "*",
                    "ref_tip": t["name"], "insdc": insdc,
                    "reason": f"same_outbreak:{sample_meta.get('outbreak','')}",
                })
        print(f"    {len(outbreak_matches)} tips share outbreak, {len(outbreak_insdc)} new INSDC accessions", file=sys.stderr)

        # --- Strategy 3: Mutation matching (per gene) ---
        print(f"\n  [Strategy 3] Mutation matching (per gene):", file=sys.stderr)
        mut_matches = discover_mutation_matches(gene_mutations, all_ref_tips, nc_gene_muts)
        mut_match_insdc = set()
        seen_match = set()
        for t, gene, mut_label in mut_matches:
            key = (t["name"], gene, mut_label)
            if key in seen_match:
                continue
            seen_match.add(key)
            insdc = get_insdc(t)
            if insdc:
                mut_match_insdc.add(insdc)
                all_discoveries.append({
                    "sample": sample_name, "gene": gene,
                    "ref_tip": t["name"], "insdc": insdc,
                    "reason": f"shared_mutation:{mut_label}",
                })
        print(f"    {len(seen_match)} gene:mutation matches, {len(mut_match_insdc)} unique INSDC accessions", file=sys.stderr)

        # --- Collect all unique INSDC accessions for this sample ---
        sample_insdc = set()
        for d in all_discoveries:
            if d["sample"] == sample_name and d["insdc"]:
                sample_insdc.add(d["insdc"])
        print(f"\n  Total unique INSDC accessions to map: {len(sample_insdc)}", file=sys.stderr)

        # --- Strategy 4: GenBank → UniProt xref mapping ---
        print(f"\n  [Strategy 4] GenBank → UniProt xref mapping:", file=sys.stderr)
        xref_results = genbank_to_uniprot(sorted(sample_insdc))

        # Enrich discovery rows with UniProt accessions
        for d in all_discoveries:
            if d["sample"] == sample_name and d["insdc"]:
                entries = xref_results.get(d["insdc"], [])
                # If discovery is gene-specific, filter to that gene
                if d["gene"] != "*":
                    gene_lower = d["gene"].lower()
                    matching = [e for e in entries if gene_lower in e.get("gene", "").lower()]
                    d["uniprot_accessions"] = ",".join(e["accession"] for e in matching) if matching else ""
                    for e in matching:
                        all_uniprot_accessions.add(e["accession"])
                else:
                    d["uniprot_accessions"] = ",".join(e["accession"] for e in entries)
                    for e in entries:
                        all_uniprot_accessions.add(e["accession"])

        # --- Strategy 5 (Fallback): Search by organism + gene ---
        # If xref mapping yielded few accessions, fallback to organism-based search
        genes_covered = set()
        for d in all_discoveries:
            if d["sample"] == sample_name and d.get("uniprot_accessions"):
                if d["gene"] != "*":
                    genes_covered.add(d["gene"])
                else:
                    # Wildcard entries cover all genes
                    if d["uniprot_accessions"]:
                        genes_covered.update(genome_genes)

        uncovered_genes = [g for g in genome_genes if g not in genes_covered]

        # Common patterns: "Bundibugyo virus", "Ebola virus", etc.
        organism_name = ORGANISM_MAP.get(args.species, "")

        if uncovered_genes and organism_name:
            print(f"\n  [Strategy 5] Fallback: search UniProt by organism '{organism_name}' for {len(uncovered_genes)} uncovered genes:", file=sys.stderr)
            org_results = search_uniprot_by_organism(organism_name, uncovered_genes)
            for gene, entries in org_results.items():
                for e in entries:
                    all_uniprot_accessions.add(e["accession"])
                    all_discoveries.append({
                        "sample": sample_name, "gene": gene,
                        "ref_tip": "organism_search",
                        "insdc": "",
                        "reason": f"organism_fallback:{organism_name}:{gene}",
                        "uniprot_accessions": e["accession"],
                    })

        # --- Strategy 6: Reviewed-canonical accessions for mutagenesis/variation ---
        # Runs unconditionally for every gene: curated mutagenesis/variant
        # feature data (queried later by annotate_rbioapi.R) only exists on
        # reviewed Swiss-Prot entries, never on the unreviewed isolate-
        # specific accessions found by strategies 1-5 above.
        print(f"\n  [Strategy 6] Reviewed-canonical UniProt accessions (species='{organism_name or 'unknown'}', genus fallback='{GENUS_FALLBACK_NAME}'):", file=sys.stderr)
        reviewed_results = search_uniprot_reviewed_canonical(organism_name, GENUS_FALLBACK_NAME, genome_genes)
        n_reviewed = 0
        n_reviewed_excluded = 0
        for gene, entries in reviewed_results.items():
            for e in entries:
                is_xspecies = bool(e.get("cross_species"))
                reason = (
                    f"reviewed_canonical_crossspecies_excluded:{e.get('organism','')}:{gene}"
                    if is_xspecies
                    else f"reviewed_canonical:{e.get('organism','')}:{gene}"
                )
                # Cross-species fallback hits are recorded in discovery.tsv for
                # transparency/debugging, but deliberately NOT added to the
                # accession set queried for annotation output. Substituting a
                # related species' entry here is what caused near-identical
                # rows to appear across different species' uniprotr/rbioapi
                # results (e.g. the same Zaire ebolavirus accession showing up
                # for bdbv and sudv alike). A gene with no species-specific
                # reviewed accession is left uncovered rather than faked.
                if is_xspecies:
                    n_reviewed_excluded += 1
                else:
                    n_reviewed += 1
                    all_uniprot_accessions.add(e["accession"])
                all_discoveries.append({
                    "sample": sample_name, "gene": gene,
                    "ref_tip": "reviewed_canonical_search",
                    "insdc": "",
                    "reason": reason,
                    "uniprot_accessions": "" if is_xspecies else e["accession"],
                })

        print(f"\n  Total UniProt accessions after all strategies: {len(all_uniprot_accessions)} "
              f"({n_reviewed} reviewed-canonical, {n_reviewed_excluded} cross-species excluded)", file=sys.stderr)

        # --- Reconstruct protein sequences ---
        for gene in genome_genes:
            if gene in root_seqs:
                muts_for_gene = gene_mutations.get(gene, [])
                seq = apply_mutations(root_seqs[gene], muts_for_gene)
                reconstructed_seqs[f"{sample_name}_{gene}"] = seq

        # --- Build sample summary ---
        sample_summaries.append({
            "sample": sample_name,
            "species": args.species,
            "title": meta.get("title", ""),
            "outbreak": sample_meta.get("outbreak", ""),
            "host": sample_meta.get("host", ""),
            "country": sample_meta.get("country", ""),
            "strain": sample_meta.get("strain", ""),
            "date": sample_meta.get("date", ""),
            "genome_coverage": sample_meta.get("genome_coverage", ""),
            "nuc_substitution_count": nuc_count,
            "aa_mutation_count": sum(len(v) for v in gene_mutations.values()),
            "genes_with_mutations": sorted(gene_mutations.keys()),
            "phylo_neighbors": len(phylo_neighbors),
            "outbreak_matches": len(outbreak_matches),
            "mutation_matches": len(seen_match),
            "unique_insdc_discovered": len(sample_insdc),
            "unique_uniprot_discovered": len(all_uniprot_accessions),
            "reviewed_canonical_accessions": n_reviewed,
        })

    # -----------------------------------------------------------------------
    # Write output files
    # -----------------------------------------------------------------------

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Writing outputs (prefix={args.prefix})", file=sys.stderr)
    print(f"Total UniProt accessions discovered: {len(all_uniprot_accessions)}", file=sys.stderr)

    # 1. Discovery table
    disc_file = f"{args.prefix}_discovery.tsv"
    disc_fields = ["sample", "gene", "ref_tip", "insdc", "reason", "uniprot_accessions"]
    with open(disc_file, "w") as f:
        writer = csv.DictWriter(f, fieldnames=disc_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_discoveries)
    print(f"  discovery.tsv: {len(all_discoveries)} rows → {disc_file}", file=sys.stderr)

    # 2. All accessions (for UniprotR)
    acc_file = f"{args.prefix}_all_accessions.txt"
    sorted_acc = sorted(all_uniprot_accessions)
    with open(acc_file, "w") as f:
        for acc in sorted_acc:
            f.write(acc + "\n")
    print(f"  all_accessions.txt: {len(sorted_acc)} accessions → {acc_file}", file=sys.stderr)

    # 3. Query mutations TSV
    mut_file = f"{args.prefix}_query_mutations.tsv"
    with open(mut_file, "w") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "gene", "position", "ref_aa", "alt_aa", "mutation_label"], delimiter="\t")
        writer.writeheader()
        writer.writerows(all_mutations)
    print(f"  query_mutations.tsv: {len(all_mutations)} mutations → {mut_file}", file=sys.stderr)

    # 4. Download UniProtKB TSV (for UniProtExtractR)
    tsv_file = f"{args.prefix}_uniprot_download.tsv"
    if sorted_acc:
        print(f"  Downloading UniProtKB TSV for {len(sorted_acc)} accessions...", file=sys.stderr)
        tsv_content = download_uniprot_tsv(sorted_acc)
        if tsv_content:
            with open(tsv_file, "w") as f:
                f.write(tsv_content)
            row_count = len(tsv_content.strip().split("\n")) - 1
            print(f"  uniprot_download.tsv: {row_count} entries → {tsv_file}", file=sys.stderr)
        else:
            with open(tsv_file, "w") as f:
                f.write("Entry\n")
            print(f"  uniprot_download.tsv: empty placeholder → {tsv_file}", file=sys.stderr)
    else:
        with open(tsv_file, "w") as f:
            f.write("Entry\n")
        print(f"  uniprot_download.tsv: no accessions discovered → {tsv_file}", file=sys.stderr)

    # 5. Reconstructed protein sequences FASTA
    fasta_file = f"{args.prefix}_query_proteins.fasta"
    with open(fasta_file, "w") as f:
        for header, seq in sorted(reconstructed_seqs.items()):
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    print(f"  query_proteins.fasta: {len(reconstructed_seqs)} sequences → {fasta_file}", file=sys.stderr)

    # 6. Discovery summary JSON
    summary_file = f"{args.prefix}_discovery_summary.json"
    summary = {
        "species": args.species,
        "title": meta.get("title", ""),
        "genome_genes": genome_genes,
        "total_ref_tips": len(all_ref_tips),
        "total_uniprot_accessions": len(all_uniprot_accessions),
        "accessions": sorted_acc,
        "samples": sample_summaries,
    }
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  discovery_summary.json → {summary_file}", file=sys.stderr)

    print("\nDone.", file=sys.stderr)


if __name__ == "__main__":
    main()
