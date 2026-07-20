#!/usr/bin/env python3
"""Generate a synthetic Newick tree for intelligence-engine development.

This is a placeholder for the real IQ-TREE/TreeTime output. It selects
contextual genomes from the curated genome_metadata.csv, adds the focal
sample, and writes a dated phylogeny to
output/bioinformatics/EBOV-UGA-2027-001/tree.nwk.
The tree is random but tip labels match real accessions so the engine's
tip-annotation and R figures can be developed before the bioinformatics
pipeline is finished.
"""

import argparse
import json
import random
import sys
from pathlib import Path

import dendropy
import pandas as pd

# Ensure project root is importable when running as a script
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from intelligence_engine.evidence_integration.engine import (  # noqa: E402
    _match_lineage,
    _safe_str,
)


def _lineage_id_for(lineage: str, lineages_df: pd.DataFrame) -> str:
    row = _match_lineage(lineage, lineages_df)
    return _safe_str(row.get("lineage_id")) if row is not None else ""


def generate_synthetic_tree(
    bio_output_path: str,
    genome_metadata_csv: str,
    lineages_csv: str,
    output_tree_path: str,
    n_contextual: int = 30,
    random_seed: int = 42,
) -> str:
    """Build a synthetic tree and write it as Newick. Returns the written path."""
    random.seed(random_seed)

    with open(bio_output_path) as f:
        bio = json.load(f)
    stage9 = bio.get("stage9_normalised_output", bio)
    sample_id = _safe_str(stage9.get("sample_id"))
    sample_metadata = stage9.get("metadata", {})
    sample_country = _safe_str(sample_metadata.get("country"))
    sample_date = _safe_str(sample_metadata.get("collection_date"))
    lineage = _safe_str(stage9.get("lineage"))

    genome_metadata = pd.read_csv(genome_metadata_csv, low_memory=False)
    lineages = pd.read_csv(lineages_csv, low_memory=False)

    lineage_id = _lineage_id_for(lineage, lineages)
    if not lineage_id:
        raise ValueError(f"Could not resolve lineage_id for lineage {lineage!r}")

    # Select contextual genomes for the same lineage, preferring diversity in country/date
    subset = genome_metadata[genome_metadata["lineage_id"].astype(str) == lineage_id].copy()
    if subset.empty:
        raise ValueError(f"No genomes found for lineage_id {lineage_id!r}")

    subset["collection_year"] = pd.to_numeric(
        subset["collection_date"].astype(str).str[:4], errors="coerce"
    )
    # Pick up to n_contextual/2 diverse countries, then backfill with other records
    country_counts = subset["collection_country"].value_counts()
    selected = []
    for country, _ in country_counts.items():
        if len(selected) >= n_contextual:
            break
        rows = subset[subset["collection_country"] == country]
        if rows.empty:
            continue
        # prefer a record with a valid date
        dated = rows.dropna(subset=["collection_year"])
        pick = dated.sample(1, random_state=random_seed) if not dated.empty else rows.sample(1, random_state=random_seed)
        selected.append(pick.iloc[0])

    # Backfill if we have fewer than n_contextual
    if len(selected) < n_contextual:
        used = {s["genome_accession"] for s in selected}
        remaining = subset[~subset["genome_accession"].isin(used)]
        n_needed = n_contextual - len(selected)
        if not remaining.empty:
            extras = remaining.sample(min(n_needed, len(remaining)), random_state=random_seed)
            selected.extend([r for _, r in extras.iterrows()])

    # Tip labels: genome_accession for contextual; sample_id for focal sample
    tip_labels = [str(s["genome_accession"]) for s in selected]
    if sample_id in tip_labels:
        tip_labels.remove(sample_id)
    tip_labels.append(sample_id)

    taxa = dendropy.TaxonNamespace(tip_labels)
    tree = dendropy.model.coalescent.pure_kingman_tree(
        taxon_namespace=taxa,
        pop_size=1.0,
        rng=random,
    )

    # Scale branch lengths so the total tree height spans ~13 years
    # (from the Makona MRCA ~2013 to the focal sample in 2027). This makes
    # the synthetic tree roughly interpretable on a time axis.
    leaf_distances = [leaf.distance_from_root() for leaf in tree.leaf_nodes()]
    max_tip_distance = max(leaf_distances) if leaf_distances else 1.0
    if max_tip_distance == 0:
        max_tip_distance = 1.0
    scale = 13.5 / max_tip_distance
    for edge in tree.postorder_edge_iter():
        if edge.length is not None:
            edge.length *= scale

    # Write Newick
    output_path = Path(output_tree_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(path=str(output_path), schema="newick", suppress_rooting=True)

    # Also write a small metadata table for the tree tips (useful for R / debugging)
    meta_rows = []
    for s in selected:
        meta_rows.append(
            {
                "tip_label": s["genome_accession"],
                "genome_accession": s["genome_accession"],
                "collection_country": _safe_str(s.get("collection_country")),
                "collection_date": _safe_str(s.get("collection_date")),
                "collection_year": pd.to_numeric(s.get("collection_year"), errors="coerce"),
                "host": _safe_str(s.get("host")),
                "lineage_id": lineage_id,
            }
        )
    meta_rows.append(
        {
            "tip_label": sample_id,
            "genome_accession": sample_id,
            "collection_country": sample_country,
            "collection_date": sample_date,
            "collection_year": 2027,
            "host": _safe_str(sample_metadata.get("host")),
            "lineage_id": lineage_id,
        }
    )
    tips_meta_path = output_path.with_suffix(".tips_metadata.csv")
    pd.DataFrame(meta_rows).to_csv(tips_meta_path, index=False)

    return str(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Generate a synthetic Newick tree for development."
    )
    parser.add_argument(
        "--bio-output",
        default="intelligence_engine/bioinformatics/templates/bioinformatics_output_template.json",
        help="Path to the bioinformatics output JSON.",
    )
    parser.add_argument(
        "--genome-metadata",
        default="database/exports/genome_metadata.csv",
        help="Path to the genome metadata CSV.",
    )
    parser.add_argument(
        "--lineages",
        default="database/exports/lineages.csv",
        help="Path to the lineages CSV.",
    )
    parser.add_argument(
        "--output",
        default="output/bioinformatics/EBOV-UGA-2027-001/tree.nwk",
        help="Where to write the synthetic Newick tree.",
    )
    parser.add_argument(
        "--n-contextual",
        type=int,
        default=30,
        help="Number of contextual genomes to include.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    args = parser.parse_args()

    path = generate_synthetic_tree(
        bio_output_path=args.bio_output,
        genome_metadata_csv=args.genome_metadata,
        lineages_csv=args.lineages,
        output_tree_path=args.output,
        n_contextual=args.n_contextual,
        random_seed=args.seed,
    )
    print(f"Wrote synthetic tree to: {path}")


if __name__ == "__main__":
    main()
