#!/usr/bin/env python3
"""
Stage 6: Phylogenetic tree construction.

Placeholder module. The planned implementation will:
  1. Fetch contextual genomes from NCBI (accessions from Stage 2).
  2. Build a multiple sequence alignment with MAFFT.
  3. Construct a maximum-likelihood tree with IQ-TREE2 (with UFboot support).
  4. Optionally time-scale the tree with TreeTime/LSD2 when enough dates exist.

For now, the Nextclade placement tree (`nextclade.nwk`) is reused.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo

log = logging.getLogger(__name__)


def _count_tips(tree_path: Path) -> int:
    """Return the number of terminal nodes in a Newick tree."""
    try:
        tree = Phylo.read(tree_path, "newick")
        return len(tree.get_terminals())
    except Exception as exc:
        log.warning(f"Could not parse placement tree {tree_path}: {exc}")
        return 0


def run_phylogenetics(
    sample_fasta: Path,
    reference_fasta: Path,
    context_genomes: list[dict[str, Any]],
    output_dir: Path,
    lineage_result: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Return the Nextclade placement tree path as a temporary stand-in."""
    tree_path: Optional[Path] = None
    if lineage_result:
        tree_path = Path(lineage_result["nextclade"]["tree"])
        if not tree_path.exists():
            tree_path = None

    if tree_path and tree_path.exists():
        dest = output_dir / "tree.nwk"
        import shutil
        shutil.copy(tree_path, dest)
        tree_path = dest

    n_tips = _count_tips(tree_path) if tree_path else 0

    return {
        "method": "nextclade_placement_stub",
        "tree_method": "Nextclade phylogenetic placement",
        "model": None,
        "support_metric": None,
        "bootstrap_replicates": None,
        "alignment_file": None,
        "tree_file": str(tree_path) if tree_path else None,
        "sequences_in_tree": n_tips,
        "time_scaled_tree": None,
        "note": "Full ML phylogenetics with IQ-TREE2 is a planned enhancement. "
                "Currently reuses the Nextclade placement tree.",
    }
