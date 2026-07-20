#!/usr/bin/env python3
"""
Stage 3: Genome alignment & annotation.

This module currently relies on Nextclade's aligned output (produced during
Stage 5). A future enhancement is to run MAFFT directly against the curated
reference genome and produce a standalone alignment + GFF annotation, which is
required for the full phylogenetics module and for non-Nextclade-supported
pathogens.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional


def get_alignment_path(lineage_result: dict[str, Any]) -> Optional[Path]:
    """Return the path to the aligned FASTA produced by Nextclade, if available."""
    nextclade = lineage_result.get("nextclade", {})
    aligned = nextclade.get("aligned_fasta")
    if aligned:
        p = Path(aligned)
        if p.exists():
            return p
    return None


def run_alignment(
    sample_fasta: Path,
    reference_fasta: Path,
    output_dir: Path,
    method: str = "nextclade",
) -> dict[str, Any]:
    """Stub: return the Nextclade-aligned FASTA path when available."""
    return {
        "method": method,
        "aligned_fasta": None,
        "gff": None,
        "note": "Nextclade alignment is reused from lineage_assignment stage. "
                "MAFFT-based standalone alignment is a planned enhancement.",
    }
