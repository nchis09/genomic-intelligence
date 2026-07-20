#!/usr/bin/env python3
"""
Stage 8: Comparative genomics.

Placeholder module. The planned implementation will compute:
  - Gene presence / truncation vs reference
  - Leader/trailer and gene-order integrity
  - GC content
  - dN/dS ratio (via codeml or similar)

For now, only basic sequence statistics are provided.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio import SeqIO


def _gc_content(sequence: str) -> float:
    seq = sequence.upper()
    gc = seq.count("G") + seq.count("C")
    total = len(seq)
    return round(100.0 * gc / total, 2) if total else 0.0


def run_comparative_genomics(
    sample_fasta: Path,
    reference_fasta: Path,
    gene_coordinates: dict[str, Any],
) -> dict[str, Any]:
    """Compute minimal comparative statistics."""
    sample_seq = str(SeqIO.read(sample_fasta, "fasta").seq)
    return {
        "method": "basic",
        "gene_content": "annotation pending",
        "gc_content": _gc_content(sample_seq),
        "dn_ds": None,
        "selection_pressure": None,
        "note": "Full comparative genomics (gene truncation, dN/dS) is a planned enhancement.",
    }
