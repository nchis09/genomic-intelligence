#!/usr/bin/env python3
"""
Stage 4: Variant calling.

The variant caller is currently powered by Nextclade's mutation calls (from
Stage 5), but the emphasis is deliberately on **amino-acid variants** in the
context of pathogen proteins. Nucleotide substitutions are kept as supporting
evidence for each protein variant where a direct mapping is available, and as a
separate list of genome-level changes.

Future work: replace the Nextclade-derived call-set with an independent caller
that consumes the Stage 3 MAFFT alignment and the curated reference proteome,
so non-Nextclade pathogens can be handled and reading-frame validation is
explicit.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any


def _parse_mutation_token(token: str) -> dict[str, Any] | None:
    """Parse a Nextclade mutation token.

    Handles:
      - substitutions:  GP:A82V
      - deletions:      GP:T30- or GP:T30del
      - insertions:     GP:336:KVTPTSFANNQTSKNHEDLVP

    Returns a dict with gene, position, ref_aa, alt_aa, hgvs_p or None if the
    token cannot be parsed.
    """
    if not token or token.strip() in {"", "nan", "NA"}:
        return None
    token = token.strip()

    # Substitution or deletion: gene:refPOSalt  (alt may be -, del, *, or AA)
    m = re.match(r"^([A-Za-z0-9_]+):([A-Z*]?)(\d+)([A-Z*\-]|del|DEL)?$", token)
    if m:
        gene, ref, pos, alt = m.groups()
        alt = alt or "-"  # deletion if missing
        return {
            "gene": gene,
            "position": int(pos),
            "ref_aa": ref or None,
            "alt_aa": alt,
            "hgvs_p": f"{gene}:{ref or ''}{pos}{alt}",
        }

    # Insertion: gene:POS:insertedSeq
    m = re.match(r"^([A-Za-z0-9_]+):(\d+):([A-Z*\-]+)$", token)
    if m:
        gene, pos, inserted = m.groups()
        return {
            "gene": gene,
            "position": int(pos),
            "ref_aa": "-",
            "alt_aa": inserted,
            "hgvs_p": f"{gene}:{pos}ins{inserted}",
        }

    return {"raw": token}


def _parse_mutation_list(value: str) -> list[dict[str, Any]]:
    """Parse a comma-separated list of Nextclade mutation tokens."""
    results: list[dict[str, Any]] = []
    for token in value.split(","):
        parsed = _parse_mutation_token(token)
        if parsed:
            results.append(parsed)
    return results


def _attach_nt_support(aa_variants: list[dict[str, Any]], nt_substitutions: list[str]) -> None:
    """Attach summary nucleotide evidence to AA variants.

    Nextclade does not provide a per-codon mapping in the TSV, so we record the
    count and (when small) the first few genome-level substitutions. A proper
    codon-level mapping will be added with the independent caller.
    """
    nt_list = [t.strip() for t in nt_substitutions if t.strip()]
    for variant in aa_variants:
        variant["nucleotide_change_supporting_count"] = len(nt_list)
        variant["nucleotide_change_supporting_sample"] = nt_list[:10] if len(nt_list) <= 10 else nt_list[:5]


def _to_tokens(value: Any) -> list[str]:
    """Accept a list of mutation tokens or a comma-separated string."""
    if isinstance(value, str):
        return [t.strip() for t in value.split(",") if t.strip()]
    if isinstance(value, list):
        tokens = []
        for item in value:
            if isinstance(item, str):
                tokens.append(item)
            elif isinstance(item, dict) and "raw" in item:
                tokens.append(str(item["raw"]))
        return tokens
    return []


def extract_amino_acid_variants(lineage_result: dict[str, Any]) -> list[dict[str, Any]]:
    """Return amino-acid-centric variant records from a parsed Nextclade result.

    Each record contains: gene, position, ref_aa, alt_aa, hgvs_p, variant_type,
    and nucleotide_change_supporting.
    """
    mutations = lineage_result.get("mutations", {})

    variants: list[dict[str, Any]] = []
    for token in _to_tokens(mutations.get("aa_substitutions")):
        parsed = _parse_mutation_token(token)
        if parsed and parsed.get("position") is not None:
            parsed["variant_type"] = "substitution"
            variants.append(parsed)

    for token in _to_tokens(mutations.get("aa_deletions")):
        parsed = _parse_mutation_token(token)
        if parsed and parsed.get("position") is not None:
            parsed["variant_type"] = "deletion"
            variants.append(parsed)

    for token in _to_tokens(mutations.get("aa_insertions")):
        parsed = _parse_mutation_token(token)
        if parsed and parsed.get("position") is not None:
            parsed["variant_type"] = "insertion"
            variants.append(parsed)

    nt_substitutions = [
        t.strip() if isinstance(t, str) else str(t)
        for t in _to_tokens(mutations.get("nt_substitutions"))
    ]
    _attach_nt_support(variants, nt_substitutions)

    return variants


def extract_nucleotide_variants(lineage_result: dict[str, Any]) -> list[dict[str, Any]]:
    """Return genome-level nucleotide substitutions as simple records."""
    mutations = lineage_result.get("mutations", {})
    return [{"raw": t.strip() if isinstance(t, str) else str(t)} for t in _to_tokens(mutations.get("nt_substitutions"))]


def run_variant_calling(
    lineage_result: dict[str, Any],
    reference_fasta: Path | None = None,
    reference_proteome: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Produce the Stage 4 variant call-set for a sample.

    Currently this is a thin wrapper around Nextclade output. The returned
    structure is intentionally protein-centric because downstream phenotype
    interpretation (vaccine escape, virulence, etc.) is driven by amino-acid
    changes.
    """
    aa_variants = extract_amino_acid_variants(lineage_result)
    nt_variants = extract_nucleotide_variants(lineage_result)

    return {
        "method": "nextclade",
        "variant_caller": "nextclade",
        "primary_unit": "amino_acid",
        "aa_variants": aa_variants,
        "nt_variants": nt_variants,
        "total_aa_substitutions": len([v for v in aa_variants if v["variant_type"] == "substitution"]),
        "total_aa_deletions": len([v for v in aa_variants if v["variant_type"] == "deletion"]),
        "total_aa_insertions": len([v for v in aa_variants if v["variant_type"] == "insertion"]),
        "note": "Amino-acid variants are primary; nucleotide changes are supporting evidence. "
                "Independent variant caller from alignment is planned.",
    }
