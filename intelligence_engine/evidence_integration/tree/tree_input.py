"""Tree input handling for the Genomic Intelligence Engine.

Supports a Newick tree file as the fourth engine input. When a tree is not
available or cannot be parsed, the loader falls back to genome metadata so that
decision-oriented analyzers still produce an answer (with lower confidence).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import pandas as pd

log = logging.getLogger(__name__)


@dataclass
class TipAnnotation:
    """Annotation for a single tree tip or metadata-derived surrogate."""

    name: str
    is_sample: bool = False
    accession: str = ""
    country: str = ""
    date: str = ""
    host: str = ""
    lineage_id: str = ""
    branch_length: float = 0.0
    distance_from_root: float = 0.0


@dataclass
class TreeInput:
    """Parsed tree plus metadata-derived annotations.

    If ``tree`` is None, the object acts as a metadata-only placeholder and
    ``has_tree`` returns False.
    """

    tips: list[TipAnnotation] = field(default_factory=list)
    sample_tip: Optional[TipAnnotation] = None
    tree: Any = None
    time_scaled: dict = field(default_factory=dict)
    file_path: Optional[Path] = None
    _tip_by_name: dict[str, TipAnnotation] = field(default_factory=dict, repr=False)

    @property
    def has_tree(self) -> bool:
        return self.tree is not None and bool(self.tips)

    def get_tip(self, name: str) -> Optional[TipAnnotation]:
        return self._tip_by_name.get(name)

    def get_nearest_tips(self, n: int = 3) -> list[TipAnnotation]:
        """Return the n nearest non-sample tips sorted by genetic distance."""
        if not self.sample_tip:
            return []
        candidates = [t for t in self.tips if not t.is_sample]
        candidates.sort(key=lambda t: t.distance_from_root)
        return candidates[:n]

    def get_dated_neighbors(self, window_days: int = 180) -> list[TipAnnotation]:
        """Return nearest tips collected within a date window of the sample."""
        if not self.sample_tip or not self.sample_tip.date:
            return []
        try:
            sample_date = pd.to_datetime(self.sample_tip.date)
        except (ValueError, TypeError):
            return []
        neighbors = []
        for tip in self.tips:
            if tip.is_sample or not tip.date:
                continue
            try:
                tip_date = pd.to_datetime(tip.date)
            except (ValueError, TypeError):
                continue
            if abs((tip_date - sample_date).days) <= window_days:
                neighbors.append(tip)
        return neighbors


def _safe_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    return str(value).strip()


def _parse_date(value: Any) -> str:
    s = _safe_str(value)
    if not s:
        return ""
    # Keep YYYY-MM-DD or YYYY-MM or YYYY
    m = re.match(r"(\d{4})(?:-\d{2}(?:-\d{2})?)?", s)
    return m.group(0) if m else s


def _match_tip_to_metadata(
    tip_name: str, metadata_df: pd.DataFrame, sample_id: str, sample_accession: str
) -> dict[str, Any]:
    """Match a tree tip name to a metadata row.

    Matching order:
      1. exact genome_accession
      2. exact genome_name
      3. contains sample_id or sample_accession
      4. tip name is a known alias in lineages.csv (not implemented here)
    """
    if metadata_df.empty:
        return {}
    name = tip_name.strip()
    lower = name.lower()

    for col in ["genome_accession", "genome_name"]:
        if col in metadata_df.columns:
            mask = metadata_df[col].astype(str).str.strip().str.lower() == lower
            rows = metadata_df[mask]
            if not rows.empty:
                return rows.iloc[0].to_dict()

    # Fuzzy: tip name contains the sample accession/id
    for col in ["genome_accession", "genome_name"]:
        if col in metadata_df.columns:
            mask = (
                metadata_df[col]
                .astype(str)
                .str.lower()
                .str.contains(re.escape(lower), na=False)
            )
            rows = metadata_df[mask]
            if not rows.empty:
                return rows.iloc[0].to_dict()

    # If the tip name is the sample identifier itself, return an empty dict
    # so the caller can annotate it from stage9 metadata.
    if name == sample_id or name == sample_accession:
        return {"_is_sample": True}

    return {}


def _annotate_tips_from_tree(
    tree: Any,
    metadata_df: pd.DataFrame,
    sample_id: str,
    sample_accession: str,
    sample_metadata: dict[str, Any],
) -> list[TipAnnotation]:
    """Walk tree tips and annotate each with metadata."""
    tips: list[TipAnnotation] = []

    try:
        from Bio.Phylo.BaseTree import Clade

        terminals = list(tree.get_terminals())
    except Exception as exc:
        log.warning("Could not enumerate tree terminals: %s", exc)
        return tips

    # Biopython does not populate clade.parent; assign it explicitly.
    def _assign_parents(clade, parent=None):
        clade.parent = parent
        for child in getattr(clade, "clades", []):
            _assign_parents(child, parent=clade)

    _assign_parents(tree.root)

    # Compute distances from root if branch lengths exist
    def dist_to_root(clade):
        total = 0.0
        while clade:
            total += clade.branch_length or 0.0
            if clade == tree.root:
                break
            clade = getattr(clade, "parent", None)
            if clade is None:
                break
        return total

    for terminal in terminals:
        name = terminal.name or ""
        branch_length = terminal.branch_length or 0.0
        distance = dist_to_root(terminal)
        row = _match_tip_to_metadata(name, metadata_df, sample_id, sample_accession)

        is_sample = bool(row.pop("_is_sample", False)) or _is_sample_name(
            name, sample_id, sample_accession
        )

        if is_sample and sample_metadata:
            country = _safe_str(sample_metadata.get("country"))
            date = _parse_date(sample_metadata.get("collection_date"))
            host = _safe_str(sample_metadata.get("host"))
            lineage_id = _safe_str(sample_metadata.get("lineage"))
        else:
            country = _safe_str(row.get("collection_country"))
            date = _parse_date(row.get("collection_date"))
            host = _safe_str(row.get("host"))
            lineage_id = _safe_str(row.get("lineage_id"))

        accession = _safe_str(row.get("genome_accession", name))

        tips.append(
            TipAnnotation(
                name=name,
                is_sample=is_sample,
                accession=accession,
                country=country,
                date=date,
                host=host,
                lineage_id=lineage_id,
                branch_length=branch_length,
                distance_from_root=distance,
            )
        )

    return tips


def _is_sample_name(name: str, sample_id: str, sample_accession: str) -> bool:
    if not name:
        return False
    lower = name.strip().lower()
    return lower == (sample_id or "").strip().lower() or lower == (
        sample_accession or ""
    ).strip().lower()


def _build_metadata_only_tree(
    metadata_df: pd.DataFrame,
    sample_id: str,
    sample_accession: str,
    sample_metadata: dict[str, Any],
    lineage_id: Optional[str] = None,
) -> list[TipAnnotation]:
    """When no tree is provided, use historical metadata tips as surrogates."""
    tips: list[TipAnnotation] = []
    if metadata_df.empty:
        return tips

    mask = pd.Series(True, index=metadata_df.index)
    if lineage_id and "lineage_id" in metadata_df.columns:
        mask = metadata_df["lineage_id"].astype(str).str.strip().str.lower() == lineage_id.lower()

    # Limit to a reasonable number of nearest contextual tips (same lineage or all)
    rows = metadata_df[mask].head(50)
    if rows.empty:
        rows = metadata_df.head(50)

    for _, row in rows.iterrows():
        name = _safe_str(row.get("genome_accession")) or _safe_str(row.get("genome_name"))
        if not name:
            continue
        is_sample = _is_sample_name(name, sample_id, sample_accession)
        tips.append(
            TipAnnotation(
                name=name,
                is_sample=is_sample,
                accession=_safe_str(row.get("genome_accession", name)),
                country=_safe_str(row.get("collection_country")),
                date=_parse_date(row.get("collection_date")),
                host=_safe_str(row.get("host")),
                lineage_id=_safe_str(row.get("lineage_id")),
                branch_length=0.0,
                distance_from_root=0.0,
            )
        )

    # Ensure sample tip exists even if it is not in metadata
    if not any(t.is_sample for t in tips) and sample_id:
        tips.insert(
            0,
            TipAnnotation(
                name=sample_id,
                is_sample=True,
                accession=sample_accession,
                country=_safe_str(sample_metadata.get("country")),
                date=_parse_date(sample_metadata.get("collection_date")),
                host=_safe_str(sample_metadata.get("host")),
                lineage_id=lineage_id or _safe_str(sample_metadata.get("lineage")),
            ),
        )

    return tips


def load_tree_input(
    tree_path: Optional[str],
    metadata_df: pd.DataFrame,
    sample_id: str,
    sample_accession: str,
    sample_metadata: dict[str, Any],
    time_scaled: Optional[dict] = None,
    lineage_id: Optional[str] = None,
) -> TreeInput:
    """Load a tree and annotate tips, or fall back to metadata-only context."""
    time_scaled = time_scaled or {}
    file_path = Path(tree_path) if tree_path else None

    tree = None
    tips: list[TipAnnotation] = []

    if file_path and file_path.exists():
        try:
            from Bio import Phylo

            trees = list(Phylo.parse(str(file_path), "newick"))
            if trees:
                tree = trees[0]
                if not tree.rooted:
                    tree.root_at_midpoint()
                tips = _annotate_tips_from_tree(
                    tree,
                    metadata_df,
                    sample_id,
                    sample_accession,
                    sample_metadata,
                )
        except Exception as exc:
            log.warning("Failed to parse tree %s: %s", file_path, exc)
            tree = None

    if not tips:
        tips = _build_metadata_only_tree(
            metadata_df,
            sample_id,
            sample_accession,
            sample_metadata,
            lineage_id=lineage_id,
        )

    sample_tip = next((t for t in tips if t.is_sample), None)
    if sample_tip is None and tips:
        # If no tip matched the sample, create a virtual sample tip
        sample_tip = TipAnnotation(
            name=sample_id,
            is_sample=True,
            accession=sample_accession,
            country=_safe_str(sample_metadata.get("country")),
            date=_parse_date(sample_metadata.get("collection_date")),
            host=_safe_str(sample_metadata.get("host")),
            lineage_id=lineage_id or _safe_str(sample_metadata.get("lineage")),
        )
        tips.append(sample_tip)

    tip_by_name = {t.name: t for t in tips}
    return TreeInput(
        tips=tips,
        sample_tip=sample_tip,
        tree=tree,
        time_scaled=time_scaled,
        file_path=file_path,
        _tip_by_name=tip_by_name,
    )
