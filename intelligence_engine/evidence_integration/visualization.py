"""
evidence_integration/visualization.py — Evidence relationship & summary visualizations.

Pure-Python (networkx + matplotlib) visualizations of the evidence package:
an evidence-relationship network (variant <-> phenotype <-> lineage <-> country),
geographic distribution bar charts, and temporal trend line charts. These are
descriptive/statistical visualizations only -- no risk or conclusion framing.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import networkx as nx

from intelligence_engine.evidence_integration.harmonization import EvidenceObject

log = logging.getLogger(__name__)


def build_evidence_network(evidence_objects: list[EvidenceObject]) -> nx.Graph:
    """Build a graph linking variants, phenotype categories, lineages, and
    outbreak countries, for network visualization/analysis (e.g. centrality)."""
    graph = nx.Graph()

    for obj in evidence_objects:
        graph.add_node(obj.key, kind=obj.level)

        if obj.level == "variant" and obj.lineage:
            lineage_key = f"lineage:{obj.lineage.get('lineage_id') or obj.lineage.get('lineage_name')}"
            graph.add_node(lineage_key, kind="lineage")
            graph.add_edge(obj.key, lineage_key, relation="detected_in_lineage")

        for p in obj.phenotype_associations:
            cat = p.get("phenotype_category")
            if not cat:
                continue
            cat_key = f"phenotype:{cat}"
            graph.add_node(cat_key, kind="phenotype_category")
            graph.add_edge(obj.key, cat_key, relation="associated_with", evidence_strength=p.get("evidence_strength"))

        for ob in obj.historical_outbreaks:
            country = ob.get("country")
            if not country:
                continue
            country_key = f"country:{country}"
            graph.add_node(country_key, kind="country")
            graph.add_edge(obj.key, country_key, relation="reported_in")

    return graph


def save_evidence_network(
    evidence_objects: list[EvidenceObject], output_dir: str, filename: str = "evidence_network.png"
) -> Optional[str]:
    """Render the evidence network to a PNG and also export GraphML for
    downstream tools (e.g. Gephi, igraph/ggraph if that pipeline is added later)."""
    graph = build_evidence_network(evidence_objects)
    if graph.number_of_nodes() == 0:
        log.info("Evidence network is empty; skipping visualization.")
        return None

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        color_map = {
            "variant": "#d62728", "lineage": "#1f77b4",
            "phenotype_category": "#2ca02c", "country": "#9467bd",
        }
        node_colors = [color_map.get(graph.nodes[n].get("kind"), "#7f7f7f") for n in graph.nodes]

        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(graph, seed=42, k=0.6)
        nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=500, alpha=0.9)
        nx.draw_networkx_edges(graph, pos, alpha=0.4)
        nx.draw_networkx_labels(graph, pos, font_size=7)
        plt.title("Evidence Relationship Network")
        plt.axis("off")
        png_path = out_dir / filename
        plt.savefig(png_path, dpi=150, bbox_inches="tight")
        plt.close()
    except Exception as e:  # noqa: BLE001 - visualization is supplementary, never block the pipeline
        log.warning(f"Failed to render evidence network PNG: {e}")
        png_path = None

    graphml_path = out_dir / filename.replace(".png", ".graphml")
    try:
        nx.write_graphml(graph, graphml_path)
    except Exception as e:  # noqa: BLE001
        log.warning(f"Failed to write evidence network GraphML: {e}")
        graphml_path = None

    return str(png_path) if png_path else (str(graphml_path) if graphml_path else None)


def save_geographic_distribution_chart(
    country_counts: dict, output_dir: str, filename: str = "geographic_distribution.png", ylabel: str = "count"
) -> Optional[str]:
    """Render a simple bar chart of counts per country."""
    if not country_counts:
        return None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        items = sorted(country_counts.items(), key=lambda kv: kv[1], reverse=True)[:20]
        labels, values = zip(*items)
        plt.figure(figsize=(10, max(4, 0.35 * len(labels))))
        plt.barh(labels, values, color="#4c72b0")
        plt.xlabel(ylabel)
        plt.title("Geographic Distribution")
        plt.gca().invert_yaxis()
        path = out_dir / filename
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        return str(path)
    except Exception as e:  # noqa: BLE001
        log.warning(f"Failed to render geographic distribution chart: {e}")
        return None


def save_temporal_trend_chart(
    year_value_pairs: list[tuple], output_dir: str, filename: str = "temporal_trend.png", ylabel: str = "count"
) -> Optional[str]:
    """Render a simple line chart of a metric over time."""
    if not year_value_pairs:
        return None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        pairs = sorted(year_value_pairs, key=lambda p: p[0])
        years, values = zip(*pairs)
        plt.figure(figsize=(8, 4))
        plt.plot(years, values, marker="o", color="#c44e52")
        plt.xlabel("year")
        plt.ylabel(ylabel)
        plt.title("Temporal Trend")
        plt.tight_layout()
        path = out_dir / filename
        plt.savefig(path, dpi=150)
        plt.close()
        return str(path)
    except Exception as e:  # noqa: BLE001
        log.warning(f"Failed to render temporal trend chart: {e}")
        return None
