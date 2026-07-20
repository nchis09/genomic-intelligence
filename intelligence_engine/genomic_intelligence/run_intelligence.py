#!/usr/bin/env python3
"""
Entry point for the Genomic Epidemic Intelligence System.

Subjects a newly sequenced genome to the curated Reference Library (PGIRL) and
assembles a Genomic Epidemic Intelligence Brief.

The engine performs LOOKUPS and deterministic RULE APPLICATION only.
No AI / ML, no biological inference. Every statement in the brief traces to a
curated reference_library entry (or a flagged external source pending curation).

Pipeline stages (see intelligence_engine/README.md):
    1. bioinformatics              new genome -> species / lineage / called mutations
    2. data_engine                 query reference_library / external sources
    3. data_engine.online_querying enrich from NCBI / GISAID / Pathoplexus / WHO / CDC
    4. evidence_integration      harmonize evidence + run cross-evidence statistical analysis
    5. genomic_intelligence      LLM-synthesized, evidence-grounded assessment (no new stats, no recommendations)

NOTE: engine stages are not yet implemented. This is the orchestration scaffold.
"""

import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
REFERENCE_LIBRARY = PROJECT_ROOT / "reference_library"

SUPPORTED_PATHOGENS = ["ebola", "dengue", "influenza", "rvf", "mpox"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Genomic Epidemic Intelligence System - run a genome against the PGIRL."
    )
    parser.add_argument("--input-fasta", required=True,
                        help="Consensus genome (FASTA).")
    parser.add_argument("--metadata", required=True,
                        help="Sample metadata (collection date, location, host, ...).")
    parser.add_argument("--pathogen", choices=SUPPORTED_PATHOGENS + ["auto"],
                        default="auto", help="Pathogen, or auto-detect (default).")
    parser.add_argument("--output-dir", default="output",
                        help="Where to write the Intelligence Brief.")
    return parser.parse_args()


def main():
    args = parse_args()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("GENOMIC EPIDEMIC INTELLIGENCE SYSTEM")
    print("=" * 60)

    # Stage 1 - bioinformatics
    print("\n[1/5] Input processing: species / lineage / mutation calling ...")
    # TODO: from intelligence_engine.bioinformatics import process_genome

    # Stage 2 - data_engine
    print("[2/5] Evidence lookup: querying reference_library ...")
    # TODO: from intelligence_engine.data_engine import query_library

    # Stage 3 - data_engine.online_querying
    print("[3/5] External linking: NCBI / GISAID / Pathoplexus / WHO / CDC ...")
    # TODO: from intelligence_engine.data_engine.online_querying import enrich

    # Stage 4 - evidence_integration
    print("[4/5] Evidence integration: harmonizing evidence + cross-evidence statistics ...")
    # TODO: from intelligence_engine.evidence_integration import assess

    # Stage 5 - genomic_intelligence
    print("[5/5] Genomic Intelligence: synthesizing evidence-grounded assessment ...")
    # TODO: from intelligence_engine.genomic_intelligence import synthesize

    print("\nScaffold only - engine stages not yet implemented.")
    print(f"Reference library: {REFERENCE_LIBRARY}")


if __name__ == "__main__":
    main()
