#!/usr/bin/env python3
"""DEPRECATED entry point.

Superseded by ``genomic_intelligence/synthesize.py``, which contextualizes
and synthesizes evidence via the LLM (grounded strictly in the
evidence_integration output) instead of rendering a deterministic,
risk-tier-based template. Kept only so old invocations fail loudly with a
pointer to the new entry point, rather than silently producing a
recommendation-flavored brief that no longer matches this layer's contract.
"""

import sys

if __name__ == "__main__":
    sys.exit(
        "generate_intelligence_brief.py is deprecated. "
        "Use: python3 -m intelligence_engine.genomic_intelligence.synthesize"
    )
