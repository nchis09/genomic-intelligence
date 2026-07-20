#!/usr/bin/env python3
"""
genomic_intelligence/synthesize.py — Genomic Intelligence Assessment entry point.

Acting as an expert genomic intelligence analyst, this stage contextualizes
and synthesizes the molecular and epidemiological evidence produced by
``evidence_integration`` (intelligence_object.json + evidence_package.json)
into a coherent, evidence-grounded narrative -- without performing new
statistics or making public-health recommendations. Every statement in the
narrative is expected to cite the specific finding/metric/source behind it.

Usage:
    python3 -m intelligence_engine.genomic_intelligence.synthesize \\
        --evidence-integration-dir output/evidence_integration/EBOV-UGA-2027-001 \\
        --output-dir output/genomic_intelligence/EBOV-UGA-2027-001

    (--evidence-integration-dir is where intelligence_pipeline.py wrote
    intelligence_object.json + evidence_package/evidence_package.json;
    --output-dir is where this stage writes its own assessment + context.
    They may be the same directory, or different ones. You can also pass
    exact file paths with --intelligence-object / --evidence-package.)
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from intelligence_engine.data_engine.llm_querying.llm_client import LLMClient  # noqa: E402
from intelligence_engine.genomic_intelligence.context_builder import build_context  # noqa: E402

log = logging.getLogger(__name__)


def _load_json(path: Path) -> Optional[dict]:
    if not path or not Path(path).exists():
        return None
    with open(path) as f:
        return json.load(f)


def _deterministic_fallback(document_name: str, context: str) -> str:
    """If no LLM provider is available, fall back to presenting the curated
    evidence context directly (still evidence-only, just without narrative
    synthesis/prose)."""
    return (
        f"{document_name} (deterministic fallback -- no LLM provider available)\n\n"
        "No LLM synthesis could be performed. Below is the curated evidence context that would "
        "have been passed to the model; review it directly.\n\n"
        + context + "\n"
    )


def _strip_markdown_bold(text: str) -> str:
    """Remove ** markdown bold markers from the text."""
    return re.sub(r"\*\*", "", text)


def _wrap_at_width(text: str, width: int = 85) -> str:
    """Wrap lines at *width* characters for readability without horizontal
    scrolling.  Preserves divider lines (===, ---), section headers, and
    leading indentation / bullet markers."""
    lines = text.split("\n")
    out: list[str] = []
    for line in lines:
        stripped = line.rstrip()
        # Keep divider lines as-is
        if re.match(r"^={10,}$", stripped) or re.match(r"^-{10,}$", stripped):
            out.append(stripped)
            continue
        # Keep blank lines
        if not stripped:
            out.append("")
            continue
        # Detect leading whitespace and bullet markers
        lead_match = re.match(r"^(\s*(?:-\s+)?)(.*)$", line)
        prefix = lead_match.group(1) if lead_match else ""
        content = lead_match.group(2) if lead_match else line

        # If the line fits, keep it
        if len(line.rstrip()) <= width:
            out.append(line.rstrip())
            continue

        # Wrap the content portion
        wrapped = textwrap.fill(
            content,
            width=width - len(prefix),
            break_long_words=False,
            break_on_hyphens=False,
        )
        # Re-apply prefix to first line, continuation indent to rest
        wrapped_lines = wrapped.split("\n")
        cont_indent = " " * len(prefix)
        for i, wl in enumerate(wrapped_lines):
            if i == 0:
                out.append(prefix + wl)
            else:
                out.append(cont_indent + wl)
    return "\n".join(out)


def _make_ids(sample_id: str) -> tuple[str, str, str]:
    """Derive Brief/Report IDs and a generation timestamp from the sample ID,
    following the GIB-/GIR-<date>-<sample> convention used in the README
    templates."""
    now = datetime.now(timezone.utc)
    date_tag = now.strftime("%Y%m%d")
    generated_at = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    brief_id = f"GIB-{date_tag}-{sample_id}"
    report_id = f"GIR-{date_tag}-{sample_id}"
    return brief_id, report_id, generated_at


def _copy_figures(figures: list, output_dir: Path) -> list:
    """Copy the figure images referenced in the context into this stage's
    own output directory so the Brief/Report and their images travel
    together, and update each figure's path to the copied location."""
    if not figures:
        return []
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    for fig in figures:
        src = Path(fig["source_path"])
        if not src.exists():
            log.warning(f"Figure source not found, skipping copy: {src}")
            copied.append(fig)
            continue
        dest = figures_dir / fig["filename"]
        try:
            shutil.copy2(src, dest)
        except Exception as e:  # noqa: BLE001
            log.warning(f"Could not copy figure {src} -> {dest}: {e}")
        fig = dict(fig)
        fig["copied_path"] = str(dest)
        copied.append(fig)
    return copied


def synthesize(
    intelligence_object: dict,
    evidence_package: Optional[dict],
    output_dir: Path,
    bio_output: Optional[dict] = None,
    data_query_results: Optional[dict] = None,
    epi_output: Optional[dict] = None,
    analysis_outputs_dir: Optional[Path] = None,
    evidence_integration_dir: Optional[Path] = None,
) -> dict:
    """Build the grounding context, call the LLM, and write the Brief + full
    Report (+ context + figures) to disk."""
    context = build_context(
        intelligence_object,
        evidence_package,
        bio_output=bio_output,
        data_query_results=data_query_results,
        epi_output=epi_output,
        analysis_outputs_dir=analysis_outputs_dir,
        evidence_integration_dir=evidence_integration_dir,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    figures = _copy_figures(context.get("figures") or [], output_dir)

    context_str = json.dumps(context, indent=2, default=str)
    context_path = output_dir / "context_used.json"
    context_path.write_text(context_str)
    log.info(f"Wrote grounding context (for full traceability) to {context_path}")

    sample = intelligence_object.get("sample", {})
    sample_id = sample.get("sample_id", "unknown")
    brief_id, report_id, generated_at = _make_ids(sample_id)

    try:
        client = LLMClient()
        brief = client.synthesize_intelligence_brief(
            context_str, brief_id=brief_id, report_id=report_id,
            sample_id=sample_id, generated_at=generated_at, figures=figures,
        )
    except Exception as e:  # noqa: BLE001 - degrade gracefully, never crash the pipeline for a missing LLM
        log.warning(f"LLM brief synthesis unavailable ({e}); falling back to raw evidence context.")
        brief = _deterministic_fallback("Genomic Intelligence Brief", context_str)

    try:
        client = LLMClient()
        report = client.synthesize_intelligence_report(
            context_str, report_id=report_id, brief_id=brief_id,
            sample_id=sample_id, generated_at=generated_at, figures=figures,
        )
    except Exception as e:  # noqa: BLE001
        log.warning(f"LLM report synthesis unavailable ({e}); falling back to raw evidence context.")
        report = _deterministic_fallback("GENOMIC INTELLIGENCE REPORT", context_str)

    brief = _strip_markdown_bold(brief)
    brief = _wrap_at_width(brief, width=85)

    brief_path = output_dir / "genomic_intelligence_brief.txt"
    brief_path.write_text(brief)
    log.info(f"Wrote genomic intelligence brief to {brief_path}")

    report = _strip_markdown_bold(report)
    report = _wrap_at_width(report, width=85)

    report_path = output_dir / "genomic_intelligence_report.txt"
    report_path.write_text(report)
    log.info(f"Wrote genomic intelligence report to {report_path}")

    return {
        "brief_path": str(brief_path),
        "report_path": str(report_path),
        "context_path": str(context_path),
        "figures_dir": str(output_dir / "figures") if figures else None,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Synthesize a Genomic Intelligence Assessment from evidence_integration outputs."
    )
    parser.add_argument(
        "--evidence-integration-dir",
        default=None,
        help="Directory containing the evidence_integration outputs to read from "
             "(i.e. the --output-dir you passed to intelligence_pipeline.py). "
             "Looks for intelligence_object.json and "
             "evidence_package/evidence_package.json inside it. "
             "Defaults to --output-dir if not given (for backward compatibility).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write this stage's own outputs "
             "(genomic_intelligence_assessment.md, context_used.json). "
             "Can be different from --evidence-integration-dir.",
    )
    parser.add_argument(
        "--intelligence-object",
        default=None,
        help="Explicit path to intelligence_object.json (overrides --evidence-integration-dir).",
    )
    parser.add_argument(
        "--evidence-package",
        default=None,
        help="Explicit path to evidence_package.json (overrides --evidence-integration-dir).",
    )
    parser.add_argument(
        "--bioinformatics-dir",
        default=None,
        help="Directory containing bio_output.json and tree.nwk for the sample "
             "(e.g. output/bioinformatics/EBOV-UGA-2027-001). "
             "Defaults to output/bioinformatics/<sample_id> inferred from the intelligence object.",
    )
    parser.add_argument(
        "--data-query-dir",
        default=None,
        help="Directory containing db_query_results.json and epi_output.json "
             "(e.g. output/data_query/EBOV-UGA-2027-001). "
             "Defaults to output/data_query/<sample_id> inferred from the intelligence object.",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    output_dir = Path(args.output_dir)
    input_dir = Path(args.evidence_integration_dir) if args.evidence_integration_dir else output_dir
    io_path = Path(args.intelligence_object) if args.intelligence_object else input_dir / "intelligence_object.json"
    ep_path = Path(args.evidence_package) if args.evidence_package else input_dir / "evidence_package" / "evidence_package.json"

    intelligence_object = _load_json(io_path)
    if intelligence_object is None:
        sys.exit(
            f"No intelligence_object.json found at {io_path}. "
            "Run intelligence_engine.evidence_integration.pipeline.intelligence_pipeline first."
        )
    evidence_package = _load_json(ep_path)
    if evidence_package is None:
        log.warning(f"No evidence_package.json found at {ep_path}; synthesizing from intelligence_object.json only.")

    sample_id = (intelligence_object.get("sample") or {}).get("sample_id", "unknown")
    bio_dir = Path(args.bioinformatics_dir) if args.bioinformatics_dir else (PROJECT_ROOT / "output" / "bioinformatics" / sample_id)
    dq_dir = Path(args.data_query_dir) if args.data_query_dir else (PROJECT_ROOT / "output" / "data_query" / sample_id)
    analysis_outputs_dir = input_dir / "analysis_outputs"

    bio_output = _load_json(bio_dir / "bio_output.json")
    data_query_results = _load_json(dq_dir / "db_query_results.json")
    epi_output = _load_json(dq_dir / "epi_output.json")

    result = synthesize(
        intelligence_object,
        evidence_package,
        output_dir,
        bio_output=bio_output,
        data_query_results=data_query_results,
        epi_output=epi_output,
        analysis_outputs_dir=analysis_outputs_dir if analysis_outputs_dir.exists() else None,
        evidence_integration_dir=input_dir,
    )
    print(f"Genomic intelligence brief written to: {result['brief_path']}")
    print(f"Genomic intelligence report written to: {result['report_path']}")
    if result.get("figures_dir"):
        print(f"Figures copied to: {result['figures_dir']}")
    print(f"Grounding context (for audit) written to: {result['context_path']}")


if __name__ == "__main__":
    main()
