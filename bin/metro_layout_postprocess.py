#!/usr/bin/env python3
"""
Post-process the Nextflow-generated Mermaid DAG into a hierarchical nf-metro
visualisation.

Goals:
- Rename top-level stages to plain, unnumbered labels (Classification,
  <Pathogen> Workflow, Reporting).
- Discover which pathogen workflow subgraph actually ran generically, by the
  `*_WORKFLOW` naming convention (e.g. EBOLA_WORKFLOW, future MPOX_WORKFLOW)
  -- nothing pathogen-specific is hardcoded here. The pathogen itself is
  determined live from species_assignments.tsv, which in turn reflects
  whichever real Nextclade dataset scored best for the samples in this run.
- Rename the PATHOGEN_ROUTER box directly to that pathogen's human-readable
  name (e.g. "Ebola Workflow"), resolved from species_assignments.tsv when
  it matches, else derived generically from the subgraph's own name, and
  flatten the discovered workflow's children directly into it alongside the
  ROUTE_PATHOGEN node (nf-metro renders every subgraph level as a sibling
  station regardless of nesting depth, so there is no benefit to nesting).
- Rename child subworkflows under the pathogen workflow (Bioinformatics
  Analysis, Phenotype Annotation, Epidemiological Data, Evidence Synthesis).
- Reorganise so child modules stack vertically (direction TB) while their
  internal process chains extend horizontally (direction LR) for a compact,
  readable layout.
"""

import argparse
import csv
import re
import sys
from pathlib import Path


PATHOGEN_DISPLAY_NAMES = {
    "orthoebolavirus": "Ebola Workflow",
    "ebola": "Ebola Workflow",
    "ebov": "Ebola Workflow",
    "sars-cov-2": "SARS-CoV-2 Workflow",
    "sarscov2": "SARS-CoV-2 Workflow",
    "flu": "Influenza Workflow",
    "influenza": "Influenza Workflow",
    "mpox": "Mpox Workflow",
    "dengue": "Dengue Workflow",
    "measles": "Measles Workflow",
    "rsv": "RSV Workflow",
    "yellow-fever": "Yellow Fever Workflow",
    "hmpv": "HMPV Workflow",
    "wnv": "WNV Workflow",
}

CHILD_RENAMES = {
    "BIOINFORMATICS_ANALYSIS": "Bioinformatics Analysis",
    "PHENOTYPE_ANNOTATION": "Phenotype Annotation",
    "EPIDEMIOLOGICAL_DATA": "Epidemiological Data",
    "KNOWLEDGE_WAREHOUSE": "Knowledge Warehouse",
    "EVIDENCE_SYNTHESIS": "Evidence Synthesis",
}


def detect_pathogen(assignments_tsv: Path) -> tuple[str, str]:
    """Return (raw_pathogen, display_name) from species_assignments.tsv."""
    if not assignments_tsv.exists():
        return "unknown", "Pathogen Workflow"

    raw_pathogens = set()
    with open(assignments_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pathogen = row.get("pathogen", "").strip().lower()
            if pathogen:
                raw_pathogens.add(pathogen)

    if not raw_pathogens:
        return "unknown", "Pathogen Workflow"

    chosen = sorted(raw_pathogens)[0]
    display = PATHOGEN_DISPLAY_NAMES.get(chosen)
    if not display:
        for key, name in PATHOGEN_DISPLAY_NAMES.items():
            if key in chosen or chosen in key:
                display = name
                break
    if not display:
        display = f"{chosen.title()} Workflow"

    return chosen, display


def derive_workflow_display_name(workflow_label: str, pathogen_raw: str, pathogen_display: str) -> str:
    """
    Resolve a human-readable name for a discovered `<X>_WORKFLOW` subgraph.

    Prefers the display name already resolved from species_assignments.tsv
    when the workflow's own key matches the detected pathogen (handles
    nuanced naming, e.g. species-level detail). Falls back to a generic
    name derived purely from the subgraph's own name -- this is what makes
    discovery of future pathogen workflows (e.g. MPOX_WORKFLOW) work
    without any code changes here.
    """
    workflow_key = re.sub(r'_WORKFLOW$', '', workflow_label, flags=re.IGNORECASE).lower()
    if pathogen_raw and pathogen_raw != "unknown" and (
        workflow_key in pathogen_raw or pathogen_raw in workflow_key
    ):
        return pathogen_display
    return f"{workflow_key.title()} Workflow"


def extract_label(subgraph_line: str) -> str:
    """Return the bracketed label from a `subgraph "... [label]"` line."""
    m = re.search(r'\[(.+?)\]', subgraph_line)
    return m.group(1) if m else ""


def set_label(subgraph_line: str, new_label: str) -> str:
    """Replace the bracketed label in a subgraph declaration."""
    return re.sub(r'\[(.+?)\]', f'[{new_label}]', subgraph_line, count=1)


def parse_block(lines: list[str], start: int) -> tuple[int, list[str]]:
    """
    Extract a balanced subgraph/statement block starting at `start`.
    Returns (end_index_exclusive, block_lines).
    """
    depth = 0
    i = start
    while i < len(lines):
        line = lines[i]
        if re.match(r'\s*subgraph\s', line):
            depth += 1
        if re.match(r'\s*end\s*$', line):
            depth -= 1
            if depth == 0:
                return i + 1, lines[start : i + 1]
        i += 1
    return len(lines), lines[start:]


def find_subgraph_by_label(block: list[str], label: str) -> tuple[int, list[str]] | None:
    """Find a direct child subgraph whose label matches `label`."""
    i = 0
    while i < len(block):
        line = block[i]
        m = re.match(r'(\s*)subgraph\s+"([^"]+)"', line)
        if m:
            current_label = extract_label(line)
            if current_label == label:
                end, child_block = parse_block(block, i)
                return i, child_block
            # The first line may be the parent wrapper itself; step into it
            # rather than skipping the entire block.
            if i == 0:
                i += 1
                continue
            end, _ = parse_block(block, i)
            i = end
        else:
            i += 1
    return None


def flatten_router_block(router_block: list[str], pathogen_raw: str, pathogen_display: str) -> list[str]:
    """
    Rename the outer PATHOGEN_ROUTER subgraph directly to the detected
    pathogen workflow's display name, and flatten the (generically
    discovered) `<X>_WORKFLOW` child subgraph's children directly into it,
    alongside the ROUTE_PATHOGEN node. nf-metro renders every subgraph level
    as a sibling station regardless of nesting depth, so a flat box with
    real children is the layout that actually renders correctly.
    """
    inner_lines = router_block[1:-1]

    # Discover the (first) direct child subgraph matching `*_WORKFLOW` --
    # no hardcoded pathogen list, so this scales to future pathogens.
    workflow_block = None
    wf_start = None
    i = 0
    while i < len(inner_lines):
        line = inner_lines[i]
        m = re.match(r'(\s*)subgraph\s+"([^"]+)"', line)
        if m:
            label = extract_label(line)
            end, sub_block = parse_block(inner_lines, i)
            if re.search(r'_WORKFLOW$', label, re.IGNORECASE):
                workflow_block = sub_block
                wf_start = i
                wf_end = end
                wf_label = label
                break
            i = end
        else:
            i += 1

    router_indent = len(re.match(r'(\s*)', router_block[0]).group(1))
    base_indent = " " * (router_indent + 4)

    if workflow_block is None:
        # No pathogen workflow subgraph found; leave the router as-is
        # (just rename it to the detected pathogen for consistency).
        router_block[0] = set_label(router_block[0], pathogen_display)
        return router_block

    wf_display = derive_workflow_display_name(wf_label, pathogen_raw, pathogen_display)
    router_block[0] = set_label(router_block[0], wf_display)

    # Rename+restructure each child subworkflow inside the discovered
    # pathogen workflow block (direction LR for its internal process chain).
    child_blocks = []
    wf_inner_lines = workflow_block[1:-1]
    j = 0
    while j < len(wf_inner_lines):
        line = wf_inner_lines[j]
        m = re.match(r'(\s*)subgraph\s+"([^"]+)"', line)
        if m:
            child_label = extract_label(line)
            end, child_block = parse_block(wf_inner_lines, j)
            if child_label in CHILD_RENAMES:
                child_block[0] = set_label(child_block[0], CHILD_RENAMES[child_label])
                child_indent = len(m.group(1))
                child_block.insert(1, " " * (child_indent + 4) + "direction LR")
            child_blocks.append((child_label, child_block))
            j = end
        else:
            child_blocks.append((None, [line]))
            j += 1

    # Cosmetic-only reorder for the metro map: move KNOWLEDGE_WAREHOUSE to
    # render after EVIDENCE_SYNTHESIS. This does not reflect (or change) the
    # actual Nextflow call order -- it exists purely so nf-metro's section
    # ordering treats Knowledge Warehouse as non-adjacent to its predecessor
    # sections (Bioinformatics Analysis/Phenotype Annotation/Epidemiological
    # Data), which classifies all of its inbound edges as proper bypass
    # lines instead of misclassifying one as the "main" line. That
    # misclassification previously caused two lines to collide on the same
    # entry port and crash rendering with a CurveInvariantError.
    kw_idx = next((idx for idx, (lbl, _) in enumerate(child_blocks) if lbl == "KNOWLEDGE_WAREHOUSE"), None)
    es_idx = next((idx for idx, (lbl, _) in enumerate(child_blocks) if lbl == "EVIDENCE_SYNTHESIS"), None)
    if kw_idx is not None and es_idx is not None and kw_idx < es_idx:
        kw_entry = child_blocks.pop(kw_idx)
        es_idx = next(idx for idx, (lbl, _) in enumerate(child_blocks) if lbl == "EVIDENCE_SYNTHESIS")
        child_blocks.insert(es_idx + 1, kw_entry)

    new_lines = [router_block[0], base_indent + "direction TB"]
    # Keep ROUTE_PATHOGEN and any other loose lines before the workflow block.
    for line in inner_lines[:wf_start]:
        stripped = line.strip()
        if stripped and stripped != "end":
            new_lines.append(line)
    for _, child_block in child_blocks:
        new_lines.extend(child_block)
    # Keep any loose lines that were after the workflow block.
    for line in inner_lines[wf_end:]:
        stripped = line.strip()
        if stripped and stripped != "end":
            new_lines.append(line)
    new_lines.append(router_block[-1])
    return new_lines


def rewrite_mermaid(mmd_text: str, pathogen_raw: str, pathogen_display: str) -> str:
    lines = mmd_text.splitlines()

    # Locate the main workflow block.
    gi_index = None
    for i, line in enumerate(lines):
        if re.match(r'\s*subgraph\s+"[^"]*GENOMIC_INTELLIGENCE\s+\[GENOMIC_INTELLIGENCE\]"', line):
            gi_index = i
            break

    if gi_index is None:
        return simple_rewrite(mmd_text, pathogen_display)

    gi_start = gi_index
    gi_end, gi_block = parse_block(lines, gi_start)

    # Rename Classification.
    classification_match = find_subgraph_by_label(gi_block, "CLASSIFICATION")
    if classification_match:
        start_idx, class_block = classification_match
        class_block[0] = set_label(class_block[0], "Classification")
        gi_block[start_idx : start_idx + len(class_block)] = class_block

    # Flatten Pathogen Router directly into the detected pathogen workflow
    # box (discovered generically) -- nf-metro doesn't render nested boxes.
    router_match = find_subgraph_by_label(gi_block, "PATHOGEN_ROUTER")
    if router_match:
        router_start_in_gi, router_block = router_match
        new_router_lines = flatten_router_block(router_block, pathogen_raw, pathogen_display)
        gi_block[router_start_in_gi : router_start_in_gi + len(router_block)] = new_router_lines

    # Rename Reporting.
    reporting_match = find_subgraph_by_label(gi_block, "REPORTING")
    if reporting_match:
        rep_idx, rep_block = reporting_match
        rep_block[0] = set_label(rep_block[0], "Reporting")
        gi_block[rep_idx : rep_idx + len(rep_block)] = rep_block

    new_lines = lines[:gi_start] + gi_block + lines[gi_end:]
    return "\n".join(new_lines) + "\n"


def simple_rewrite(mmd_text: str, pathogen_display: str) -> str:
    """Fallback regex-based rewrite when structural parsing fails."""
    mmd_text = re.sub(r'\[CLASSIFICATION\]', '[Classification]', mmd_text)
    mmd_text = re.sub(r'\[REPORTING\]', '[Reporting]', mmd_text)
    mmd_text = re.sub(r'\[PATHOGEN_ROUTER\]', f'[{pathogen_display}]', mmd_text)
    mmd_text = re.sub(r'\[(\w+_WORKFLOW)\]', f'[{pathogen_display}]', mmd_text)
    for old, new in CHILD_RENAMES.items():
        mmd_text = re.sub(rf'\[{old}\]', f'[{new}]', mmd_text)
    return mmd_text


def main():
    parser = argparse.ArgumentParser(
        description="Post-process Nextflow Mermaid DAG for hierarchical nf-metro layout"
    )
    parser.add_argument("--input", required=True, type=Path, help="Input .mmd file")
    parser.add_argument(
        "--assignments", required=True, type=Path, help="species_assignments.tsv"
    )
    parser.add_argument("--output", required=True, type=Path, help="Output .mmd file")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    mmd_text = args.input.read_text()
    pathogen_raw, pathogen_display = detect_pathogen(args.assignments)

    new_mmd = rewrite_mermaid(mmd_text, pathogen_raw, pathogen_display)
    args.output.write_text(new_mmd)
    print(f"Wrote hierarchical metro DAG to {args.output}")


if __name__ == "__main__":
    main()
