#!/usr/bin/env python3
"""
Post-process the Nextflow-generated Mermaid DAG into a hierarchical nf-metro
visualisation.

Goals:
- Keep the nf-metro output faithful to the real Nextflow DAG and data
  dependencies; do not force conceptual ordering into the technical DAG.
- Rename top-level stages to plain, unnumbered labels (Classification,
  <Pathogen> Workflow, Reporting).
- Optionally generate a separate conceptual Genomic Intelligence architecture
  diagram showing Input/QC -> Pathogen Identification -> Pathogen-Specific
  Genomic Analysis -> evidence integration -> Genomic Intelligence ->
  Dashboard + Intelligence Brief.
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
  Analysis, Phenotype Annotation, Epidemiological Data, Figures).
- Reorganise so child modules stack vertically (direction TB) while their
  internal process chains extend horizontally (direction LR) for a compact,
  readable layout.
- Strip the pipeline-lifetime shared-DB lifecycle nodes (START_KNOWLEDGE_DB,
  STOP_KNOWLEDGE_DB) entirely -- they are infrastructure plumbing with no
  scientific output, not a real pipeline stage.
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
    "BIOINFORMATICS_AND_EPIDEMIOLOGICAL": "Bioinformatics + Epidemiology",
    "BIOINFORMATICS_ANALYSIS": "Bioinformatics Analysis",
    "LITERATURE_RETRIEVAL": "Literature Retrieval",
    "PHENOTYPE_ANNOTATION": "Phenotype Annotation",
    "EPIDEMIOLOGICAL_DATA": "Epidemiological Data",
    "KNOWLEDGE_WAREHOUSE": "Knowledge Warehouse",
}

# Plain-process nodes that are pipeline-lifetime infrastructure (not a real
# analysis stage) and should be stripped from the rendered diagram entirely.
HIDDEN_NODE_LABELS = {"START_KNOWLEDGE_DB", "STOP_KNOWLEDGE_DB"}


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

    # Cosmetic-only reorder (same rationale as above): move PHENOTYPE_ANNOTATION
    # to render after EPIDEMIOLOGICAL_DATA. PHENOTYPE_ANNOTATION's entry port
    # receives edges from both CLASSIFICATION (species_assignments) and
    # BIOINFORMATICS_ANALYSIS. When PHENOTYPE_ANNOTATION sits immediately after
    # BIOINFORMATICS_ANALYSIS, nf-metro classifies one of these two inbound
    # edges as the "main" line and the other as a bypass line landing on the
    # same entry port, which collide and crash rendering with a
    # CurveInvariantError. Making it non-adjacent to BIOINFORMATICS_ANALYSIS
    # causes nf-metro to route both edges as distinct bypass lines instead.
    pa_idx = next((idx for idx, (lbl, _) in enumerate(child_blocks) if lbl == "PHENOTYPE_ANNOTATION"), None)
    epi_idx = next((idx for idx, (lbl, _) in enumerate(child_blocks) if lbl == "EPIDEMIOLOGICAL_DATA"), None)
    if pa_idx is not None and epi_idx is not None and pa_idx < epi_idx:
        pa_entry = child_blocks.pop(pa_idx)
        epi_idx = next(idx for idx, (lbl, _) in enumerate(child_blocks) if lbl == "EPIDEMIOLOGICAL_DATA")
        child_blocks.insert(epi_idx + 1, pa_entry)

    # Cosmetic-only reorder (same rationale as above): move LITERATURE_RETRIEVAL
    # to render last among the flattened children. Left in its raw call-order
    # slot (2nd of 3, between BIOINFORMATICS_AND_EPIDEMIOLOGICAL and
    # PHENOTYPE_ANNOTATION), its line lands on the same vertical channel as
    # the workflow's own exit edge to the top-level KNOWLEDGE_WAREHOUSE
    # sibling, producing a collinear-overlap CurveInvariantError
    # ("line 'ebola_workflow_knowledge_warehouse' ... and line
    # 'ebola_workflow_literature_retrieval' ... coincide on the V channel").
    # Moving it to the last slot changes its vertical position enough for
    # nf-metro to route both edges as distinct lines instead.
    lit_idx = next((idx for idx, (lbl, _) in enumerate(child_blocks) if lbl == "LITERATURE_RETRIEVAL"), None)
    if lit_idx is not None and lit_idx != len(child_blocks) - 1:
        lit_entry = child_blocks.pop(lit_idx)
        child_blocks.append(lit_entry)

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

    # Flatten Pathogen Identification's subworkflow subgraph
    # (PATHOGEN_IDENTIFICATION_WF, formerly PATHOGEN_CHARACTERIZATION_WF)
    # into Knowledge Warehouse.
    #
    # NOTE: A standalone "Pathogen Characterization"/"Pathogen Identification"
    # section reproducibly triggers the identical nf-metro CurveInvariantError
    # at the *same* coordinates (182.0,610.0) in every configuration tested,
    # including as recently as re-tested after the LITERATURE_RETRIEVAL
    # reorder fix (see CHILD_RENAMES/flatten_router_block above) -- top-level
    # sibling of Knowledge Warehouse, nested inside it, last child of
    # GENOMIC_INTELLIGENCE, this subworkflow's own *native*
    # Nextflow-generated subgraph (not fabricated via text surgery), and
    # even after fully decoupling it from the live DB so its only input is
    # a plain file-based DuckDB-export combine (see
    # subworkflows/local/pathogen_identification/main.nf) -- a
    # structurally different edge than every earlier test. This rules out
    # DAG structure/ordering/edge-count/edge-type as the cause -- it is a
    # defect in nf-metro's own fold/entry-port layout engine (confirmed
    # present in the latest installed version, nf_metro==1.1.0) for this
    # specific section shape. Flattening its single child node into
    # Knowledge Warehouse is the only configuration that renders
    # successfully. Do not attempt to un-flatten this again without first
    # re-verifying against a fresh raw DAG via the fast local render loop
    # (see docs/architecture_overview.mmd for the intended conceptual
    # positioning instead -- a separate, hand-authored diagram not subject
    # to this constraint).
    pc_match = find_subgraph_by_label(gi_block, "PATHOGEN_IDENTIFICATION_WF")
    kw_match = find_subgraph_by_label(gi_block, "KNOWLEDGE_WAREHOUSE")
    if pc_match and kw_match:
        pc_idx, pc_block = pc_match
        pc_inner = [line for line in pc_block[1:-1]]
        del gi_block[pc_idx : pc_idx + len(pc_block)]
        kw_idx, kw_block = find_subgraph_by_label(gi_block, "KNOWLEDGE_WAREHOUSE")
        for j in range(len(kw_block) - 1, -1, -1):
            if re.match(r'^\s*end\s*$', kw_block[j]):
                kw_block[j:j] = pc_inner
                break
        gi_block[kw_idx : kw_idx + len(kw_block)] = kw_block

    new_lines = lines[:gi_start] + gi_block + lines[gi_end:]
    new_text = "\n".join(new_lines) + "\n"
    new_text = new_text.replace('["KNOWLEDGE_WAREHOUSE"]', '["Knowledge Warehouse"]')
    return new_text


def strip_hidden_nodes(mmd_text: str, labels: set[str]) -> str:
    """
    Remove plain-process nodes whose label is in `labels` (e.g.
    START_KNOWLEDGE_DB/STOP_KNOWLEDGE_DB), along with every edge line that
    references them, from the final Mermaid text. These are pipeline-
    lifetime infrastructure steps with no scientific output, not a real
    analysis stage, so they're hidden from the diagram entirely rather than
    given a station.

    Any neighbour left without an edge on one side (e.g. the channel marker
    nodes feeding/fed by these processes) simply renders as an origin/
    terminal marker, the same harmless pattern already produced elsewhere
    in these DAGs for unconsumed channel outputs.
    """
    node_decl_re = re.compile(r'^\s*(\w+)\(\["([^"]+)"\]\)\s*$')
    edge_re = re.compile(r'^\s*(\w+)\s*-->\s*(\w+)\s*$')

    node_ids = set()
    kept_lines = []
    for line in mmd_text.splitlines():
        m = node_decl_re.match(line)
        if m and m.group(2) in labels:
            node_ids.add(m.group(1))
            continue
        kept_lines.append(line)

    if not node_ids:
        return mmd_text

    final_lines = [
        line
        for line in kept_lines
        if not (
            (m := edge_re.match(line))
            and (m.group(1) in node_ids or m.group(2) in node_ids)
        )
    ]
    return "\n".join(final_lines) + "\n"


def simple_rewrite(mmd_text: str, pathogen_display: str) -> str:
    """Fallback regex-based rewrite when structural parsing fails."""
    mmd_text = re.sub(r'\[CLASSIFICATION\]', '[Classification]', mmd_text)
    mmd_text = re.sub(r'\[REPORTING\]', '[Reporting]', mmd_text)
    mmd_text = re.sub(r'\[PATHOGEN_ROUTER\]', f'[{pathogen_display}]', mmd_text)
    mmd_text = re.sub(r'\[(\w+_WORKFLOW)\]', f'[{pathogen_display}]', mmd_text)
    for old, new in CHILD_RENAMES.items():
        mmd_text = re.sub(rf'\[{old}\]', f'[{new}]', mmd_text)
    return mmd_text


def write_conceptual_architecture(output: Path) -> None:
    """
    Write a separate conceptual architecture diagram.

    This diagram is intentionally not passed through nf-metro. It represents
    the framework's conceptual information flow rather than the literal
    Nextflow execution DAG.
    """
    conceptual = """flowchart LR

    A["Input / QC"]
    B["Pathogen Identification"]
    C["Pathogen-Specific<br/>Genomic Analysis"]
    D["Epidemiological<br/>Context"]
    E["Literature &<br/>Functional Evidence"]
    F["Knowledge Warehouse"]
    G["Genomic Intelligence"]
    H["Dashboard"]
    I["Intelligence Brief"]

    A --> B
    B --> C

    C --> D
    C --> E
    C --> F

    D --> F
    E --> F

    F --> G

    G --> H
    G --> I

    classDef stage fill:#ffffff,stroke:#4A6C8C,stroke-width:2px,color:#1F2937;
    classDef intelligence fill:#4A6C8C,stroke:#4A6C8C,stroke-width:2px,color:#ffffff;

    class A,B,C,D,E,F stage;
    class G,H,I intelligence;
"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(conceptual)
    print(f"Wrote conceptual architecture to {output}")



def main():
    parser = argparse.ArgumentParser(
        description="Post-process Nextflow Mermaid DAG for hierarchical nf-metro layout"
    )
    parser.add_argument("--input", required=True, type=Path, help="Input .mmd file")
    parser.add_argument(
        "--assignments", required=True, type=Path, help="species_assignments.tsv"
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="Output technical nf-metro .mmd file"
    )
    parser.add_argument(
        "--conceptual-output",
        type=Path,
        default=None,
        help="Optional path for a separate conceptual Genomic Intelligence architecture .mmd file",
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    mmd_text = args.input.read_text()
    pathogen_raw, pathogen_display = detect_pathogen(args.assignments)

    new_mmd = rewrite_mermaid(mmd_text, pathogen_raw, pathogen_display)
    new_mmd = strip_hidden_nodes(new_mmd, HIDDEN_NODE_LABELS)
    args.output.write_text(new_mmd)
    print(f"Wrote hierarchical metro DAG to {args.output}")

    if args.conceptual_output:
        write_conceptual_architecture(args.conceptual_output)


if __name__ == "__main__":
    main()