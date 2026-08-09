#!/usr/bin/env python3
"""Wrapper script: bridges pipeline .txt files and AgentSLR's custom extraction.

Responsibilities:
1. Auto-clone AgentSLR repo if tools/AgentSLR/ doesn't exist.
2. Convert pipeline .txt files → AgentSLR's fulltext_screening.csv format.
3. Set up AgentSLR workspace directory structure.
4. Call AgentSLR main.py --stage data_extraction_custom --domain <domain>.
5. Collect JSONL output and convert to per-PMID JSON files.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

AGENTSLR_REPO_URL = "https://github.com/OxRML/AgentSLR.git"

SPECIES_TO_PATHOGEN = {
    "bdbv": "ebola",
    "ebov": "ebola",
    "sudv": "ebola",
    "tafv": "ebola",
    "restv": "ebola",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run AgentSLR custom evidence extraction on pipeline .txt files."
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing .txt files")
    parser.add_argument("--outdir", required=True, help="Output directory for per-paper JSON files")
    parser.add_argument("--species", required=True, help="Pipeline species key (e.g. bdbv, ebov)")
    parser.add_argument("--domain", required=True, help="Extraction domain (e.g. diagnostic, seroprevalence)")
    parser.add_argument("--agentslr-dir", required=True, help="Path to tools/AgentSLR directory")
    parser.add_argument("--templates-yml", required=True, help="Path to evidence_templates.yml")
    parser.add_argument("--model-name", default="openai/gpt-oss-120b")
    parser.add_argument("--base-url", default="http://localhost:6767/v1")
    parser.add_argument("--api-key", default="6767")
    parser.add_argument("--config-json", default=None)
    parser.add_argument("--metadata-dir", default=None, help="Directory containing PubMed metadata JSON files to merge into output")
    parser.add_argument("--concurrency", type=int, default=1)
    parser.add_argument("--max-completion-tokens", type=int, default=65536)
    parser.add_argument("--reasoning", default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument("--reasoning-effort", choices=["low", "medium", "high"], default="high")
    return parser.parse_args()


def ensure_agentslr_repo(agentslr_dir: Path) -> None:
    """Clone AgentSLR repo if it doesn't exist."""
    if agentslr_dir.exists() and (agentslr_dir / "main.py").exists():
        return
    print(f"[run_agentslr_extraction] AgentSLR not found at {agentslr_dir}. Cloning...")
    agentslr_dir.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["git", "clone", AGENTSLR_REPO_URL, str(agentslr_dir)],
        check=True,
    )
    print(f"[run_agentslr_extraction] AgentSLR cloned to {agentslr_dir}")


def build_fulltext_csv(txt_dir: Path, csv_path: Path) -> pd.DataFrame:
    """Convert .txt files to AgentSLR's fulltext_screening.csv format."""
    rows = []
    for txt_file in sorted(txt_dir.glob("*.txt")):
        article_id = txt_file.stem
        fulltext = txt_file.read_text(encoding="utf-8", errors="replace")
        rows.append({
            "article_id": article_id,
            "fulltext": fulltext,
            "perg_fulltext_result": "INCLUDE",
        })
    df = pd.DataFrame(rows)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    return df


def setup_agentslr_workspace(workspace_dir: Path, pathogen: str, csv_path: Path) -> None:
    """Create the directory structure AgentSLR's Config expects."""
    harvest_root = workspace_dir / "harvests" / pathogen
    client_root = workspace_dir / "client" / "oss" / pathogen
    screening_root = client_root / "screening"
    extractions_root = client_root / "extractions"
    logs_root = client_root / "logs"

    for d in [harvest_root, screening_root, extractions_root, logs_root]:
        d.mkdir(parents=True, exist_ok=True)

    # AgentSLR's Config.fulltext_screening_path points to screening_root / "fulltext_screening.csv"
    target_csv = screening_root / "fulltext_screening.csv"
    shutil.copy2(str(csv_path), str(target_csv))


def run_agentslr(args, agentslr_dir: Path, workspace_dir: Path, pathogen: str) -> Path:
    """Call AgentSLR main.py with the custom extraction stage."""
    cmd = [
        sys.executable,
        str(agentslr_dir / "main.py"),
        "--stage", "data_extraction_custom",
        "--pathogen", pathogen,
        "--data-dir", str(workspace_dir),
        "--domain", args.domain,
        "--evidence-templates-yml", str(Path(args.templates_yml).resolve()),
        "--model-name", args.model_name,
        "--base-url", args.base_url,
        "--api-key", args.api_key,
        "--data-extraction-concurrency", str(args.concurrency),
        "--max-completion-tokens", str(args.max_completion_tokens),
        "--reasoning-effort", args.reasoning_effort,
    ]
    if args.config_json:
        cmd.extend(["--config-json", str(Path(args.config_json).resolve())])
    if not args.reasoning:
        cmd.append("--no-reasoning")

    print(f"[run_agentslr_extraction] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(agentslr_dir))
    if result.returncode != 0:
        # Check if partial output exists before failing
        client_root_check = workspace_dir / "client" / "oss" / pathogen
        output_check = client_root_check / "extractions" / "data_extraction_custom.jsonl"
        if not output_check.exists():
            logs_root_check = client_root_check / "logs"
            candidates = list(logs_root_check.rglob("data_extraction_custom.jsonl"))
            if candidates:
                output_check = candidates[0]
        if output_check.exists():
            print(f"[run_agentslr_extraction] WARNING: AgentSLR exited with code {result.returncode} but partial output exists. Continuing with partial results.")
        else:
            raise RuntimeError(f"AgentSLR extraction failed with exit code {result.returncode}")

    # Find the output JSONL file
    client_root = workspace_dir / "client" / "oss" / pathogen
    extractions_root = client_root / "extractions"
    output_jsonl = extractions_root / "data_extraction_custom.jsonl"

    if not output_jsonl.exists():
        # Search for it in logs directory as fallback
        logs_root = client_root / "logs"
        candidates = list(logs_root.rglob("data_extraction_custom.jsonl"))
        if candidates:
            output_jsonl = candidates[0]
        else:
            raise FileNotFoundError(f"AgentSLR output JSONL not found at {output_jsonl}")

    return output_jsonl


def load_metadata_map(metadata_dir: Path) -> dict:
    """Load all metadata JSON files from a directory, keyed by PMID (filename stem)."""
    metadata_map = {}
    if metadata_dir is None or not metadata_dir.exists():
        return metadata_map
    for json_file in metadata_dir.glob("*.json"):
        try:
            with open(json_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            metadata_map[json_file.stem] = data
        except (json.JSONDecodeError, OSError) as e:
            print(f"[run_agentslr_extraction] Warning: could not read metadata file {json_file}: {e}")
    return metadata_map


METADATA_FIELDS = [
    "title", "authors", "year", "doi", "journal",
    "publication_date", "keywords", "pmcid",
    "cited_by_count", "is_oa",
]


def convert_jsonl_to_per_paper_json(
    jsonl_path: Path,
    outdir: Path,
    species: str,
    domain: str,
    metadata_map: dict | None = None,
) -> list[dict]:
    """Convert AgentSLR's JSONL output to per-PMID JSON files.

    Returns a list of output dicts (for post-extraction checks).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    if metadata_map is None:
        metadata_map = {}

    results = []

    with open(jsonl_path, "r", encoding="utf-8") as f:
        for line in f:
            entry = json.loads(line)
            article_id = entry.get("article_id", "")
            extraction = entry.get("extraction", [])
            not_found = entry.get("not_found", [])
            status = entry.get("status", "unknown")

            output = {
                "pmid": article_id,
                "species": species,
                "domain": domain,
                "extraction": extraction,
                "not_found": not_found,
                "status": status,
            }

            # Merge metadata if available
            meta = metadata_map.get(str(article_id))
            if meta:
                for field in METADATA_FIELDS:
                    if field in meta:
                        output[field] = meta[field]
                # Replace sequential article_id with real PMID from metadata
                if "pmid" in meta:
                    output["pmid"] = meta["pmid"]
            else:
                if metadata_map:
                    print(f"[run_agentslr_extraction] Warning: no metadata found for PMID {article_id}")
                for field in METADATA_FIELDS:
                    output[field] = None

            out_path = outdir / f"{article_id}.json"
            with open(out_path, "w", encoding="utf-8") as out_f:
                json.dump(output, out_f, indent=2, ensure_ascii=False)

            results.append(output)

    print(f"[run_agentslr_extraction] Wrote per-paper JSON files to {outdir}")
    return results


def preflight_llm_check(base_url: str, api_key: str, model_name: str) -> None:
    """Test LLM connection before running extraction. Fail fast with a clear message."""
    from openai import OpenAI
    client = OpenAI(base_url=base_url, api_key=api_key)
    try:
        client.models.list()
    except Exception as e:
        raise RuntimeError(
            f"LLM server at {base_url} is not reachable. "
            f"If using Ollama, ensure 'ollama serve' is running and model '{model_name}' is pulled. "
            f"Error: {e}"
        ) from e


def main():
    args = parse_args()

    agentslr_dir = Path(args.agentslr_dir).resolve()
    input_dir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    templates_path = Path(args.templates_yml).resolve()
    metadata_dir = Path(args.metadata_dir).resolve() if args.metadata_dir else None

    # 1. Ensure AgentSLR repo exists
    ensure_agentslr_repo(agentslr_dir)

    # 2. Pre-flight LLM connection check
    print(f"[run_agentslr_extraction] Checking LLM connection at {args.base_url}...")
    preflight_llm_check(args.base_url, args.api_key, args.model_name)
    print(f"[run_agentslr_extraction] LLM connection OK.")

    # 3. Map species to pathogen
    pathogen = SPECIES_TO_PATHOGEN.get(args.species, "ebola")

    # 4. Build fulltext CSV from .txt files
    workspace_dir = input_dir.parent / f"agentslr_workspace_{args.species}_{args.domain}"
    csv_path = workspace_dir / "fulltext_screening_input.csv"
    df = build_fulltext_csv(input_dir, csv_path)

    if len(df) == 0:
        print("[run_agentslr_extraction] No .txt files found. Writing empty output.")
        outdir.mkdir(parents=True, exist_ok=True)
        return

    print(f"[run_agentslr_extraction] Found {len(df)} .txt files for species={args.species}, domain={args.domain}")

    # 5. Setup AgentSLR workspace
    setup_agentslr_workspace(workspace_dir, pathogen, csv_path)

    # 6. Run AgentSLR extraction
    output_jsonl = run_agentslr(args, agentslr_dir, workspace_dir, pathogen)

    # 7. Load metadata for merging
    metadata_map = load_metadata_map(metadata_dir) if metadata_dir else {}
    if metadata_map:
        print(f"[run_agentslr_extraction] Loaded {len(metadata_map)} metadata records from {metadata_dir}")

    # 8. Convert JSONL to per-paper JSON (with metadata merge)
    results = convert_jsonl_to_per_paper_json(output_jsonl, outdir, args.species, args.domain, metadata_map)

    # 9. Post-extraction check: fail if all papers failed
    failed_count = sum(1 for r in results if r.get("status") == "failed")
    if results and failed_count == len(results):
        raise RuntimeError(
            f"All {len(results)} papers failed extraction. Check LLM server logs. "
            f"This usually means the model is not running or not responding correctly."
        )

    success_count = len(results) - failed_count
    print(f"[run_agentslr_extraction] Extraction complete: {success_count} succeeded, {failed_count} failed out of {len(results)} total.")

    print(f"[run_agentslr_extraction] Done. Output: {outdir}")


if __name__ == "__main__":
    main()
