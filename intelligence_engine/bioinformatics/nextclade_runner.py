#!/usr/bin/env python3
"""
Shared Nextclade runner for Ebola clade/lineage assignment.

Can be used both for:
  - database/ebola/protein_variants/assign_lineages.py   (reference DB genomes)
  - intelligence_engine/bioinformatics/taxonomic_classification/clade_assignment.py  (new FASTA samples)

Requires Nextclade CLI to be installed:
    /Users/christianndekezi/anaconda3/bin/conda install -c bioconda nextclade -y
    nextclade dataset get --name ebola --output-dir ~/.nextclade/datasets/ebola

Run with anaconda python.
"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# Species → Nextclade dataset name (update if Nextclade changes naming)
SPECIES_DATASET = {
    "EBOV": "ebola",
    "SUDV": "ebola-sudan",
    "BDBV": "ebola-bundibugyo",
    # RESTV, TAFV, BOMV may not have dedicated Nextclade datasets yet.
}

DEFAULT_DATASET_DIR = Path.home() / ".nextclade" / "datasets"


def check_nextclade() -> Optional[str]:
    """Return path to nextclade executable, or None if not found."""
    path = shutil.which("nextclade")
    if path:
        return path
    # Also check anaconda bin
    conda_path = Path("/Users/christianndekezi/anaconda3/bin/nextclade")
    if conda_path.exists():
        return str(conda_path)
    return None


def ensure_dataset(dataset_name: str, dataset_dir: Path = DEFAULT_DATASET_DIR,
                   nextclade_bin: Optional[str] = None) -> Path:
    """Download Nextclade dataset if not already present."""
    dataset_path = dataset_dir / dataset_name
    if dataset_path.exists():
        return dataset_path

    log.info(f"Downloading Nextclade dataset: {dataset_name}")
    dataset_path.mkdir(parents=True, exist_ok=True)
    nextclade = nextclade_bin or check_nextclade() or "nextclade"
    cmd = [
        nextclade,
        "dataset", "get",
        "--name", dataset_name,
        "--output-dir", str(dataset_path),
    ]
    subprocess.run(cmd, check=True)
    return dataset_path


def fetch_fasta(accession: str, email: str = "pgirl_pipeline@local") -> SeqRecord:
    """Fetch a single genome FASTA record from NCBI by accession."""
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record


def run_nextclade_on_fasta(
    input_fasta: Path,
    output_dir: Path,
    dataset_name: str = "ebola",
    nextclade_bin: Optional[str] = None,
) -> Path:
    """Run Nextclade on a FASTA file and return the path to the output TSV."""
    nextclade = nextclade_bin or check_nextclade()
    if not nextclade:
        raise RuntimeError(
            "Nextclade CLI not found. Install with:\n"
            "  /Users/christianndekezi/anaconda3/bin/conda install -c bioconda nextclade -y"
        )

    dataset_path = ensure_dataset(dataset_name, nextclade_bin=nextclade)
    output_tsv = output_dir / "nextclade.tsv"

    cmd = [
        nextclade,
        "run",
        "--input-dataset", str(dataset_path),
        "--output-all", str(output_dir),
        "--output-tsv", str(output_tsv),
        str(input_fasta),
    ]
    log.info(f"Running Nextclade on {input_fasta}")
    subprocess.run(cmd, check=True)
    return output_tsv


def parse_nextclade_tsv(tsv_path: Path) -> dict[str, dict]:
    """Parse Nextclade output TSV and return {seq_name: {clade, qc_status, ...}}."""
    results: dict[str, dict] = {}
    with open(tsv_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split("\t")
            row = dict(zip(header, cols))
            seq_name = row.get("seqName", row.get("seq_name"))
            if not seq_name:
                continue
            results[seq_name] = {
                "clade": row.get("clade", ""),
                "qc_status": row.get("qc.overallStatus", row.get("qc_overallStatus", "")),
                "total_substitutions": row.get("totalSubstitutions", ""),
                "total_deletions": row.get("totalDeletions", ""),
                "total_insertions": row.get("totalInsertions", ""),
                "total_missing": row.get("totalMissing", ""),
                "alignment_score": row.get("alignmentScore", ""),
                "nearest_node_id": row.get("nearestNodeId", ""),
            }
    return results


def assign_lineages_from_accessions(
    accessions: list[str],
    species_id: str,
    email: str = "pgirl_pipeline@local",
) -> dict[str, dict]:
    """Fetch FASTA for accessions, run Nextclade, and return clade assignments."""
    dataset_name = SPECIES_DATASET.get(species_id)
    if not dataset_name:
        raise ValueError(f"No Nextclade dataset configured for species {species_id}")

    with tempfile.TemporaryDirectory(prefix="pgirl_nextclade_") as tmpdir:
        tmp_path = Path(tmpdir)
        fasta_path = tmp_path / "input.fasta"

        log.info(f"Fetching {len(accessions)} FASTA records for {species_id}")
        records = []
        for acc in accessions:
            try:
                rec = fetch_fasta(acc, email=email)
                rec.id = acc
                rec.description = ""
                records.append(rec)
            except Exception as e:
                log.warning(f"Could not fetch {acc}: {e}")

        if not records:
            return {}

        SeqIO.write(records, fasta_path, "fasta")
        output_tsv = run_nextclade_on_fasta(fasta_path, tmp_path, dataset_name=dataset_name)
        return parse_nextclade_tsv(output_tsv)


def assign_lineages_from_fasta(
    fasta_path: Path,
    species_id: str,
) -> dict[str, dict]:
    """Run Nextclade on an existing FASTA file and return clade assignments."""
    dataset_name = SPECIES_DATASET.get(species_id)
    if not dataset_name:
        raise ValueError(f"No Nextclade dataset configured for species {species_id}")

    with tempfile.TemporaryDirectory(prefix="pgirl_nextclade_") as tmpdir:
        output_dir = Path(tmpdir)
        output_tsv = run_nextclade_on_fasta(fasta_path, output_dir, dataset_name=dataset_name)
        return parse_nextclade_tsv(output_tsv)


if __name__ == "__main__":
    # Simple self-test with a single known accession
    import sys
    if len(sys.argv) < 2:
        print("Usage: nextclade_runner.py <accession> [species_id]")
        sys.exit(1)

    acc = sys.argv[1]
    sp = sys.argv[2] if len(sys.argv) > 2 else "EBOV"
    print(assign_lineages_from_accessions([acc], sp))
