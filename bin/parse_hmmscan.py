#!/usr/bin/env python3
"""
parse_hmmscan.py

Parse HMMER3 hmmscan --tblout output and write a tidy TSV.

The FASTA headers produced by extract_query_proteins.py are of the form
{sample_name}_{gene} (e.g. sample01_NP). This script splits the query_id
on the last underscore to recover sample and gene.
"""

import argparse
import csv
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tblout", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    return parser.parse_args()


def split_sample_gene(query_id):
    if "_" in query_id:
        sample, gene = query_id.split("_", 1)
        return sample, gene
    return query_id, ""


def main():
    args = parse_args()
    out_path = f"{args.prefix}_hmm_annotations.tsv"

    with open(args.tblout) as infile, open(out_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "query_id", "sample", "gene", "hmm_id", "hmm_accession",
            "tlen", "evalue", "bit_score", "bias", "num_domains",
            "best_domain_evalue", "best_domain_score", "best_domain_bias",
            "exp", "description",
        ])

        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # hmmscan --tblout: 18 fixed fields then the description.
            parts = line.split(None, 18)
            if len(parts) < 18:
                continue

            query_id = parts[2]
            evalue = float(parts[4])
            if evalue > args.evalue:
                continue

            sample, gene = split_sample_gene(query_id)
            description = parts[18] if len(parts) > 18 else ""

            writer.writerow([
                query_id,
                sample,
                gene,
                parts[0],          # hmm_id
                parts[1],          # hmm_accession
                "",                # tlen not reported by hmmscan
                evalue,
                parts[5],          # bit_score
                parts[6],          # bias
                parts[15],         # num reported domains
                parts[7],          # best_domain_evalue
                parts[8],          # best_domain_score
                parts[9],          # best_domain_bias
                parts[10],         # exp
                description,
            ])

    print(f"Wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
