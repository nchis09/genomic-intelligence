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
        sample, gene = query_id.rsplit("_", 1)
        return sample, gene
    return query_id, ""


def main():
    args = parse_args()
    out_path = f"{args.prefix}_hmm_annotations.tsv"

    with open(args.tblout) as infile, open(out_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "query_id", "sample", "gene", "vog_id", "vog_accession",
            "tlen", "evalue", "bit_score", "bias", "num_domains",
            "best_domain_evalue", "best_domain_score", "best_domain_bias",
            "exp", "description",
        ])

        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # First 23 hmmscan --tblout fields, then the description.
            parts = line.split(None, 23)
            if len(parts) < 23:
                continue

            query_id = parts[3]
            evalue = float(parts[6])
            if evalue > args.evalue:
                continue

            sample, gene = split_sample_gene(query_id)
            description = parts[23] if len(parts) > 23 else ""

            writer.writerow([
                query_id,
                sample,
                gene,
                parts[0],          # vog_id
                parts[1],          # vog_accession
                parts[2],          # tlen
                parts[6],          # evalue
                parts[7],          # bit_score
                parts[8],          # bias
                parts[9],          # num_domains
                parts[12],         # best_domain_evalue (i-Evalue)
                parts[13],         # best_domain_score
                parts[14],         # best_domain_bias
                parts[15],         # exp
                description,
            ])

    print(f"Wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
