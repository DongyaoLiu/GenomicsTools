#!/usr/bin/env python3
import os
import re
import argparse
from Bio import SeqIO

def extract_sequences(list_file, fasta_dir, output_file="extracted_sequences.fasta"):
    # Read the list file: {fasta_filename}\t{sequence_header}
    target_sequences = {}
    with open(list_file, 'r') as f:
        for line in f:
            if line.strip():
                fasta_file, seq_header = line.strip().split('\t')
                if fasta_file not in target_sequences:
                    target_sequences[fasta_file] = []
                target_sequences[fasta_file].append(seq_header)

    # Extract sequences from FASTA files
    extracted_seqs = []
    for fasta_file, headers in target_sequences.items():
        species = fasta_file.split(".")[0]
        fasta_path = os.path.join(fasta_dir, fasta_file)
        if not os.path.exists(fasta_path):
            print(f"Warning: {fasta_file} not found in {fasta_dir}")
            continue

        for record in SeqIO.parse(fasta_path, "fasta"):
            for header in headers:
                if re.search(header, record.description):
                    record.id = f"{species}_{header}"
                    extracted_seqs.append(record)
                    break  # Avoid duplicates if multiple patterns match

    # Write extracted sequences to output
    SeqIO.write(extracted_seqs, output_file, "fasta")
    print(f"Extracted {len(extracted_seqs)} sequences to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences from multiple FASTA files based on a list.")
    parser.add_argument("--list", required=True, help="TSV file: fasta_filename <tab> sequence_header")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing FASTA files")
    parser.add_argument("--output", default="extracted_sequences.fasta", help="Output FASTA file")
    args = parser.parse_args()

    extract_sequences(args.list, args.fasta_dir, args.output)
