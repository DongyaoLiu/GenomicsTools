#!/usr/bin/env python
from Bio import SeqIO
import os
import sys

def split_fasta_biopython(input_file, output_dir="single_records"):
    """Split a multi-record FASTA file into individual files using Biopython.

    Args:
        input_file (str): Path to input FASTA file
        output_dir (str): Directory to write individual records (default: "single_records")
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Count records processed
    record_count = 0

    try:
        for record in SeqIO.parse(input_file, "fasta"):
            # Create a safe filename from the record ID
            filename = f"{record.id}.fasta"
            output_path = os.path.join(output_dir, filename)

            # Write the single record to file
            SeqIO.write(record, output_path, "fasta")
            record_count += 1

        print(f"Successfully split {input_file} into {record_count} individual records in {output_dir}/")

    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} input.fasta [output_dir]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "single_records"

    split_fasta_biopython(input_file, output_dir)
