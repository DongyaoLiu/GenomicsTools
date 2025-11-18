#!/usr/bin/env python3
from Bio import SeqIO
import argparse


def format_codons(seq, codon_length=3):
    """Split sequence into codon triplets, replace stop codons with gaps (...), and handle existing gaps."""
    seq = seq.replace('-', '.')  # Convert dashes to dots (PAML format)
    codons = []
    for i in range(0, len(seq), codon_length):
        codon = seq[i:i+codon_length].upper().ljust(codon_length, '.')
        # Replace stop codons (TAA, TAG, TGA) with gaps
        if codon in {'TAA', 'TAG', 'TGA'}:
            continue
        codons.append(codon)
    return ' '.join(codons)


def convert_to_paml_phylip(fasta_file, output_file):
    """Convert FASTA to PAML PHYLIP format with codon spacing."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) != 2:
        raise ValueError("Input FASTA must contain exactly 2 sequences for pairwise yn00 analysis.")

    # Ensure alignment length is divisible by 3 (codon alignment)
    if len(records[0].seq) % 3 != 0:
        raise ValueError("Alignment length must be divisible by 3 (codon-aligned).")

    # Write to PHYLIP
    with open(output_file, 'w') as f:
        f.write(f"{len(records)} {len(records[0].seq) - 3}\n")  # Header line
        for record in records:
            name = record.id[:10].ljust(10)  # Truncate/pad to 10 chars
            codon_seq = format_codons(str(record.seq)).upper()
            f.write(f"{name} {codon_seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert pairwise codon-aligned FASTA to PAML PHYLIP format for yn00."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", default="output.phy", help="Output PHYLIP file")
    args = parser.parse_args()

    try:
        convert_to_paml_phylip(args.input, args.output)
        print(f"Converted {args.input} to PAML PHYLIP format: {args.output}")
    except Exception as e:
        print(f"Error: {e}")
        exit(1)
