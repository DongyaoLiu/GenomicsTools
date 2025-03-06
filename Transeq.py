#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_nucleotide_to_amino_acid(input_fasta, output_fasta):
    # Read the input FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Remove dashes from the sequence
        sequence = str(record.seq).replace("-", "")

        # Create a Seq object from the cleaned sequence
        nucleotide_seq = Seq(sequence)

        # Translate the nucleotide sequence to an amino acid sequence
        amino_acid_seq = nucleotide_seq.translate()

        # Create a new SeqRecord for the amino acid sequence
        amino_acid_record = SeqRecord(amino_acid_seq, id=record.id, description="translated amino acid sequence")

        # Write the amino acid sequence to the output FASTA file
        with open(output_fasta, "a") as output_handle:
            SeqIO.write(amino_acid_record, output_handle, "fasta")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Translate a nucleotide sequence from a FASTA file into an amino acid sequence.")
    parser.add_argument("input_fasta", help="Path to the input FASTA file containing the nucleotide sequence.")
    parser.add_argument("output_fasta", help="Path to the output FASTA file to save the translated amino acid sequence.")

    # Parse arguments
    args = parser.parse_args()

    # Call the translation function
    translate_nucleotide_to_amino_acid(args.input_fasta, args.output_fasta)
    print(f"Translation complete. Amino acid sequence saved to {args.output_fasta}.")

if __name__ == "__main__":
    main()
