#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_nucleotide_to_amino_acid(args):
    # Read the input FASTA file
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        rawSeq = str(record.seq)
        if args.skip5 is not None:
            rawSeq = rawSeq[args.skip5::]
        if args.skip3 is not None:
            rawSeq = rawSeq[::-args.skip3]
        # Remove dashes from the sequence
        sequence = rawSeq.replace("-", "")

        # Create a Seq object from the cleaned sequence
        nucleotide_seq = Seq(sequence)

        # Translate the nucleotide sequence to an amino acid sequence
        amino_acid_seq = nucleotide_seq.translate()

        # Create a new SeqRecord for the amino acid sequence
        amino_acid_record = SeqRecord(amino_acid_seq, id=record.id, description="")

        # Write the amino acid sequence to the output FASTA file
        with open(args.output_fasta, "a") as output_handle:
            SeqIO.write(amino_acid_record, output_handle, "fasta")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Translate a nucleotide sequence from a FASTA file into an amino acid sequence.")
    parser.add_argument("-i", "--input_fasta", type=str, help="Path to the input FASTA file containing the nucleotide sequence.")
    parser.add_argument("-o", "--output_fasta", type=str, help="Path to the output FASTA file to save the translated amino acid sequence.")
    parser.add_argument("-s5", "--skip5", type=int, help="skip the first 5' n nuc")
    parser.add_argument("-s3", "--skip3", type=int, help="skip the last 3' n nuc")

    # Parse arguments
    args = parser.parse_args()

    # Call the translation function
    translate_nucleotide_to_amino_acid(args)
    print(f"Translation complete. Amino acid sequence saved to {args.output_fasta}.")

if __name__ == "__main__":
    main()
