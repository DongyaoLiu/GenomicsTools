#!/usr/bin/env python3
from Bio import AlignIO
import argparse

def convert_fasta_to_phylip(fasta_file, phylip_file, interleaved=False):
    """
    Converts a FASTA file to PHYLIP format.

    Args:
        fasta_file (str): Input FASTA file path.
        phylip_file (str): Output PHYLIP file path.
        interleaved (bool): If True, uses interleaved PHYLIP format.
    """
    try:
        # Determine the PHYLIP format subtype
        phylip_format = "phylip-relaxed" if interleaved else "phylip-sequential"

        # Convert using BioPython
        AlignIO.convert(fasta_file, "fasta", phylip_file, phylip_format)

        print(f"Successfully converted {fasta_file} to PHYLIP format: {phylip_file}")
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

if __name__ == "__main__":
    # Set up command-line arguments
    parser = argparse.ArgumentParser(
        description="Convert a FASTA alignment to PHYLIP format using BioPython."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "-o", "--output",
        default="output.phy",
        help="Output PHYLIP file (default: output.phy)"
    )
    parser.add_argument(
        "--interleaved",
        action="store_true",
        help="Use interleaved PHYLIP format (default: sequential)"
    )

    args = parser.parse_args()

    # Run conversion
    convert_fasta_to_phylip(args.input, args.output, args.interleaved)
