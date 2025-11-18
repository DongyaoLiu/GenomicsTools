#!/usr/bin/env python3
"""
Script to install ALT alleles from VCF into reference genome FASTA file.
Uses efficient FASTA indexing for random access.
"""

import pysam
import sys
import os
from collections import defaultdict
import argparse

def read_vcf(vcf_file):
    """
    Read VCF file and extract variants grouped by chromosome.

    Returns:
        dict: {chrom: [(pos, ref, alt, id)]} sorted by position
    """
    variants = defaultdict(list)

    try:
        vcf = pysam.VariantFile(vcf_file)
        for record in vcf:
            # Skip multi-allelic sites for simplicity
            if len(record.alts) == 1:
                chrom = record.chrom
                pos = record.pos - 1  # Convert to 0-based indexing
                ref = record.ref
                alt = record.alts[0]
                var_id = record.id if record.id else f"{chrom}_{record.pos}"

                variants[chrom].append((pos, ref, alt, var_id))

        # Sort variants by position for each chromosome
        for chrom in variants:
            variants[chrom].sort(key=lambda x: x[0])

        vcf.close()
        return variants

    except Exception as e:
        print(f"Error reading VCF file: {e}")
        sys.exit(1)

def install_alt_alleles(reference_fasta, variants, output_fasta):
    """
    Install ALT alleles into reference genome and write modified FASTA.

    Args:
        reference_fasta (str): Path to reference FASTA file
        variants (dict): Variants from read_vcf()
        output_fasta (str): Path for output FASTA file
    """
    try:
        # Open reference FASTA with indexing
        fasta = pysam.FastaFile(reference_fasta)

        # Get all reference sequences
        references = fasta.references
        reference_lengths = {ref: fasta.get_reference_length(ref) for ref in references}

        # Open output file
        with open(output_fasta, 'w') as out_fh:
            for chrom in references:
                print(f"Processing chromosome: {chrom}")

                # Get reference sequence
                ref_seq = fasta.fetch(chrom)

                if chrom in variants:
                    # Apply variants to this chromosome
                    modified_seq = apply_variants_to_sequence(ref_seq, variants[chrom])
                else:
                    modified_seq = ref_seq

                # Write modified sequence to output FASTA
                out_fh.write(f">{chrom}\n")

                # Write sequence in lines of 80 characters (standard FASTA format)
                for i in range(0, len(modified_seq), 80):
                    out_fh.write(modified_seq[i:i+80] + "\n")

        fasta.close()
        print(f"Modified genome written to: {output_fasta}")

    except Exception as e:
        print(f"Error processing sequences: {e}")
        sys.exit(1)

def apply_variants_to_sequence(sequence, chrom_variants):
    """
    Apply variants to a single chromosome sequence.

    Args:
        sequence (str): Original reference sequence
        chrom_variants (list): List of variants for this chromosome

    Returns:
        str: Modified sequence with ALT alleles installed
    """
    # Convert to list for modification
    seq_list = list(sequence)

    # Track offset changes due to indels
    offset = 0

    for pos, ref, alt, var_id in chrom_variants:
        adjusted_pos = pos + offset

        # Check if reference matches expected
        ref_in_seq = sequence[pos:pos + len(ref)]
        if ref_in_seq.upper() != ref.upper():
            print(f"Warning: Reference mismatch at {var_id}. Expected '{ref}', found '{ref_in_seq}'")
            continue

        # Handle different variant types
        if len(ref) == len(alt):
            # SNP or MNP
            for i in range(len(ref)):
                if adjusted_pos + i < len(seq_list):
                    if i < len(alt):
                        seq_list[adjusted_pos + i] = alt[i]
                    else:
                        # This shouldn't happen for same length variants
                        pass

        elif len(ref) < len(alt):
            # Insertion
            insertion = alt[len(ref):]
            seq_list[adjusted_pos:adjusted_pos + len(ref)] = alt[:len(ref)]

            # Insert remaining bases
            for i, base in enumerate(insertion):
                seq_list.insert(adjusted_pos + len(ref) + i, base)

            offset += len(insertion)

        elif len(ref) > len(alt):
            # Deletion
            deletion_length = len(ref) - len(alt)
            seq_list[adjusted_pos:adjusted_pos + len(alt)] = alt

            # Remove deleted bases
            for i in range(deletion_length):
                if adjusted_pos + len(alt) < len(seq_list):
                    seq_list.pop(adjusted_pos + len(alt))

            offset -= deletion_length

    return ''.join(seq_list)

def main():
    parser = argparse.ArgumentParser(
        description='Install ALT alleles from VCF into reference genome FASTA'
    )
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--output', required=True, help='Output modified FASTA file')

    args = parser.parse_args()

    # Check if input files exist
    if not os.path.exists(args.vcf):
        print(f"Error: VCF file {args.vcf} not found")
        sys.exit(1)

    if not os.path.exists(args.reference):
        print(f"Error: Reference FASTA file {args.reference} not found")
        sys.exit(1)

    # Check if reference FASTA is indexed, create index if not
    if not os.path.exists(args.reference + '.fai'):
        print("Creating FASTA index...")
        pysam.faidx(args.reference)

    print("Reading VCF file...")
    variants = read_vcf(args.vcf)

    print(f"Found variants for chromosomes: {list(variants.keys())}")
    total_variants = sum(len(v) for v in variants.values())
    print(f"Total variants to install: {total_variants}")

    print("Installing ALT alleles into reference...")
    install_alt_alleles(args.reference, variants, args.output)

    print("Done!")

if __name__ == "__main__":
    main()
