#!/usr/bin/env python3
import sys
import re
from collections import defaultdict

def process_bed_simple(bed_file):
    """Simple version without pandas."""
    # Read and parse file
    entries = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom, start, end, family  = parts[0], int(parts[1]), int(parts[2]), parts[3]
                    length, strand = int(parts[4]), parts[5]
                    mean = (start + end)/2
                    entries.append((family, chrom, start, end, mean, length, strand))

    # Sort by family, chrom, start
    entries.sort(key=lambda x: (x[0], x[1], x[2]))

    # Process and output
    prev_family = None
    prev_chrom = None
    prev_mean=None

    for family, chrom, start, end, mean, length, strand in entries:
        if family == prev_family and chrom == prev_chrom:
            # Calculate distance
            distance = mean - prev_mean
            location_mean = (prev_mean +  mean)/2
            ID = re.match("R=([0-9]+)", family).group(1)
            print(f"{ID}\t{chrom}\t{distance}\t{location_mean}")

        prev_family = family
        prev_chrom = chrom
        prev_mean = mean

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <bed_file>")
        sys.exit(1)

    process_bed_simple(sys.argv[1])
