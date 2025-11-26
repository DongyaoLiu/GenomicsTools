#!/usr/bin/env python3

import re
import argparse

def parse_three_prime_utrs(gtf_path):
    """
    Parses a GTF file to extract 3' UTR regions, aggregating by WormBase gene ID (WBGeneXXXXXXXX).
    
    Parameters:
        gtf_path (str): Path to the GTF file.
    
    Returns:
        dict: A dictionary where:
              - Keys are WormBase gene IDs (e.g., 'WBGene00000001')
              - Values are dicts with keys: 'chr', 'strand', 'start', 'end'
    """
    wb_pattern = re.compile(r'(WBGene\d{8})')
    utr_dict = {}

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip comment lines

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # Skip malformed lines

            feature_type = fields[2]
            if feature_type != 'three_prime_utr':
                continue

            chrom = fields[0]
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue  # Skip if coordinates are not integers

            strand = fields[6]
            attributes = fields[8]

            match = wb_pattern.search(attributes)
            if not match:
                continue

            wb_id = match.group(1)

            if wb_id not in utr_dict:
                utr_dict[wb_id] = {
                    'chr': chrom,
                    'strand': strand,
                    'start': start,
                    'end': end
                }
            else:
                # Expand the region to cover all 3' UTRs for this gene
                existing = utr_dict[wb_id]
                # Optionally: you could warn if chr/strand mismatch, but we merge anyway
                existing['start'] = min(existing['start'], start)
                existing['end'] = max(existing['end'], end)

    return utr_dict


def write_utrs_to_bed(utr_dict, bed_path):
    """
    Writes the aggregated 3' UTR dictionary to a BED file (0-based, half-open).

    Parameters:
        utr_dict (dict): Dictionary from parse_three_prime_utrs.
        bed_path (str): Output BED file path.
    """
    with open(bed_path, 'w') as bed_file:
        for wb_id, info in utr_dict.items():
            chrom = info['chr']
            start_1based = info['start']
            end_1based = info['end']
            strand = info['strand']

            chrom_start = start_1based - 1  # Convert to 0-based
            chrom_end = end_1based          # End stays the same (1-based inclusive → BED exclusive)

            bed_file.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{wb_id}\t0\t{strand}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract and merge three_prime_utr regions from a GTF file by WormBase gene ID, and output as BED."
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input GTF file path.'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output BED file path.'
    )

    args = parser.parse_args()

    # Parse GTF
    print(f"Parsing GTF file: {args.input}")
    utr_data = parse_three_prime_utrs(args.input)

    if not utr_data:
        print("Warning: No three_prime_utr features with valid WBGene IDs found.")
    else:
        print(f"Found {len(utr_data)} genes with three_prime_utr regions.")

    # Write BED
    print(f"Writing BED file: {args.output}")
    write_utrs_to_bed(utr_data, args.output)
    print("Done.")


if __name__ == '__main__':
    main()
