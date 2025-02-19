#!/usr/bin/env python3
import argparse
from Bio.SeqFeature import FeatureLocation
import gzip
from multiprocessing import Pool, cpu_count

def parse_bed(bed_file):
    """Parse BED file and return a list of intervals."""
    intervals = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            intervals.append((chrom, start, end))
    return intervals

def parse_name_list(name_file):
    """Parse a file with gene names and return a set of names."""
    gene_names = set()
    with open(name_file, 'r') as f:
        for line in f:
            gene_names.add(line.strip())
    return gene_names

def parse_gff(gff_file):
    """Parse GFF/GTF file and yield lines with their coordinates and attributes."""
    with gzip.open(gff_file, 'rt') if gff_file.endswith('.gz') else open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            attributes = fields[8]
            yield chrom, start, end, attributes, line

def overlaps(interval, feature):
    """Check if a BED interval overlaps with a GFF feature."""
    chrom1, start1, end1 = interval
    chrom2, start2, end2 = feature[:3]
    if chrom1 != chrom2:
        return False
    return max(start1, start2) < min(end1, end2)

def subset_gff_by_bed(gff_file, intervals, output_file):
    """Subset GFF/GTF file based on overlapping coordinates from a BED file."""
    with open(output_file, 'w') as out:
        for chrom, start, end, attributes, line in parse_gff(gff_file):
            feature = (chrom, start, end)
            for interval in intervals:
                if overlaps(interval, feature):
                    out.write(line)
                    break

def subset_gff_by_name(gff_file, gene_names, output_file):
    """Subset GFF/GTF file based on gene names."""
    with open(output_file, 'w') as out:
        for chrom, start, end, attributes, line in parse_gff(gff_file):
            for gene_name in gene_names:
                if f'gene_name "{gene_name}"' in attributes or f'gene_id "{gene_name}"' in attributes:
                    out.write(line)
                    break

def worker(args):
    """Worker function for multi-processing."""
    mode, gff_file, data, output_file = args
    if mode == 'bed':
        subset_gff_by_bed(gff_file, data, output_file)
    elif mode == 'name':
        subset_gff_by_name(gff_file, data, output_file)

def main():
    parser = argparse.ArgumentParser(description="Subset GFF/GTF file based on BED coordinates or gene names.")
    parser.add_argument('-m', '--mode', required=True, choices=['bed', 'name'], help="Mode: 'bed' for coordinate overlap, 'name' for gene names.")
    parser.add_argument('-i', '--input', required=True, help="Input BED file (for 'bed' mode) or gene name list (for 'name' mode).")
    parser.add_argument('-g', '--gff', required=True, help="Input GFF/GTF file.")
    parser.add_argument('-o', '--output', required=True, help="Output file to save the subset GFF/GTF.")
    parser.add_argument('-t', '--threads', type=int, default=cpu_count(), help="Number of threads to use.")
    args = parser.parse_args()

    # Parse input data based on mode
    if args.mode == 'bed':
        data = parse_bed(args.input)
    elif args.mode == 'name':
        data = parse_name_list(args.input)

    # Split data into chunks for multi-processing
    chunk_size = len(data) // args.threads
    chunks = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]

    # Create arguments for each worker
    worker_args = [(args.mode, args.gff, chunk, f"{args.output}_part{i}.txt") for i, chunk in enumerate(chunks)]

    # Use multi-processing to process chunks in parallel
    with Pool(args.threads) as pool:
        pool.map(worker, worker_args)

    # Combine the output files
    with open(args.output, 'w') as outfile:
        for i in range(len(chunks)):
            part_file = f"{args.output}_part{i}.txt"
            with open(part_file, 'r') as infile:
                outfile.write(infile.read())
            import os
            os.remove(part_file)

if __name__ == "__main__":
    main()
