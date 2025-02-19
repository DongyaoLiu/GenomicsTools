#!/usr/bin/env python3
import argparse
import gzip
from multiprocessing import Pool, cpu_count
import os

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
    """Parse a file with gene names and return a list of names."""
    gene_names = []
    with open(name_file, 'r') as f:
        for line in f:
            gene_name = line.strip()
            if gene_name in gene_names:
                print(f"Warning: Redundant gene name '{gene_name}' found in the name list. Duplicates will be ignored.")
            else:
                gene_names.append(gene_name)
    return gene_names

def index_gff_by_name(gff_file):
    """Index the GFF file by gene names and return a dictionary of file offsets."""
    index = {}
    with gzip.open(gff_file, 'rt') if gff_file.endswith('.gz') else open(gff_file, 'r') as f:
        while True:
            start_offset = f.tell()
            line = f.readline()
            if not line:
                break
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attributes = fields[8]
            if 'gene_id' in attributes:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
                if gene_id not in index:
                    index[gene_id] = []
                index[gene_id].append((start_offset, f.tell()))
    return index

def index_gff_by_coordinate(gff_file):
    """Index the GFF file by coordinates and return a list of (chrom, start, end, offset) tuples."""
    index = []
    with gzip.open(gff_file, 'rt') if gff_file.endswith('.gz') else open(gff_file, 'r') as f:
        while True:
            offset = f.tell()
            line = f.readline()
            if not line:
                break
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            index.append((chrom, start, end, offset, f.tell()))
    return index

def subset_gff_by_bed(gff_file, intervals, index, output_file):
    """Subset GFF/GTF file based on overlapping coordinates from a BED file."""
    with gzip.open(gff_file, 'rt') if gff_file.endswith('.gz') else open(gff_file, 'r') as gff, \
         open(output_file, 'w') as out:
        for chrom, start, end, start_offset, end_offset in index:
            feature = (chrom, start, end)
            for interval in intervals:
                if overlaps(interval, feature):
                    gff.seek(start_offset)
                    out.write(gff.read(end_offset - start_offset))
                    break

def subset_gff_by_name(gff_file, gene_names, index, output_file):
    """Subset GFF/GTF file based on gene names using the pre-built index."""
    missing_names = [name for name in gene_names if name not in index]
    if missing_names:
        raise ValueError(f"The following gene names were not found in the GFF file: {missing_names}")

    with gzip.open(gff_file, 'rt') if gff_file.endswith('.gz') else open(gff_file, 'r') as gff, \
         open(output_file, 'w') as out:
        for gene_name in gene_names:
            for start_offset, end_offset in index[gene_name]:
                gff.seek(start_offset)
                out.write(gff.read(end_offset - start_offset))

def overlaps(interval, feature):
    """Check if a BED interval overlaps with a GFF feature."""
    chrom1, start1, end1 = interval
    chrom2, start2, end2 = feature[:3]
    if chrom1 != chrom2:
        return False
    return max(start1, start2) < min(end1, end2)

def worker(args):
    """Worker function for multi-processing."""
    mode, gff_file, data, index, output_file = args
    if mode == 'bed':
        subset_gff_by_bed(gff_file, data, index, output_file)
    elif mode == 'name':
        subset_gff_by_name(gff_file, data, index, output_file)

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
        index = index_gff_by_coordinate(args.gff)
    elif args.mode == 'name':
        data = parse_name_list(args.input)
        index = index_gff_by_name(args.gff)

    # Split data into chunks for multi-processing
    chunk_size = len(data) // args.threads
    chunks = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]

    # Create arguments for each worker
    worker_args = [(args.mode, args.gff, chunk, index, f"{args.output}_part{i}.txt") for i, chunk in enumerate(chunks)]

    # Use multi-processing to process chunks in parallel
    with Pool(args.threads) as pool:
        pool.map(worker, worker_args)

    # Combine the output files
    with open(args.output, 'w') as outfile:
        for i in range(len(chunks)):
            part_file = f"{args.output}_part{i}.txt"
            with open(part_file, 'r') as infile:
                outfile.write(infile.read())
            os.remove(part_file)

if __name__ == "__main__":
    main()
