#!/usr/bin/env python3

import re
import argparse
from Bio import SeqIO
from collections import defaultdict
from psl_process2 import extract_truncated_sequence, find_polyA

def parse_three_prime_utrs(gtf_path=None, genome=None):
    """
    Parses a GTF file to extract 3' UTR regions, aggregating by WormBase gene ID (WBGeneXXXXXXXX).
    This version is little different from the orignal one. It returns the coordinates with from ATG to the putative 3'utr 
    
    Parameters:
        gtf_path (str): Path to the GTF file.
    
    Returns:
        dict: A dictionary where:
              - Keys are WormBase gene IDs (e.g., 'WBGene00000001')
              - Values are dicts with keys: 'chr', 'strand', 'start', 'end'
    """
    if genome:
        genome_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
        print("for Gene doesn't have 3utr use AAT/AAA to predict")
    else:
        print("Error, for expand mode pls provide the reference genome")
    
    wb_pattern = re.compile(r'(WBGene\d{8})')
    utr_dict = defaultdict(lambda: defaultdict())
    gene_dict = defaultdict(lambda: defaultdict())
    gene_list = []
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip comment lines

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # Skip malformed lines

            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            if strand == "+":
                strand_int = 1
            else:
                strand_int = -1
            attributes = fields[8]
            length = end - start
            
            match = wb_pattern.search(attributes)
            if not match:
                continue
            wb_id = match.group(1)        

    
            if feature_type == "gene":
                # only consider the protein_coding gene.
                if not re.search("protein_coding", fields[8]):
                    continue
                else:
                    gene_list.append(wb_pattern.search(attributes).group(1))
                    gene_dict[wb_id] = {'chrom': chrom,
                                        'strand':strand,
                                        'start': start,
                                        'strand_int': strand_int,
                                        'end': end}
        
            elif feature_type == 'three_prime_utr':
            
                if not utr_dict[wb_id] or utr_dict[wb_id]['length'] < length:
                    utr_dict[wb_id] = {
                        'chrom': chrom,
                        'strand': strand,
                        'start': start,
                        'end': end,
                        'length': length
                    }
                
    for gene, coord in gene_dict.items():
        
        if gene not in utr_dict:
        # use ployA finder here 
            if coord['strand_int']  == 1:
                utr_start = coord['end'] + 1
                utr_end = utr_start + 1000 
            else:
                utr_end = coord['start'] - 1
                utr_start = utr_end - 1000 

            seq = extract_truncated_sequence(genome_dict, coord['chrom'], utr_start, utr_end, coord['strand_int'])
            index = find_polyA(seq) 
            
            if coord['strand_int'] == 1:
                update_gene_dict[gene] = {"chrom": coord["chrom"],
                                  "start": coord['start'],
                                  "end": utr_start + index,
                                  "strand": "+"}
            else:
                update_gene_dict[gene] = {"chrom": coord["chrom"],
                                  "start": utr_end - index,
                                  "end": coord['end'],
                                  "strand": "-"
                                 }


    return update_gene_dict


def write_utrs_to_bed(utr_dict, bed_path):
    """
    Writes the aggregated 3' UTR dictionary to a BED file (0-based, half-open).

    Parameters:
        utr_dict (dict): Dictionary from parse_three_prime_utrs.
        bed_path (str): Output BED file path.
    """
    rename_dict = {"I": "BX284601.5", "II":"BX284602.5", "III":"BX284603.4", "IV":"BX284604.4",
                   "V": "BX284605.5", "X": "BX284606.5"}
    
    
    with open(bed_path, 'w') as bed_file:
        for wb_id, info in utr_dict.items():
            if info['chrom'] == "MtDNA":
                continue
            chrom = rename_dict[info['chrom']]
            start_1based = info['start']
            end_1based = info['end']
            strand = info['strand']

            chrom_start = start_1based - 1  # Convert to 0-based
            chrom_end = end_1based          # End stays the same (1-based inclusive → BED exclusive)

            bed_file.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{wb_id}\t0\t{strand}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract and merge three_prime_utr regions from a GTF file by WormBase gene ID, and output as BED, if multiple transcript, use the longest long 3utr."
    )
    parser.add_argument(
        '-i', '--gtf',
        required=True,
        help='Input GTF file path.'
    )
    parser.add_argument(
        '-g', '--genome',
        default = None,
        required=False,
        help='Input genome fasta file path'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output BED file path.'
    )

    args = parser.parse_args()

    # Parse GTF
    print(f"Parsing GTF file: {args.gtf}")
    utr_data = parse_three_prime_utrs(args.gtf, args.genome)

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
