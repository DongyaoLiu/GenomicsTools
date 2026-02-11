#!/usr/bin/env python3

import psutil
import re
import os
import bisect
import argparse
import pandas as pd
import glob
from typing import List, Dict, Tuple, Set, Any, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from typing import Dict, Optional
from intervaltree import IntervalTree
from psl_process2 import parse_psl, overlaps_3prime_region, find_polyA, build_genomic_interval_tree, extract_truncated_sequence, parse_celegans_utr_bed, target_genes_3utr

#----------------------------------
# take gff/gtf as in put build the intervals by chromosomes.
#----------------------------------

def parse_annotation(psl_file, interval_level = ["mRNA", "transcript"]):
    output_dict = defaultdict(list)
    with open(psl_file) as f:
        for line in f:
            if line.startswith(('track', '#')) or not line.strip():
                    continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            if cols[2] in interval_level:
                for attr in cols[8].split(';'):
                    if attr.startswith('ID='):
                        name = attr[3:]
                        break
                    if re.search(r"transcript_id \"(\S+)\"", attr):
                        #This is for gtf format extraction
                        name = re.search(r"transcript_id \"(\S+)\"", attr).group(1)
                        break
                
                if cols[6] == "+":
                    strand = 1
                else:
                    strand = -1
                output_dict[cols[0]].append({"gff_start" : int(cols[3]),
                                             "gff_end" : int(cols[4]),
                                             "gene_name" : name, 
                                             "strand" : strand})
    return output_dict


def overlaps_gene_region(gene, index, hit_info, hit_chr, min_overlap = 20):
    # Check if chromosome exists in the trees
    if hit_chr not in index:
        return []

    chrom_tree = index[hit_chr]
    overlapping_intervals = chrom_tree.overlap(hit_info["tStart"], hit_info["tEnd"])
    result = []
    
    length = []
    for interval in overlapping_intervals:
        # Unpack stored data: (feature_index, feature_strand)
        hit_idx, hit_strand = interval.data

        # Calculate overlap
        overlap_start = max(hit_info["tStart"], interval.begin)
        overlap_end = min(hit_info["tEnd"], interval.end)
        overlap_length = overlap_end - overlap_start

        # Skip if overlap is below threshold
        if overlap_length < min_overlap:
            continue

        # Check strand matching
        strand_match = (hit_strand == hit_info["strand"])


        # Skip if strand doesn't match when required
        if not strand_match:
            continue

        # Get the full feature data
        overlapped_hits = gene[hit_chr][hit_idx]

        # Add to results
        result.append(overlapped_hits)
        length.append(overlap_length)
    
    if len(length) == 1:
        return result[0]
    elif result: 
        max_index = max(enumerate(length), key=lambda x: x[1])[0]
        return result[max_index]
    else:
        return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psl", required=True, help="Custom PSL: [celegans_gene_id] + 21 standard cols")
    parser.add_argument("--psl3utr", required=True, help="Custom PSL: [celegans_gene_id] + 21 standard cols")
    parser.add_argument("--celegans-utr-bed", required=True, help="C. elegans 3' UTR BED (col4=geneID, col6=strand)")
    parser.add_argument("--target-gff", required=True, type=str)
    parser.add_argument("--target-genome", required=True)
    #parser.add_argument("--orthofinder", required=True)
    parser.add_argument("--target-species", required=True, help="Target species name")
    #parser.add_argument("--rename-file", required=False)
    parser.add_argument("--output-fasta", default="target_3utrs.fa")
    #parser.add_argument("--stats", default="target_stats.tsv")
    args = parser.parse_args()


    # Load data
    # header_dict = header_mapper(args.rename_file, args.target_species)
    celegans_gene_strand = parse_celegans_utr_bed(args.celegans_utr_bed)
    target_genome = SeqIO.to_dict(SeqIO.parse(args.target_genome, "fasta"))
    target_genes = target_genes_3utr(args.target_gff)

    #load the target species gff file or gtf file
    target_annotation = parse_annotation(args.target_gff)


    #build the interval tree for target_annotation
    tree_index = build_genomic_interval_tree(target_annotation)
    

    # psl import and build the tree
    psl_hits = parse_psl(args.psl, celegans_gene_strand, longest = True)
    psl_hits_3utr = parse_psl(args.psl3utr, celegans_gene_strand)

    #build the interval tree for 3utr_alignment
    tree_index2 = build_genomic_interval_tree(psl_hits_3utr)
    
    predict_3utr_dict = defaultdict(str)
     

    predicted_name_list = []
    # start loop the psl file.
    for Chr, gene_list in psl_hits.items():
        for hit in gene_list:
        
            # The hit gene should be the best homology in the target speceis.
            hited_gene = overlaps_gene_region(target_annotation, tree_index, hit, Chr)        
            
            if not hited_gene:
                continue                        

            # search the 3utr information
            
            gene_3utr = target_genes[hited_gene["gene_name"]]
            
            # predict use AATAAA
            seq = extract_truncated_sequence(target_genome, gene_3utr["chrom"], gene_3utr["utr_start"], 
                                                gene_3utr["utr_end"], gene_3utr["strand"])
            end = find_polyA(seq=seq)
            seq2 = seq[:end]
            if seq == seq2:
                id = hited_gene["gene_name"]  +  "_1kb"
            else:
                id = hited_gene["gene_name"]  + "_AAUAAA"
                seq = seq2
            
            if id not in predicted_name_list:

                with open(args.output_fasta, "a") as out:
                    out.write(f">{id}\n{seq}\n")
                out.close()
            
            predicted_name_list.append(id)

            # find psl overlap overlap the 3utr alignment with this target gene and do polish. 

            hits_3utrs = overlaps_3prime_region(psl_hits_3utr, tree_index2, gene_3utr)

            length = []
            select_hits = []
            if hits_3utrs:
                for hit_3utr in hits_3utrs:
                    # if they hit the same gene which means this is one to one relationship.
                    # here should pick the longest 3utr hits from the same query hit. 
                    
                    if hit_3utr["celegans_gene_id"] == hit["celegans_gene_id"]:
                        select_hits.append(hit_3utr)

                        align_length = abs(hit_3utr["tStart"] - hit_3utr["tEnd"])
                        length.append(align_length)
    
            if select_hits:
                max_index = max(enumerate(length), key = lambda x:x[1])[0]
                select_hit = select_hits[max_index]
                
                
                #write the out put

                if gene_3utr["strand"] == 1:
                    if select_hit["gff_end"] <= gene_3utr["utr_end"]:
                        seq2 = extract_truncated_sequence(target_genome, gene_3utr["chrom"], gene_3utr["utr_start"], select_hit["gff_end"], gene_3utr["strand"])
                elif gene_3utr["strand"] == -1: 
                     if select_hit["gff_start"] >= gene_3utr["utr_start"]:
                        seq2 = extract_truncated_sequence(target_genome, gene_3utr["chrom"], select_hit["gff_start"], gene_3utr["utr_end"], gene_3utr["strand"])

                try:
                    id = hited_gene["gene_name"] + "_polishedby_" + hit["celegans_gene_id"]
                
                    # write the out put for polished one. 
                    with open(args.output_fasta, "a") as out:
                        out.write(f">{id}\n{seq2}\n")
                    out.close() 
                except NameError:
                    continue


main()
            
    

