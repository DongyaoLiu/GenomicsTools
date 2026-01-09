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

# ------------------------------------------------
# memory monitor
# ------------------------------------------------

def print_memory_usage(label=""):
    process = psutil.Process(os.getpid())
    print(f"{label}: {process.memory_info().rss / 1024 / 1024:.2f} MB")

# ------------------------------------------------
# extract the sequence inclusive bi-side base on the 1-based index
# ------------------------------------------------

def extract_truncated_sequence(
    genome_dict: Dict[str, SeqRecord],
    sequence_id: str,
    start: int,
    end: int,
    strand: int,
    use_zero_based: bool = False
) -> Optional[str]:
    """
    Extract a truncated sequence from a genome dictionary.
    If coordinates exceed boundaries, adjust to use the boundaries instead.

    Parameters:
    -----------
    genome_dict : Dict[str, SeqRecord]
        Dictionary containing SeqRecord objects (e.g., from Bio.SeqIO.to_dict())
    sequence_id : str
        ID of the sequence in the genome dictionary
    start : int
        Start position (1-based by default)
    end : int
        End position (inclusive)
    strand : int
        1 for forward strand, -1 for reverse strand
    use_zero_based : bool
        If True, positions are treated as 0-based (like Python indexing)
        If False, positions are treated as 1-based (like biological coordinates)

    Returns:
    --------
    Optional[str]
        The truncated DNA sequence as a string, or None if sequence_id not found
    """

    # Check if sequence_id exists in the dictionary
    if sequence_id not in genome_dict:
        print(f"Warning: Sequence ID '{sequence_id}' not found in genome dictionary")
        return None

    # Get the sequence record
    seq_record = genome_dict[sequence_id]
    seq_length = len(seq_record.seq)

    # Adjust coordinates if using 0-based indexing
    if use_zero_based:
        adj_start = start
        adj_end = end
    else:
        adj_start = start - 1  # Convert 1-based to 0-based
        adj_end = end  # End position is inclusive in slice

    # Adjust coordinates to fit within boundaries
    # Ensure start is at least 0 (or 1-based equivalent)
    if adj_start < 0:
        print(f"#Note: Start coordinate {start} exceeds lower boundary, using boundary instead")
        adj_start = 0

    # Ensure end doesn't exceed sequence length
    if adj_end > seq_length:
        print(f"#Note: End coordinate {end} exceeds sequence length {seq_length}, using boundary instead")
        adj_end = seq_length

    # Check if coordinates are valid after adjustment
    if adj_start >= adj_end:
        print(f"Warning: Invalid coordinates after boundary adjustment: start={adj_start}, end={adj_end}")
        # Return empty string if start >= end after adjustment
        return ""

    # Extract the subsequence
    subsequence = seq_record.seq[adj_start:adj_end]

    # Handle strand
    if strand == 1:
        result_seq = subsequence
    elif strand == -1:
        result_seq = subsequence.reverse_complement()
    else:
        print(f"#Warning: Invalid strand value {strand}. Must be 1 or -1")
        return None

    return str(result_seq)


#--------------------------------------------------
# Custmized function. header mapper
# Chr name in gff file is a little different from psl file
#--------------------------------------------------
#
#def header_mapper(rename_file, species):
#    names_dict = {}
#    with open(rename_file, "r") as f:
#        for line in f:
#            elements = line.strip().split()
#            if len(elements) != 3:
#                print(f"{line} doesn't has three elements")
#                continue
#            elif re.search(species, elements[0]):
#                names_dict[elements[2]] = elements[1]
#    return names_dict
#

# ---------------------------------------------------
# Parse C. elegans 3' UTR BED → gene_id -> strand
# ---------------------------------------------------
def parse_celegans_utr_bed(bed_file):
    """
    BED format:
      chr, start, end, gene_id, ..., strand
    Returns: dict {gene_id: '+' or '-'}
    """
    gene_strand = {}
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 6:
                continue
            gene_id = cols[3]
            strand = cols[5]
            if strand in ('+', '-'):
                gene_strand[gene_id] = strand
    return gene_strand

# ---------------------------------------------------
# Parse target GFF → gene info
# ---------------------------------------------------
def target_genes_3utr(gff_file, header_dict=None, max_len=1000):
    if not re.search(r"(gtf|gff)", gff_file):
        gff_file = glob.glob(f"{gff_file}*")[0]

    genes = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')

            # add regular expression for both gtf and gff file.
            if cols[2] not in ["transcript", "mRNA"] or len(cols) < 9:
                continue
                                
            chrom = cols[0]
            if header_dict:
                chrom2 = header_dict[chrom]
            else:
                chrom2 = chrom

            start = int(cols[3]) # 1-based
            end = int(cols[4]) # 1-based
            strand_char = cols[6]
            strand = 1 if strand_char == '+' else -1

            gene_id = None
            for attr in cols[8].split(';'):
                if attr.startswith('ID='):
                    gene_id = attr[3:]
                    break
                if re.search(r"transcript_id \"(\S+)\"", attr):
                    #This is for gtf format extraction
                    gene_id = re.search(r"transcript_id \"(\S+)\"", attr).group(1)
                    break
                
            if strand == 1:
                utr_start = end + 1
                utr_end = utr_start + max_len
            else:
                utr_end = start - 1
                utr_start = utr_end - max_len

            genes[gene_id] = {
                "chrom": chrom2,
                "strand": strand,
                "strand_char": strand_char,
                "gene_start": start,
                "gene_end": end,
                "utr_start": utr_start,
                "utr_end": utr_end
            }
    return genes

# ---------------------------------------------------
# Parse PSL with custom gene ID and strand validation
# ---------------------------------------------------
def parse_psl(psl_file, celegans_gene_strand, stitch=True, min_match=20):
    """
    PSL format: [0] celegans_gene_id, [1-21] standard PSL
    Validate: PSL query strand == BED strand
    Output includes GFF coordinates (1-based, inclusive) in hit dictionary

    The hits dictionary,
        The key is the tName and the value is a list of records follow the psl
        file order. And each record is a inner dictionary which has following keys,


    The dual-index system,
        the dual-index system consists of two dictionary. the primary index and the helper
        couter for establishing the primary index.
        The primary index is a embedded dictionary the outer key is the tName same as the hits
        dictionary and the innner key is WBGeneID which is the bed file annotation when runinng the
        bed file based halliftover from cactus .hal file.

        With this dual-index system. I can quick group the records base on a C. elegans geneID and do stitch.

    The primer index is used for
    """
    hits = defaultdict(list)  # key: tName (chromosome)
    index = defaultdict(lambda: defaultdict(list)) #
    second_index = defaultdict(int) # a counter for each chromosome
    psl_line_count = 0
    hits_append = 0

    with open(psl_file) as f:
        for line in f:

            if line.startswith(('track', '#')) or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 22:
                continue

            psl_line_count += 1

            celegans_gene_id = cols[0]
            match = int(cols[1])
            strand_col = cols[9]      # e.g., "++"
            tName = cols[14]
            tStart = int(cols[16])    # 0-based start
            tEnd = int(cols[17])      # 0-based end (exclusive)
            tSize = int(cols[15])     # target size
            length = abs(tEnd - tStart) + 1
            pident = match/length

            q_str_psl = strand_col[0]   # query strand in PSL
            t_str_psl = strand_col[1]   # target strand in PSL

            # Validate query strand consistency
            # if not then change both strand for query and target
            true_strand = celegans_gene_strand[celegans_gene_id]
            if q_str_psl != true_strand:
                if q_str_psl == "-":  # strand mismatch
                    q_str_psl = "+"
                else:
                    q_str_psl = "-"
                if t_str_psl == "-":
                    t_str_psl = "+"
                else:
                    t_str_psl = "-"

            if t_str_psl == "-":
                strand = -1
            else:
                strand = 1

            if stitch:
                index[tName][celegans_gene_id].append(second_index[tName])
                second_index[tName] += 1

            hits[tName].append({
                "tStart": tStart,
                "tEnd": tEnd,
                "tStrand": t_str_psl,
                "celegans_gene_id": celegans_gene_id,
                "length": tEnd - tStart,
                "gff_start": tStart + 1, # 1-based and inclusive
                "gff_end": tEnd, #inclusive
                "gff_strand": t_str_psl,  # Same as target strand
                "strand":strand
            })
            hits_append += 1

    print(f"Processed : {str(psl_line_count)} line")
    print(f"the number of hits append into the hit dictionary {hits_append}")

    if stitch:

        # rebuild the index here
        clean_hits = X626(hits, index)

        return clean_hits

    else:
        return hits

# ---------------------------------------------------
# Check overlap in 3' region (strand-aware)
# ---------------------------------------------------
def overlaps_3prime_region(hits, index, gene_info, min_overlap = 20):
    # Check if chromosome exists in the trees
    if gene_info["chrom"] not in index:
        return []

    chrom_tree = index[gene_info["chrom"]]
    overlapping_intervals = chrom_tree.overlap(gene_info["utr_start"], gene_info["utr_end"])
    result = []
    for interval in overlapping_intervals:
        # Unpack stored data: (feature_index, feature_strand)
        hit_idx, hit_strand = interval.data

        # Calculate overlap
        overlap_start = max(gene_info["utr_start"], interval.begin)
        overlap_end = min(gene_info["utr_end"], interval.end)
        overlap_length = overlap_end - overlap_start

        # Skip if overlap is below threshold
        if overlap_length < min_overlap:
            continue

        # Check strand matching
        strand_match = (hit_strand == gene_info["strand"])


        # Skip if strand doesn't match when required
        if not strand_match:
            continue

        # Get the full feature data
        overlapped_hits = hits[gene_info["chrom"]][hit_idx]

        # Add to results
        result.append(overlapped_hits)

    return result

# ---------------------------------------------------
# part of stitch function (X626) This idea had never been
# generated by deepseek and qwen
# ---------------------------------------------------

def rolling_merge(candidate_list=None, gene_name=None, threshold = 30, strand=None):


    # init
    sorted_condidate = sorted(candidate_list, key=lambda x: x['gff_start'])
    start_list = [sorted_condidate[i]["gff_start"] for i in range(len(sorted_condidate))]
    end_list = [sorted_condidate[i]["gff_end"] for i in range(len(sorted_condidate))]

    if len(candidate_list) == 1:
        return candidate_list

    n = len(start_list)
    i = 0

    while i < (n-1):
        #print("The i" + str(i))
        #print("The n" + str(n))
        # situation 1 + 2
        if end_list[i] + threshold >= start_list[i+1]:
            start_list.pop(i+1)
            if end_list[i] >= end_list[i+1]:
                end_list.pop(i+1)
            else:
                end_list.pop(i)

            # The new merged one need compare with the next to check if need to merge it again.
            # Cause the situation 3 will constantly move forword by step=1, so here need backword by step=1
            i -= 1


        # situation 3
        i += 1

        #
        n = len(start_list)

    # Generate the clean out
    output = [{"celegans_gene_id": gene_name, "gff_start" : start_list[i], "gff_end" : end_list[i], "strand": strand} for i in range(len(start_list))]
    #print(output)
    return output

# --------------------------------------------------
# X626 function is for target sequence ligation. (stitch!)
# --------------------------------------------------

def X626(raw_dict, index, max_gap = 30):
    #print("The index dict")
    #print(index)
    #print("the raw stitich dictionary")
    #print(raw_dict)
    output_dict = defaultdict(list) # keep the same data structure as raw_dict
    for Chr, gene_list in index.items():
        for gene_name, gene_index in gene_list.items():
            Plus_stitch_candidate_list = [raw_dict[Chr][i] for i in gene_index if raw_dict[Chr][i]["tStrand"] == "+"]
            Minus_stitch_candidate_list = [raw_dict[Chr][i] for i in gene_index if raw_dict[Chr][i]["tStrand"] == "-"]
            if Plus_stitch_candidate_list:
                Plus_clean = rolling_merge(Plus_stitch_candidate_list, gene_name, threshold= max_gap, strand  = 1)
            else:
                Plus_clean = []
            if Minus_stitch_candidate_list:
                Minus_clean = rolling_merge(Minus_stitch_candidate_list, gene_name, threshold= max_gap, strand = -1)
            else:
                Minus_clean = []
            output_dict[Chr].extend(Plus_clean + Minus_clean)
    print("The merged psl hits")
    return output_dict


# ---------------------------------------------------
# PolyA finder
# ---------------------------------------------------

def find_polyA(seq=None, motifs=["AATAAA"]):
    s = seq.upper()
    i_list = []
    for m in motifs:
        i = s.find(m)
        if i != -1:
            i_list.append(i)

    if len(i_list) == 1:
        return(i + 6)
    elif len(i_list) == 2:
        i = min(i_list)
        return(i + 6)
    else:
        return len(seq)

# ---------------------------------------------------
# gene_name extraction 
# ---------------------------------------------------

def geneNameExtraction(gene):

    if re.search(r"(g[0-9]+)", gene):
        #AUGUSTUS
        return re.search(r"(g[0-9]+)", gene).group(1)
    elif re.search(r"(\S+[0-9]+)\.[tT]?[0-9]?", gene):
        #Other
        return re.search(r"(\S+[0-9]+)\.[tT]?[0-9]?", gene).group(1)
    elif re.search(r"WBGene", gene):
        return gene
    else:
        return gene 
        
    



# ---------------------------------------------------
# interval tree index for fast search
# ---------------------------------------------------

from typing import Dict, List, Optional, Any, Tuple
from intervaltree import Interval, IntervalTree


def build_genomic_interval_tree(
    features_by_chrom) -> Dict[str, IntervalTree]:
    """
    Build interval tree index for genomic features.

    Parameters:
    -----------
    features_by_chrom : Dict[str, List[Dict[str, int]]]
        Dictionary with chromosome keys, each containing a list of feature dictionaries.
        Each feature must have 'start', 'end' keys, and optionally 'strand' key.

    Returns:
    --------
    Dict[str, IntervalTree]
        Dictionary mapping chromosome names to IntervalTree objects.
        Each interval stores: (feature_index, strand)
    """


    # Dictionary to store interval trees for each chromosome
    chrom_trees = {}

    for chrom, features in features_by_chrom.items():
        # Create a new interval tree for this chromosome
        tree = IntervalTree()

        for idx, feature in enumerate(features):
            start = feature['gff_start']
            end = feature['gff_end']
            strand = feature.get('strand', 1)  # Default to forward strand if not specified

            # Validate coordinates
            if start >= end:
                print(f"Warning: Invalid feature coordinates on {chrom}: {start}-{end}")
                continue

            # Add interval to tree
            # Store both feature index and strand as interval data
            tree.addi(start, end, (idx, strand))

        # Store the tree for this chromosome
        chrom_trees[chrom] = tree

    return chrom_trees

# ---------------------------------------------------
# Main
# ---------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psl", required=True, help="Custom PSL: [celegans_gene_id] + 21 standard cols")
    parser.add_argument("--celegans-utr-bed", required=True, help="C. elegans 3' UTR BED (col4=geneID, col6=strand)")
    parser.add_argument("--target-gff", required=True, type=str)
    parser.add_argument("--target-genome", required=True)
    parser.add_argument("--orthofinder", required=True)
    parser.add_argument("--target-species", required=True, help="Target species name")
    parser.add_argument("--rename-file", required=False)
    parser.add_argument("--output-fasta", default="target_3utrs.fa")
    parser.add_argument("--stats", default="target_stats.tsv")
    args = parser.parse_args()


    # Load data
    # header_dict = header_mapper(args.rename_file, args.target_species)
    celegans_gene_strand = parse_celegans_utr_bed(args.celegans_utr_bed)
    target_genome = SeqIO.to_dict(SeqIO.parse(args.target_genome, "fasta"))
    target_genes = target_genes_3utr(args.target_gff)

    # psl import and build the tree
    psl_hits = parse_psl(args.psl, celegans_gene_strand)
    print("the total psl hits")
    tree_index = build_genomic_interval_tree(psl_hits)

    # Parse OrthoFinder: build target_gene -> celegans_gene map
    df = pd.read_csv(args.orthofinder, sep="\t", index_col=0)

    target_to_celegans = defaultdict(list)
    ortho_dict = defaultdict(str)

    for _, row in df.iterrows():

        # for each row, they are in the same ortho group
        # should be WBGene\d+
        cegs = str(row.get("Caenorhabditis_elegans.proteins")).split(", ")
        
        # target genes
        tgs = str(row.get(f"{args.target_species}.proteins")).split(", ")
        # C. elegans genes and target species gene are in the same ortho group
        cegs = [g for g in cegs if g and g != None]
        tgs = [re.search(r"(\S+)", g).group(1) for g in tgs if g and re.search(r"(\S+)", g)]
        if tgs:
            for g in tgs:
                ortho_dict[g] = _
                if target_to_celegans[g]:
                    target_to_celegans[g] = cegs + target_to_celegans[g]
                else:
                    target_to_celegans[g] = cegs
    #print(ortho_dict)
    #print(target_to_celegans)    
    # Stats
    stats = defaultdict(int)

    polished = defaultdict(str)
    polished_count = defaultdict(int)


    # Begin the loop from the target species parsed gff file
    for gene_id, gene_info in target_genes.items():
        print_memory_usage(f"#loop {stats["total"]}")
        
        gene_id2 = geneNameExtraction(gene_id)
        print(gene_id)

        # init
        stats["total"] += 1
        print(gene_info)
        valid_hits = overlaps_3prime_region(psl_hits, tree_index, gene_info)

        if valid_hits:
            # homology with N2 protein
            # no homology with N2 protein
            print("#There are validated hits")
            for single_hit in valid_hits:
                # Find the ortholog information here
                # homology mapping

                # TODO:add more code for robust handdling of not match error
                egs_gene_list = target_to_celegans[gene_id2]
                print(f"# The single hit: {single_hit}")
                print(f"# The target species gene: {gene_id}, the gene info {gene_info}")

                # If I had confired the homology
                print("#The hit is " + single_hit["celegans_gene_id"])
                print("#The ortho group " + ".".join(egs_gene_list))

                if polished_count[gene_id]:
                    polished_count[gene_id] += 1
                else:
                    polished_count[gene_id] = 1

                record_title = f"{gene_id}_P{str(polished_count[gene_id])}_{single_hit["celegans_gene_id"]}"


                if single_hit["celegans_gene_id"] in egs_gene_list:
                    print(f"{record_title}\tpsl_homology_both\t{",".join(target_to_celegans[gene_id])}")
                    stats["psl_homology_both"] += 1
                    Orthogroup = ortho_dict[gene_id2]
                    record_title = record_title + f"_{Orthogroup}"
                else:
                    stats["psl_only_protein_divergent"] += 1
                    print(f"{record_title}\tpsl_only_protein_divergent\tNA")

                if gene_info["strand"] == 1:
                    polished[record_title] = {"chrom": gene_info["chrom"],
                                                "start": gene_info["utr_start"],
                                                "end": single_hit["gff_end"],
                                                "strand": gene_info["strand"],
                                                "predict": False}
                    if single_hit["gff_end"] > gene_info["utr_end"]:
                        polished[record_title]["predict"] = True
                else:
                    polished[record_title] = {"chrom": gene_info["chrom"],
                                              "start": single_hit["gff_start"],
                                              "end": gene_info["utr_end"],
                                              "strand": gene_info["strand"],
                                              "predict": False}
                    if single_hit["gff_start"] < gene_info["utr_start"]:
                        polished[record_title]["predict"] = True
                print(polished[record_title])
        else:
            #print(f"No psl hits match with {args.target_species} {gene_id}")
            if target_to_celegans[gene_id2]:
                Orthogroup = ortho_dict[gene_id2]
                out_id = gene_id + f"_{Orthogroup}"
                stats["protein_homology_only_3utr_divergent"] += 1
                print(f"{gene_id}\tprotein_homology_only_3utr_divergent\t{",".join(target_to_celegans[gene_id])}")
            else:
                out_id = gene_id + "_newgene"
                stats["new_gene_new_3utr"] += 1
                print(f"{gene_id}\tNOprotein_homology_only_3utr_divergent\tNA")

            polished[out_id] = {"chrom": gene_info["chrom"],
                                 "start": gene_info["utr_start"],
                                 "end": gene_info["utr_end"],
                                 "strand": gene_info["strand"],
                                 "predict": True}

    # Output
    with open(args.output_fasta, "w") as out:
        for id, coord in polished.items():
            id = args.target_species + "_" + id
            if coord["predict"]:
                seq = extract_truncated_sequence(target_genome, coord["chrom"], coord["start"], coord["end"], coord["strand"])
                end = find_polyA(seq=seq)
                seq2 = seq[:end]
                if seq == seq2:
                    id = id + "_1kb"
                else:
                    id = id + "_AAUAAA"
                    seq = seq2
            else:
                seq = extract_truncated_sequence(target_genome, coord["chrom"], coord["start"], coord["end"], coord["strand"])
            out.write(f">{id}\n{seq}\n")
    out.close()

    with open(args.stats, "w") as f:
        for k, v in stats.items():
            f.write(f"{args.target_species}\t{k}\t{v}\n")
    f.close()

if __name__ == "__main__":
    main()
