#!/usr/bin/env python3
"""
Script to parse a PSL file, extract sequences using start/end/strand info,
perform pairwise alignment with Biopython, and display in BLAST-like format.
"""
import argparse
import re
import warnings
from Bio import SeqIO
from Bio import Align
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

def parse_psl(psl_line):
    """
    Parses a single line from a PSL file.
    Returns a dictionary with field names as keys.
    PSL format: https://genome.ucsc.edu/FAQ/FAQformat.html#format2
    """
    fields = [
        'miRNA',
        'matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert',
        'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd',
        'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts'
    ]
    values = psl_line.strip().split('\t')

    if len(values) < 21:
        warnings.warn("psl file error")
        return None

    record = dict(zip(fields, values))

    int_fields = ['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert',
                  'tNumInsert', 'tBaseInsert', 'qSize', 'qStart', 'qEnd', 'tSize', 'tStart', 'tEnd', 'blockCount']
    for field in int_fields:
        record[field] = int(record[field])

    return record

def get_sequence(fasta_file, seq_id, start, end, size, strand='+'):
    """
    Extracts a sequence from a FASTA file based on ID and coordinates.
    Handles reverse complement if strand is '-'.
    """
    # WARNING: For large genomes, replace this with an indexed approach (e.g., pyfaidx)
    genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    if seq_id not in genome_dict:
        raise ValueError(f"Sequence ID {seq_id} not found in {fasta_file}")

    full_seq = genome_dict[seq_id].seq
    sub_seq = full_seq[start:end]

    if strand == '-':
        sub_seq = sub_seq.reverse_complement()

    return sub_seq



# wait for exploration
def X626(psl_dict, gap = 10):
    # no strand info be considered
    # {Chr start end strand}


    query_ID_list = []
    duplicated_ID_dict = {} # for duplicated ID

    for line_number, rec in psl_dict.items():
        miRNA_ID = rec['miRNA']
        if miRNA_ID in duplicated_ID_dict.keys():
            duplicated_ID_dict[miRNA_ID].append(rec)
        else:
            duplicated_ID_dict[miRNA_ID] = [rec]

    output_dict = {}
    #process major ligation function.
    for miRNA, rec_list in duplicated_ID_dict.items():
        if len(rec_list) >= 2:

            target_name_list = []
            target_start_list = []
            target_end_list = []
            target_strand_list = []

            query_name_list = []
            query_start_list = []
            query_end_list = []
            query_strand_list = []

            for i in range(len(rec_list)):

                target_name_list.append(rec_list[i]["tName"])
                target_start_list.append(rec_list[i]["tStart"])
                target_end_list.append(rec_list[i]["tEnd"])
                target_strand_list.append(rec_list[i]["strand"][1])

                query_name_list.append(rec_list[i]["qName"])
                query_start_list.append(rec_list[i]["qStart"])
                query_end_list.append(rec_list[i]["qEnd"])
                query_strand_list.append(rec_list[i]["strand"][0])

            uniq_tname = list(set(target_name_list))
            uniq_qname = list(set(query_name_list))

            plus_t_index = [index for index, element in enumerate(target_strand_list) if element == "+"]
            minus_t_index = [index for index, element in enumerate(target_strand_list) if element == "-"]
            plus_q_index = [index for index, element in enumerate(query_strand_list) if element == "+"]
            minus_q_index =[index for index, element in enumerate(query_strand_list) if element == "-"]



        else:
            # add a tag for Stitched record
            rec_list[0]["Stitched"] = "True"
            output_dict = rec_list[0] # remove the []

    return psl_dict


def aligment_analysis(align_obj, target, query, output, species, pattern, name, start, end, size, strand):
    alignment_str = str(align_obj)
    lines = alignment_str.strip().split('\n')

    target_line = lines[0]
    middle_line = lines[1]
    query_line = lines[2]

    matches = middle_line.count('|')
    mismatches = middle_line.count('.')  # Similar but not identical
    gaps = target_line.count('-')

    #for checking the coordinates
    target_aligned_start = align_obj.aligned[0][0][0]
    target_aligned_end = align_obj.aligned[0][-1][-1]
    target_aligned = target[target_aligned_start:target_aligned_end]

    #for checking the coordinates
    query_aligned_start = align_obj.aligned[1][0][0]
    query_aligned_end = align_obj.aligned[1][-1][1]
    query_aligned = query[query_aligned_start:query_aligned_end]


    #check the unaligned the head and tail and force to generate the full length of miRNA in the homogly region...in the target
    head = query_aligned_start
    tail = len(query) - query_aligned_end
    target_aligned_start2 = target_aligned_start - head
    target_aligned_end2 = target_aligned_end + tail
    if target_aligned_start2 < 0:
        target_aligned_start2 = 0
    if target_aligned_end2 > len(target):
        target_aligned_end2 = len(target)

    target_miRNA_predict = target[target_aligned_start2:target_aligned_end2]


    with open(output, "a") as file_out:
        file_out.write(f"{species}\t{pattern}\t{target_aligned}\t{query_aligned}\t{matches}\t{mismatches}\t{gaps}\t{target_miRNA_predict}\t{query}\t{name}\t{start}\t{end}\t{size}\t{strand}\n")
        file_out.close()



    matches = middle_line.count('|')
    positives = middle_line.count('.')  # Similar but not identical
    mismatches = len([i for i, (t, q) in enumerate(zip(target_line, query_line))
                     if t != q and t != '-' and q != '-' and middle_line[i] not in '|.'])
    gaps = target_line.count('-') + query_line.count('-')

    return matches, mismatches, gaps



def align2(psl_filename, query_fasta, target_fasta, miRNA_fasta, output):
    aligner = PairwiseAligner(match_score = 1, mismatch_score = -0.5, mode = "global")
    aligner.match_score = 2
    aligner.mismatch_score = -0.5
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2

    miRNA_dict = SeqIO.to_dict(SeqIO.parse(miRNA_fasta, "fasta"))
    print(target_fasta)
    species = re.search(r"Caenorhabditis_elegans_(.+).psl", psl_filename).group(1)


    # create the dictionary for the psl file, used for find the lines has same query region.
    psl_dict = {}
    with open(psl_filename, 'r') as psl_file:
        for line_num, line in enumerate(psl_file):
            #if line_num < 5:  # Skip PSL header lines
            #    continue
            if not line.strip():
                continue
            rec = parse_psl(line)
            if rec is None:
                print(f"Skipping malformed line: {line}")
                continue
            psl_dict[line_num] = rec

    psl_dict = X626(psl_dict, gap=10)

    for line_num, rec in psl_dict.items():

        print(f"\n{'='*80}")
        print(f"Processing PSL record #{line_num}")
        print(f"{'='*80}")
        print(f"miRNA : {rec['miRNA']}")
        print(f"Alignment: {rec['qName']} (Query) vs {rec['tName']} (Target)")
        print(f"PSL Score: {rec['matches']} matches, {rec['misMatches']} mismatches")
        print(f"Strand: {rec['strand']}")
        print(f"Query: {rec['qStart']}-{rec['qEnd']} (of {rec['qSize']})")
        print(f"Target: {rec['tStart']}-{rec['tEnd']} (of {rec['tSize']})")

        try:

            pattern = re.search(r'ID=([^_]+)_pre', rec['miRNA']).group(1)
            miRNA = {k: v for k, v in miRNA_dict.items() if re.search(f"{pattern}_[53]p", k)}


        except Exception as e:
            print(f"Error processing miRNA {e}")

        try:
            # Extract sequences using only start/end/strand info
            query_strand = rec['strand'][0]  # Gets the first character of strand
            target_strand = rec['strand'][1] # Gets the second character of strand

            #extract the sequence
            #pay attention to the strand information

            query_seq = get_sequence(query_fasta, rec['qName'], rec['qStart'], rec['qEnd'], rec['qSize'],query_strand)
            target_seq = get_sequence(target_fasta, rec['tName'], rec['tStart'], rec['tEnd'], rec['tSize'],target_strand)


            #remove too short alignment
            if len(query_seq) <= 18 or len(target_seq) <= 18:
                with open(output, "a") as file_out:
                    file_out.write(f"{species}\t{pattern}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                    file_out.close()
                continue



            print(">query_seq")
            print(f"{query_seq}")
            print(">target_seq")
            print(f"{target_seq}")

            # Perform pairwise alignment
            alignments = aligner.align(target_seq, query_seq)
            try:
                best_alignment = alignments[0]  # Take the top alignment
                print(best_alignment)
            except IndexError:
                with open(output, "a") as file_out:
                    file_out.write(f"{species}\t{pattern}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                    file_out.close()
                continue

            for k,v in miRNA.items():
                aligner.mode = "local"
                print(f"Query: {k}")
                print(f"Query aligned to the C. elegans region query region")
                print(aligner.align(query_seq, str(v.seq).replace('U', 'T').replace('u', 't'))[0])

                print(f"Query: {k} ")
                print(f"Query aligned to the Target: {target_fasta}")
                align_obj = aligner.align(target_seq, str(v.seq).replace('U', 'T').replace('u', 't'))
                try:
                    print(align_obj[0])
                    aligment_analysis(align_obj[0], target_seq, str(v.seq).replace('U', 'T').replace('u', 't'),
                                 output, species, pattern, rec['tName'], rec['tStart'], rec['tEnd'], rec['tSize'], target_strand)

                    print(f"Biopython Alignment Score: {best_alignment.score:.2f}")
                except IndexError:
                    file_out.write(f"{species}\t{pattern}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                    file_out.close()

        except Exception as e:
            print(f"Error processing record: {e}")
            import traceback
            traceback.print_exc()




def main(psl_filename, query_fasta, target_fasta, miRNA_fasta, output):
    """
    Main function to process the PSL file and generate alignments.
    """

    align2(psl_filename, query_fasta, target_fasta, miRNA_fasta, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
Script to parse a PSL file, extract sequences using start/end/strand info,
perform pairwise alignment with Biopython, and display in BLAST-like format.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''empty'''
     )
    parser.add_argument("psl_file", help='psl_file path')
    parser.add_argument("qfasta", type=str,
                        help="query fasta file")
    parser.add_argument("tfasta", type=str, help="target fasta file")
    parser.add_argument("mirfasta", type=str)
    parser.add_argument("output", type=str)

    args = parser.parse_args()
    main(args.psl_file,
         args.qfasta, args.tfasta, args.mirfasta, args.output)

