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
        start1 = size - end
        end1 = size - start
        sub_seq = sub_seq.reverse_complement()
        
    return sub_seq


def align2(psl_filename, query_fasta, target_fasta, miRNA_fasta, output):
    aligner = PairwiseAligner(match_score = 1, mismatch_score = -0.5, mode = "global")
    aligner.match_score = 2
    aligner.mismatch_score = -0.5
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2
    
    miRNA_dict = SeqIO.to_dict(SeqIO.parse(miRNA_fasta, "fasta"))
    print(target_fasta)
    species = re.search(r"(Cae[a-z_]+)", target_fasta).group(1)

    with open(psl_filename, 'r') as psl_file:
        for line_num, line in enumerate(psl_file):
            #if line_num < 5:  # Skip PSL header lines
            #    continue
            if not line.strip():
                continue
                
            print(f"\n{'='*80}")
            print(f"Processing PSL record #{line_num}")
            print(f"{'='*80}")
            
            rec = parse_psl(line)
            if rec is None:
                print(f"Skipping malformed line: {line}")
                continue
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
                        file_out.write(f"{species}\t{pattern}\tNA\n")
                    continue 
                

                
                print(">query_seq")
                print(f"{query_seq}")
                print(">target_seq")
                print(f"{target_seq}")
                
                # Perform pairwise alignment
                alignments = aligner.align(target_seq, query_seq)
                best_alignment = alignments[0]  # Take the top alignment
                print(best_alignment)
                
                for k,v in miRNA.items():
                    aligner.mode = "local"
                    print(f"Query: {k}")
                    print(f"Query aligned to the C. elegans region query region")
                    print(aligner.align(query_seq, str(v.seq).replace('U', 'T').replace('u', 't'))[0])
                    
                    print(f"Query: {k} ")
                    print(f"Query aligned to the Target: {target_fasta}")
                    print(aligner.align(target_seq, str(v.seq).replace('U', 'T').replace('u', 't'))[0])
                
            
                # Print alignment score information
                print(f"Biopython Alignment Score: {best_alignment.score:.2f}")
                
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




