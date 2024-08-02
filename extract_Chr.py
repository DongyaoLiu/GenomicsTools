import re
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='The simple script is designed for extract one or more complate Scaffold/Chr')

# Add the file arguments
parser.add_argument('-fasta', type=str, help='Path to the fasta file')
parser.add_argument('-Chr_list', type=str, help='Path to the Chr list that need to extract')
parser.add_argument('-fasta_fold', type=str, help='Path to the fasta folder')

# Parse the command-line arguments
args = parser.parse_args()

def extract_Chr(parsed_fasta, id):
	for sequence in parsed_fasta:
		if sequence.id == id:
			return sequence
			break

def write_Chr(select_sequence, filename):
	with open(filename, 'w') as f:
		SeqIO.write(select_sequence, f, 'fasta')

def parse_Chrlist(Chrlist):
	ChrDict = {}
	for Chr in Chrlist.readlines():
		line_elements = Chr.split("\t")
		Genome = f"{args.fasta_fold}/{line_elements[0]}"
		Chr = [line_elements[1]]
		if ChrDict.has_key(Genome):
			TmpChr = ChrDict[Genome] + Chr
			ChrDict[Genome] = TmpChr
		else:
			ChrDict[Genome] = Chr
	return ChrDict


"""
The Chr list for multiple genome extraction.

Col1	Col2
Genomename	Chr

"""


if args.fasta_fold is not None:
	GenomeChr = open(args.fasta_fold, "r")
	GenomeChrDict = parse_Chrlist(GenomeChr)
	for Genome in GenomeChrDict.keys():
		StranName = re.match(r"^(.*)\.(fa|fasta)$", os.path.basename(Genome))
		GenomeParser = SeqIO.parse(Genome, 'fasta')
		for Chr in GenomeChrDict[Genome]:
			SelectSeq = extract_Chr(GenomeParser, Chr)
			SeqFileName = f"{StranName.group(1)}_{Chr}.fa"
			write_Chr(SelectSeq, SeqFileName)
elif args.fasta is not None:
	sequences = SeqIO.parse(args.fasta, 'fasta')
	print("under development")
else:
	print("fasta fold path or fasta didn't specify")

