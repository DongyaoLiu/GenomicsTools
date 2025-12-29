#!/usr/bin/env python3
import os
import re
import argparse
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description='Standardize FASTA headers in a folder')
    parser.add_argument('input_folder', help='Path to folder containing FASTA files')
    parser.add_argument('output_folder', help='Path to output folder for renamed files')
    parser.add_argument('--table', default='rename_table.tsv', 
                       help='Name of rename table file (default: rename_table.tsv)')
    return parser.parse_args()

def extract_chromosome_info(header):
    """Extract chromosome information from header"""
    # Try to find chromosome pattern
    chromosome_patterns = [
        r'chromosome\s*[:\s]\s*([IVX0-9]+)',  # chromosome: I, chromosome: 1, etc.
        r'chr\s*[:\s]\s*([IVX0-9]+)',        # chr: I, chr: 1, etc.
        r'chromosome\s*=\s*([IVX0-9]+)',     # chromosome=I
        r'chr\s*=\s*([IVX0-9]+)',            # chr=I
    ]
    
    for pattern in chromosome_patterns:
        match = re.search(pattern, header, re.IGNORECASE)
        if match:
            chrom = match.group(1).strip()
            # Convert numeric chromosomes to Roman numerals if needed
            if chrom.isdigit():
                chrom_num = int(chrom)
                if chrom_num <= 10:
                    # Simple conversion for small numbers
                    roman_numerals = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']
                    return roman_numerals[chrom_num - 1]
            return chrom
    
    # Check for mitochondrion
    if re.search(r'mitochondrion|mtDNA|mitochondrial', header, re.IGNORECASE):
        return 'MtDNA'
    
    return None

def process_fasta_file(input_file, output_file):
    """Process a single FASTA file and return rename table entries"""
    rename_entries = []
    scaffold_count = 1
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Original header
                original_header = line[1:].strip()
                
                # Extract cactus_name (first word before space)
                cactus_name = original_header.split()[0] if ' ' in original_header else original_header
                
                # Extract chromosome information
                chrom_info = extract_chromosome_info(original_header)
                
                if chrom_info:
                    # Use chromosome-based naming
                    output_header = chrom_info
                else:
                    # Use scaffold naming
                    output_header = f'Scaffold{scaffold_count}'
                    scaffold_count += 1
                
                # Write new header
                outfile.write(f'>{output_header}\n')
                
                # Add to rename table
                rename_entries.append((original_header, cactus_name, output_header))
            else:
                # Write sequence lines as-is
                outfile.write(line)
    
    return rename_entries

def main():
    args = parse_arguments()
    
    # Create output folder if it doesn't exist
    Path(args.output_folder).mkdir(parents=True, exist_ok=True)
    
    # Initialize rename table
    all_rename_entries = []
    
    # Process all FASTA files in input folder
    input_path = Path(args.input_folder)
    fasta_files = list(input_path.glob('*.fasta')) + list(input_path.glob('*.fa')) + \
                  list(input_path.glob('*.fna')) + list(input_path.glob('*.fas'))
    
    if not fasta_files:
        print(f"No FASTA files found in {args.input_folder}")
        return
    
    print(f"Found {len(fasta_files)} FASTA files to process")
    
    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file.name}")
        
        # Define output file path (same name, different folder)
        output_file = Path(args.output_folder) / fasta_file.name
        
        # Process the file
        rename_entries = process_fasta_file(fasta_file, output_file)
        
        # Add entries with file context
        for entry in rename_entries:
            all_rename_entries.append((fasta_file.name, *entry))
    
    # Write rename table
    with open(args.table, 'w') as table_file:
        # Write header
        table_file.write("filename\toriginal_header\tcactus_name\toutput_header\n")
        
        # Write entries
        for entry in all_rename_entries:
            table_file.write("\t".join(entry) + "\n")
    
    print(f"\nProcessing complete!")
    print(f"Renamed files saved to: {args.output_folder}")
    print(f"Rename table saved to: {args.table}")
    print(f"Total sequences processed: {len(all_rename_entries)}")

if __name__ == "__main__":
    main()
