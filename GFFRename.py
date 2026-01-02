#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
import pandas as pd
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Update GFF files using rename table')
    parser.add_argument('rename_table', help='Path to rename table TSV file')
    parser.add_argument('gff_folder', help='Path to folder containing GFF files')
    parser.add_argument('output_folder', help='Path to output folder for updated GFF files')
    parser.add_argument('--gff_extension', default='.gff', 
                       help='GFF file extension (default: .gff)')
    return parser.parse_args()

def create_header_mapping(rename_table_path):
    """Create mapping dictionaries from rename table"""
    # Read rename table
    df = pd.read_csv(rename_table_path, sep='\t')
    
    # Create mapping: {filename: {original_header: output_header}}
    modified_mapping = defaultdict()
    
    # Create mapping: {cactus_name: output_header} for fallback
    cactus_mapping = defaultdict()
    
    for _, row in df.iterrows():
        # cactus_name	output_header	header_modified
        header_modified = row['header_modified']
        cactus_name = row['cactus_name']
        output_header = row['output_header']
        
        cactus_mapping[cactus_name] = output_header 
        # Also add cactus_name mapping
        modified_mapping[header_modified] = output_header
    
    return [cactus_mapping, modified_mapping]

def extract_seqid_from_gff_line(line):
    """Extract sequence ID from GFF line"""
    parts = line.strip().split('\t')
    if len(parts) >= 9:
        return parts[0]
    return None

def update_gff_file(gff_file, mappings, output_file):
    """Update a single GFF file using the mappings"""
    filename = gff_file.name
    
    with open(gff_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip comment/header lines
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            # Try to extract sequence ID
            seqid = extract_seqid_from_gff_line(line)
            
            if seqid:
                
                if mappings[0][seqid]:
                    updated_seqid = mappings[0][seqid]
                elif mappings[1][seqid]:
                    updated_seqid = mappings[1][seqid]
                else:
                    raise ValueError(f"The seqid :{seqid} could not found in any transfer name dictionary")
                parts = line.split('\t')
                parts[0] = updated_seqid
                line = '\t'.join(parts)
                           
            outfile.write(line)

def main():
    args = parse_arguments()
    
    # Create output folder if it doesn't exist
    Path(args.output_folder).mkdir(parents=True, exist_ok=True)
    
    # Load rename table and create mappings
    print(f"Loading rename table: {args.rename_table}")
    mappings = create_header_mapping(args.rename_table)
    
    # Get list of GFF files
    gff_folder = Path(args.gff_folder)
    gff_files = list(gff_folder.glob(f'*{args.gff_extension}'))
    
    # Also look for other common GFF extensions
    if not gff_files:
        gff_files = list(gff_folder.glob('*.gff3')) + list(gff_folder.glob('*.gff2')) + list(gff_folder.glob('*.gff'))
    
    if not gff_files:
        print(f"No GFF files found in {args.gff_folder}")
        print(f"Looking for extension: gff3, gff2, gff")
        return
    
    print(f"Found {len(gff_files)} GFF files to process")
    
    # Process each GFF file
    stats = {'updated': 0, 'skipped': 0, 'total_files': 0}
    
    for gff_file in gff_files:
        print(f"Processing: {gff_file.name}")
        
        
        # Define output file path
        output_file = Path(args.output_folder) / gff_file.name
        
        # Update the GFF file
        update_gff_file(gff_file, mappings, output_file)
        print(f"Successful: {gff_file.name}")
    print(f"\nProcessing complete!")
    print(f"Output saved to: {args.output_folder}")
    

if __name__ == "__main__":
    main()
