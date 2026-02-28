#!/usr/bin/env python3
"""
Process PSL files by grouping by query gene name and target chromosome,
and keeping only the entry with the most aligned nucleotides per group.
"""

import pandas as pd
import os
import glob
import argparse
import sys

def process_psl_file(input_file, output_folder, query_col=0, target_col=14, matches_col=1):
    """
    Process a single PSL file with customizable column indices.
    
    Parameters:
    input_file (str): Path to the input PSL file
    output_folder (str): Path to the output folder
    query_col (int): Column index for query gene name (0-based)
    target_col (int): Column index for target chromosome (0-based)
    matches_col (int): Column index for matches/aligned nucleotides (0-based)
    
    Returns:
    tuple: (input_filename, original_row_count, filtered_row_count, group_count)
    """
    
    # Define column names for PSL format (0-based, 21 columns)
    psl_columns = [
        'matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 
        'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 
        'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 
        'blockCount', 'blockSizes', 'qStarts', 'tStarts'
    ]
    
    try:
        # Read the PSL file
        df = pd.read_csv(input_file, sep='\t', header=None, names=psl_columns)
        
        # Get column names based on indices
        query_col_name = psl_columns[query_col]
        target_col_name = psl_columns[target_col]
        matches_col_name = psl_columns[matches_col]
        
        # Group by query and target columns
        grouped = df.groupby([query_col_name, target_col_name])
        
        # For each group, keep the row with the maximum matches
        filtered_rows = []
        
        for group_keys, group in grouped:
            # Find the row with maximum matches
            max_matches_idx = group[matches_col_name].idxmax()
            best_row = group.loc[max_matches_idx]
            filtered_rows.append(best_row)
        
        # Create filtered dataframe
        filtered_df = pd.DataFrame(filtered_rows)
        
        # Generate output filename
        output_file = os.path.join(output_folder, os.path.basename(input_file))
        
        # Save filtered dataframe
        filtered_df.to_csv(output_file, sep='\t', header=False, index=False)
        
        return (os.path.basename(input_file), len(df), len(filtered_df), len(grouped))
        
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}", file=sys.stderr)
        return (os.path.basename(input_file), 0, 0, 0)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Process PSL files by grouping and filtering based on aligned nucleotides',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i /path/to/input/folder -o /path/to/output/folder
  %(prog)s --input ./psl_files --output ./filtered_psl_files --query-col 0 --target-col 14 --matches-col 1
  %(prog)s -i ./data -o ./results --verbose
        
Column indices (0-based):
  Default: query_col=0 (first column), target_col=14 (15th column), matches_col=1 (second column)
  Adjust these if your PSL files have different column ordering
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', 
                        required=True,
                        help='Input folder containing PSL files')
    
    parser.add_argument('-o', '--output', 
                        required=True,
                        help='Output folder for filtered PSL files')
    
    # Optional arguments for column indices
    parser.add_argument('--query-col', 
                        type=int, 
                        default=0,
                        help='Column index for query gene name (0-based, default: 0)')
    
    parser.add_argument('--target-col', 
                        type=int, 
                        default=14,
                        help='Column index for target chromosome (0-based, default: 14)')
    
    parser.add_argument('--matches-col', 
                        type=int, 
                        default=1,
                        help='Column index for matches/aligned nucleotides (0-based, default: 1)')
    
    # Optional flags
    parser.add_argument('-v', '--verbose', 
                        action='store_true',
                        help='Print detailed processing information')
    
    parser.add_argument('--ext', 
                        default='.psl',
                        help='File extension to process (default: .psl)')
    
    parser.add_argument('--recursive', 
                        action='store_true',
                        help='Search for PSL files recursively in subfolders')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate input folder
    if not os.path.exists(args.input):
        print(f"Error: Input folder '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isdir(args.input):
        print(f"Error: '{args.input}' is not a directory.", file=sys.stderr)
        sys.exit(1)
    
    # Create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        if args.verbose:
            print(f"Created output folder: {args.output}")
    
    # Find all PSL files
    search_pattern = f"*{args.ext}"
    if args.recursive:
        search_path = os.path.join(args.input, "**", search_pattern)
        psl_files = glob.glob(search_path, recursive=True)
    else:
        search_path = os.path.join(args.input, search_pattern)
        psl_files = glob.glob(search_path)
    
    # Also check for uppercase extension if default .psl
    if args.ext == '.psl':
        uppercase_files = glob.glob(os.path.join(args.input, "*.PSL"))
        psl_files.extend(uppercase_files)
    
    if not psl_files:
        print(f"No {args.ext} files found in {args.input}")
        if args.recursive:
            print("Searched recursively including subfolders.")
        sys.exit(0)
    
    # Print configuration
    print(f"\n{'='*60}")
    print(f"PSL File Processor")
    print(f"{'='*60}")
    print(f"Input folder:  {args.input}")
    print(f"Output folder: {args.output}")
    print(f"File pattern:  *{args.ext}")
    print(f"Recursive:     {args.recursive}")
    print(f"\nColumn configuration (0-based indices):")
    print(f"  Query gene name column:    {args.query_col}")
    print(f"  Target chromosome column:  {args.target_col}")
    print(f"  Matches column:             {args.matches_col}")
    print(f"\nFound {len(psl_files)} files to process")
    print(f"{'='*60}\n")
    
    # Process files
    total_original = 0
    total_filtered = 0
    successful_files = 0
    failed_files = 0
    
    for psl_file in psl_files:
        result = process_psl_file(
            psl_file, 
            args.output,
            query_col=args.query_col,
            target_col=args.target_col,
            matches_col=args.matches_col
        )
        
        filename, original_count, filtered_count, group_count = result
        
        if original_count > 0:
            successful_files += 1
            total_original += original_count
            total_filtered += filtered_count
            
            if args.verbose:
                print(f"✓ {filename}:")
                print(f"    Original rows: {original_count}")
                print(f"    Groups found:  {group_count}")
                print(f"    Filtered rows: {filtered_count}")
                print(f"    Reduction:     {(1 - filtered_count/original_count)*100:.1f}%")
            else:
                print(f"✓ {filename} ({original_count} → {filtered_count} rows)")
        else:
            failed_files += 1
            if args.verbose:
                print(f"✗ {filename}: Failed to process")
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Processing Complete")
    print(f"{'='*60}")
    print(f"Successfully processed: {successful_files} files")
    if failed_files > 0:
        print(f"Failed to process:      {failed_files} files")
    print(f"Total original rows:    {total_original}")
    print(f"Total filtered rows:    {total_filtered}")
    if total_original > 0:
        reduction = (1 - total_filtered/total_original) * 100
        print(f"Overall reduction:      {reduction:.1f}%")
    print(f"Output files saved to: {args.output}")

if __name__ == "__main__":
    main()
