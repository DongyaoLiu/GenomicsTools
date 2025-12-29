#!/usr/bin/env python3
"""
SubsetRepeats.py
Extract specific repeat families from a FASTA file based on a list of family IDs.
"""

import sys
import re
import argparse
from typing import Set, Dict

def read_fasta_to_dict(fasta_file: str) -> Dict[str, str]:
    """Read FASTA file into dictionary."""
    sequences = {}
    current_header = None
    current_sequence = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    if current_header is not None:
                        sequences[current_header] = ''.join(current_sequence)
                    current_header = line[1:]
                    current_sequence = []
                else:
                    current_sequence.append(line)
            
            if current_header is not None:
                sequences[current_header] = ''.join(current_sequence)
                
    except FileNotFoundError:
        print(f"Error: File {fasta_file} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    
    return sequences

def parse_repeat_id(header: str) -> str:
    """
    Parse the repeat family ID from the header.
    Only returns the numeric part after R=.
    """
    # Use re.match to match from the beginning of the string
    # Pattern: R= followed by digits only
    match = re.match(r'^R=(\d+)', header)
    if match:
        return match.group(1)
    
    # Try if the pattern is not at the beginning
    match = re.search(r'(?:^|:)R=(\d+)', header)
    if match:
        return match.group(1)
    
    # If no match, return the entire header up to first colon
    return header.split(':')[0] if ':' in header else header

def read_family_list(list_file: str) -> Set[str]:
    """
    Read list of family IDs from file.
    Supports formats like: R=0, R=1, or just 0, 1
    """
    families = set()
    
    try:
        with open(list_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Extract just the numeric part
                # Remove R= prefix if present
                if line.startswith('R='):
                    # Extract digits after R=
                    match = re.match(r'^R=(\d+)$', line)
                    if match:
                        families.add(match.group(1))
                    else:
                        # Try to extract any digits
                        digits = re.search(r'\d+', line)
                        if digits:
                            families.add(digits.group())
                else:
                    # Assume it's just a number
                    match = re.match(r'^(\d+)$', line)
                    if match:
                        families.add(match.group())
                    else:
                        print(f"Warning: Unrecognized format in list: {line}")
    
    except FileNotFoundError:
        print(f"Error: File {list_file} not found.")
        sys.exit(1)
    
    return families

def extract_families(sequences: Dict[str, str], family_ids: Set[str]) -> Dict[str, str]:
    """
    Extract sequences belonging to specified families.
    
    Args:
        sequences: Dictionary of all sequences
        family_ids: Set of family IDs to extract (as strings, e.g., "0", "1")
    
    Returns:
        Dictionary of extracted sequences
    """
    extracted = {}
    found_families = set()
    
    for header, sequence in sequences.items():
        family_id = parse_repeat_id(header)
        
        if family_id in family_ids:
            extracted[header] = sequence
            found_families.add(family_id)
    
    # Report missing families
    missing_families = family_ids - found_families
    if missing_families:
        print(f"Warning: The following families were not found in the FASTA file:")
        for missing in sorted(missing_families):
            print(f"  R={missing}")
    
    return extracted

def write_fasta(sequences: Dict[str, str], output_file: str):
    """Write sequences to FASTA file."""
    with open(output_file, 'w') as f:
        for header, sequence in sequences.items():
            f.write(f'>{header}\n')
            # Write in chunks of 60 characters
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')

def main():
    parser = argparse.ArgumentParser(
        description='Extract specific repeat families from a FASTA file.'
    )
    parser.add_argument(
        'input_fasta',
        help='Input FASTA file containing repeat sequences'
    )
    parser.add_argument(
        'family_list',
        help='File containing list of family IDs to extract (one per line)'
    )
    parser.add_argument(
        '-o', '--output',
        default='subset_repeats.fasta',
        help='Output FASTA file for extracted sequences (default: subset_repeats.fasta)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    parser.add_argument(
        '--include-family-id',
        action='store_true',
        help='Include R=family_id in the output headers'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Reading FASTA file: {args.input_fasta}")
    
    # Read all sequences
    sequences = read_fasta_to_dict(args.input_fasta)
    
    if args.verbose:
        print(f"Total sequences read: {len(sequences)}")
    
    # Read family list
    if args.verbose:
        print(f"Reading family list: {args.family_list}")
    
    family_ids = read_family_list(args.family_list)
    
    if args.verbose:
        print(f"Families to extract: {len(family_ids)}")
        print(f"Family IDs: {', '.join(sorted(family_ids))}")
    
    # Extract sequences
    extracted = extract_families(sequences, family_ids)
    
    if args.verbose:
        print(f"\nExtraction results:")
        print(f"  Sequences extracted: {len(extracted)}")
        
        # Count by family
        family_counts = {}
        for header in extracted.keys():
            family_id = parse_repeat_id(header)
            family_counts[family_id] = family_counts.get(family_id, 0) + 1
        
        print(f"  Distribution by family:")
        for family_id, count in sorted(family_counts.items()):
            print(f"    R={family_id}: {count} sequences")
    
    # Write output
    write_fasta(extracted, args.output)
    
    if args.verbose:
        print(f"\nOutput written to: {args.output}")

if __name__ == "__main__":
    main()
