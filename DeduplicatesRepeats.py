#!/usr/bin/env python3
"""
DeduplicateRepeats.py
Remove identical and substring sequences from repeat FASTA file.
"""

import sys
import re
from collections import defaultdict
from typing import Dict, List, Set, Tuple
import argparse

def read_fasta_to_dict(fasta_file: str) -> Dict[str, str]:
    """Read FASTA file into dictionary."""
    sequences = {}
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_sequence).upper()
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)
        
        if current_header is not None:
            sequences[current_header] = ''.join(current_sequence).upper()
    
    return sequences

def parse_repeat_id(header: str) -> str:
    """Parse repeat family ID from header."""
    match = re.search(r'R=(\d+)', header)
    if match:
        return match.group(1)
    else:
        return header.split(':')[0] if ':' in header else header

def is_substring(seq1: str, seq2: str) -> bool:
    """
    Check if one sequence is a substring of another.
    Case-insensitive comparison.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Quick check: if lengths are equal, they can't be substrings unless identical
    if len(seq1) == len(seq2):
        return seq1 == seq2
    
    # Ensure seq1 is the longer one
    if len(seq1) < len(seq2):
        seq1, seq2 = seq2, seq1
    
    return seq2 in seq1

def deduplicate_family(sequences: Dict[str, str]) -> Dict[str, str]:
    """
    Remove identical sequences and perfect substrings within a family.
    Returns the longest unique sequences.
    """
    if not sequences:
        return {}
    
    # Convert to list for easier processing
    seq_items = list(sequences.items())
    
    # Group identical sequences
    unique_seqs = {}
    seq_to_headers = defaultdict(list)
    
    for header, seq in seq_items:
        seq_upper = seq.upper()
        seq_to_headers[seq_upper].append(header)
    
    # For each unique sequence, keep the first header
    for seq, headers in seq_to_headers.items():
        unique_seqs[headers[0]] = seq
    
    # Now remove perfect substrings
    deduplicated = {}
    seq_list = list(unique_seqs.items())
    n = len(seq_list)
    
    # Sort by length descending
    seq_list.sort(key=lambda x: len(x[1]), reverse=True)
    
    # Create a set to track sequences to keep
    keep = [True] * n
    
    for i in range(n):
        if not keep[i]:
            continue
        
        header_i, seq_i = seq_list[i]
        
        for j in range(i + 1, n):
            if not keep[j]:
                continue
            
            header_j, seq_j = seq_list[j]
            
            # Check if seq_j is a substring of seq_i
            if is_substring(seq_i, seq_j):
                keep[j] = False
    
    # Collect kept sequences
    for i, (keep_flag, (header, seq)) in enumerate(zip(keep, seq_list)):
        if keep_flag:
            deduplicated[header] = seq
    
    return deduplicated

def deduplicate_all_families(sequences: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, Dict]]:
    """
    Deduplicate sequences across all families.
    
    Returns:
        Tuple of (deduplicated_sequences, statistics_dict)
    """
    # Separate by family first
    families = defaultdict(dict)
    for header, seq in sequences.items():
        family_id = parse_repeat_id(header)
        families[family_id][header] = seq
    
    # Deduplicate each family
    deduplicated = {}
    stats = {}
    
    for family_id, family_seqs in families.items():
        original_count = len(family_seqs)
        dedup_family = deduplicate_family(family_seqs)
        dedup_count = len(dedup_family)
        
        # Add to main deduplicated dict
        deduplicated.update(dedup_family)
        
        # Record statistics
        stats[family_id] = {
            'original': original_count,
            'deduplicated': dedup_count,
            'removed': original_count - dedup_count,
            'reduction_percent': ((original_count - dedup_count) / original_count * 100) if original_count > 0 else 0
        }
    
    return deduplicated, stats

def write_deduplicated_fasta(sequences: Dict[str, str], output_file: str):
    """Write deduplicated sequences to FASTA file."""
    with open(output_file, 'w') as f:
        for header, seq in sequences.items():
            f.write(f'>{header}\n')
            # Write in chunks of 60 characters
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

def write_statistics(stats: Dict[str, Dict], output_file: str):
    """Write deduplication statistics to file."""
    with open(output_file, 'w') as f:
        f.write("Deduplication Statistics\n")
        f.write("=" * 50 + "\n\n")
        
        # Sort by original count descending
        sorted_stats = sorted(stats.items(), 
                            key=lambda x: x[1]['original'], 
                            reverse=True)
        
        total_original = sum(s['original'] for s in stats.values())
        total_dedup = sum(s['deduplicated'] for s in stats.values())
        total_removed = total_original - total_dedup
        
        f.write("Per Family Statistics:\n")
        for family_id, family_stats in sorted_stats:
            f.write(f"Family R={family_id}:\n")
            f.write(f"  Original sequences: {family_stats['original']}\n")
            f.write(f"  After deduplication: {family_stats['deduplicated']}\n")
            f.write(f"  Removed: {family_stats['removed']} ({family_stats['reduction_percent']:.1f}%)\n")
            f.write("\n")
        
        f.write("\nOverall Statistics:\n")
        f.write(f"  Total original sequences: {total_original}\n")
        f.write(f"  Total after deduplication: {total_dedup}\n")
        f.write(f"  Total removed: {total_removed} ({total_removed/total_original*100:.1f}%)\n")
        f.write(f"  Number of families: {len(stats)}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Remove identical and substring sequences from repeat FASTA file.'
    )
    parser.add_argument(
        'input_fasta',
        help='Input FASTA file containing repeat sequences'
    )
    parser.add_argument(
        '-o', '--output',
        default='deduplicated_repeats.fasta',
        help='Output FASTA file for deduplicated sequences (default: deduplicated_repeats.fasta)'
    )
    parser.add_argument(
        '-s', '--stats',
        default='deduplication_stats.txt',
        help='Output statistics file (default: deduplication_stats.txt)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output to console'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Reading input FASTA: {args.input_fasta}")
    
    # Read sequences
    sequences = read_fasta_to_dict(args.input_fasta)
    
    if args.verbose:
        print(f"Total sequences read: {len(sequences)}")
    
    # Deduplicate
    deduplicated, stats = deduplicate_all_families(sequences)
    
    if args.verbose:
        print(f"Sequences after deduplication: {len(deduplicated)}")
        print(f"Removed: {len(sequences) - len(deduplicated)} sequences")
    
    # Write output
    write_deduplicated_fasta(deduplicated, args.output)
    write_statistics(stats, args.stats)
    
    if args.verbose:
        print(f"\nOutput files created:")
        print(f"  Deduplicated sequences: {args.output}")
        print(f"  Statistics: {args.stats}")
        
        # Show top 5 families by reduction
        top_reductions = sorted(stats.items(), 
                              key=lambda x: x[1]['reduction_percent'], 
                              reverse=True)[:5]
        
        print("\nTop 5 families by reduction percentage:")
        for family_id, family_stats in top_reductions:
            print(f"  Family R={family_id}: {family_stats['reduction_percent']:.1f}% reduction")

if __name__ == "__main__":
    main()
