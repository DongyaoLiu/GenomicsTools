#!/usr/bin/env python3
"""
Protein Sequence Filter Script
Removes sequences containing stop codons (* or .) from FASTA files
"""

import os
import sys
import argparse
from pathlib import Path
from collections import Counter
from datetime import datetime
import textwrap

def read_fasta(file_path):
    """Read a FASTA file and return dictionary of sequences"""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]  # Take only the first word as header
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        # Don't forget the last sequence
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences

def contains_stop_codon(sequence):
    """Check if sequence contains stop codon symbols (* or .)"""
    return '*' in sequence or '.' in sequence

def process_fasta_file(input_path, output_folder, removed_records):
    """Process a single FASTA file"""
    filename = os.path.basename(input_path)
    output_path = os.path.join(output_folder, filename)
    
    sequences = read_fasta(input_path)
    
    kept_sequences = {}
    file_removed = []
    
    for header, seq in sequences.items():
        if contains_stop_codon(seq):
            # Record removed sequence information
            stop_positions = []
            if '*' in seq:
                stop_positions.extend([i for i, char in enumerate(seq) if char == '*'])
            if '.' in seq:
                stop_positions.extend([i for i, char in enumerate(seq) if char == '.'])
            
            removed_records.append({
                'filename': filename,
                'header': header,
                'sequence_length': len(seq),
                'stop_positions': stop_positions,
                'stop_count': len(stop_positions),
                'sequence': seq if len(seq) <= 100 else seq[:100] + "..."
            })
            file_removed.append(header)
        else:
            kept_sequences[header] = seq
    
    # Write filtered sequences to output file
    with open(output_path, 'w') as f:
        for header, seq in kept_sequences.items():
            f.write(f'>{header}\n')
            # Write sequence in chunks of 80 characters (standard FASTA format)
            for i in range(0, len(seq), 80):
                f.write(f'{seq[i:i+80]}\n')
    
    return len(sequences), len(kept_sequences), len(file_removed), file_removed

def analyze_removed_sequences(removed_records):
    """Perform statistical analysis on removed sequences"""
    if not removed_records:
        return "No sequences were removed."
    
    analysis = []
    
    # Basic statistics
    total_removed = len(removed_records)
    files_affected = len(set(rec['filename'] for rec in removed_records))
    
    analysis.append("=" * 60)
    analysis.append("REMOVED SEQUENCES STATISTICAL ANALYSIS")
    analysis.append("=" * 60)
    analysis.append(f"Total sequences removed: {total_removed}")
    analysis.append(f"Files affected: {files_affected}")
    analysis.append("")
    
    # Stop codon statistics
    stop_symbols = []
    for rec in removed_records:
        if '*' in rec['sequence']:
            stop_symbols.extend(['*'] * rec['sequence'].count('*'))
        if '.' in rec['sequence']:
            stop_symbols.extend(['.'] * rec['sequence'].count('.'))
    
    stop_counts = Counter(stop_symbols)
    analysis.append("Stop codon distribution:")
    for symbol, count in stop_counts.items():
        symbol_name = "asterisk (*)" if symbol == '*' else "period (.)"
        analysis.append(f"  {symbol_name}: {count} occurrences ({count/total_removed:.1%} of removed sequences)")
    
    # Sequence length analysis
    lengths = [rec['sequence_length'] for rec in removed_records]
    avg_length = sum(lengths) / len(lengths)
    min_length = min(lengths)
    max_length = max(lengths)
    
    analysis.append("")
    analysis.append("Sequence length statistics:")
    analysis.append(f"  Average length: {avg_length:.1f} amino acids")
    analysis.append(f"  Minimum length: {min_length} amino acids")
    analysis.append(f"  Maximum length: {max_length} amino acids")
    
    # Length distribution
    length_bins = {'<50': 0, '50-100': 0, '100-200': 0, '200-500': 0, '>500': 0}
    for length in lengths:
        if length < 50:
            length_bins['<50'] += 1
        elif length <= 100:
            length_bins['50-100'] += 1
        elif length <= 200:
            length_bins['100-200'] += 1
        elif length <= 500:
            length_bins['200-500'] += 1
        else:
            length_bins['>500'] += 1
    
    analysis.append("")
    analysis.append("Length distribution:")
    for bin_name, count in length_bins.items():
        if count > 0:
            analysis.append(f"  {bin_name}: {count} sequences ({count/total_removed:.1%})")
    
    # Files with most removals
    file_counts = Counter(rec['filename'] for rec in removed_records)
    analysis.append("")
    analysis.append("Files with most removed sequences:")
    for filename, count in file_counts.most_common(5):
        analysis.append(f"  {filename}: {count} sequences")
    
    return '\n'.join(analysis)

def generate_report(processing_results, removed_records, output_folder, analysis_text):
    """Generate a comprehensive report"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    report = []
    report.append("=" * 70)
    report.append("PROTEIN SEQUENCE FILTERING REPORT")
    report.append("=" * 70)
    report.append(f"Date: {timestamp}")
    report.append(f"Output folder: {output_folder}")
    report.append("")
    
    # Summary table
    report.append("PROCESSING SUMMARY")
    report.append("-" * 70)
    report.append(f"{'Filename':<30} {'Total':<10} {'Kept':<10} {'Removed':<10}")
    report.append("-" * 70)
    
    total_sequences = 0
    total_kept = 0
    total_removed = 0
    
    for result in processing_results:
        filename, total, kept, removed = result
        total_sequences += total
        total_kept += kept
        total_removed += removed
        report.append(f"{filename:<30} {total:<10} {kept:<10} {removed:<10}")
    
    report.append("-" * 70)
    report.append(f"{'TOTALS':<30} {total_sequences:<10} {total_kept:<10} {total_removed:<10}")
    report.append("")
    
    # Removed sequences list
    if removed_records:
        report.append("REMOVED SEQUENCES DETAILS")
        report.append("-" * 70)
        for i, rec in enumerate(removed_records, 1):
            report.append(f"{i}. File: {rec['filename']}")
            report.append(f"   Header: {rec['header']}")
            report.append(f"   Length: {rec['sequence_length']} aa, Stop codons: {rec['stop_count']}")
            if rec['stop_positions']:
                positions = ', '.join(str(pos) for pos in rec['stop_positions'][:5])
                if len(rec['stop_positions']) > 5:
                    positions += f", ... (+{len(rec['stop_positions'])-5} more)"
                report.append(f"   Positions: {positions}")
            report.append(f"   Sequence (first 100 chars): {rec['sequence']}")
            report.append("")
    
    # Add statistical analysis
    report.append(analysis_text)
    
    return '\n'.join(report)

def main():
    parser = argparse.ArgumentParser(
        description='Filter protein FASTA files by removing sequences containing stop codons (* or .)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
        Examples:
          %(prog)s /path/to/fasta/files --output filtered_proteins
          %(prog)s /path/to/fasta/files -o filtered_proteins --ext .faa
        ''')
    )
    
    parser.add_argument('input_folder', help='Folder containing input FASTA files')
    parser.add_argument('-o', '--output', default='filtered_proteins',
                       help='Output folder name (default: filtered_proteins)')
    parser.add_argument('--ext', default='.fasta',
                       help='FASTA file extension (default: .fasta, supports: .fasta, .faa, .fa, .fas)')
    
    args = parser.parse_args()
    
    # Set up paths
    input_folder = Path(args.input_folder)
    output_folder = Path(args.output)
    
    if not input_folder.exists():
        print(f"Error: Input folder '{input_folder}' does not exist!")
        sys.exit(1)
    
    # Create output folder
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Find FASTA files
    extensions = ['.fasta', '.faa', '.fa', '.fas', args.ext]
    fasta_files = []
    for ext in extensions:
        fasta_files.extend(list(input_folder.glob(f'*{ext}')))
    
    if not fasta_files:
        print(f"No FASTA files found in '{input_folder}' with extensions: {extensions}")
        sys.exit(1)
    
    print(f"Found {len(fasta_files)} FASTA file(s) to process")
    print(f"Output will be saved to: {output_folder}")
    
    # Process files
    processing_results = []
    removed_records = []
    
    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file.name}...")
        total, kept, removed, file_removed = process_fasta_file(
            fasta_file, output_folder, removed_records
        )
        processing_results.append((fasta_file.name, total, kept, removed))
        print(f"  Processed: {total} total, {kept} kept, {removed} removed")
    
    # Perform statistical analysis
    analysis_text = analyze_removed_sequences(removed_records)
    
    # Generate report
    report = generate_report(
        processing_results, 
        removed_records, 
        output_folder, 
        analysis_text
    )
    
    # Save report
    report_file = output_folder / "filtering_report.txt"
    with open(report_file, 'w') as f:
        f.write(report)
    
    # Print summary
    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE")
    print("=" * 60)
    print(f"Total files processed: {len(fasta_files)}")
    print(f"Total sequences removed: {len(removed_records)}")
    print(f"Filtered files saved to: {output_folder}")
    print(f"Detailed report saved to: {report_file}")
    print("\nStatistical Summary:")
    print("-" * 60)
    for line in analysis_text.split('\n')[:15]:  # Print first part of analysis
        if line.strip():
            print(line)

if __name__ == '__main__':
    main()
