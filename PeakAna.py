#!/usr/bin/env python3
"""
AB1 File Analysis for Double Peaks in Diploid Sequencing
Analyzes .ab1 files to identify heterozygous positions and generate consensus/alternative sequences.
"""

import os
import sys
import argparse
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class AB1DoublePeakAnalyzer:
    def __init__(self, min_quality=20, min_peak_ratio=0.2):
        self.min_quality = min_quality
        self.min_peak_ratio = min_peak_ratio
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -1
        self.aligner.extend_gap_score = -0.5

    def read_reference_sequence(self, reference_input):
        """
        Read reference sequence from various input types:
        - FASTA file (.fasta, .fa, .fna)
        - GenBank file (.gb, .gbk)
        - Plain text file (.txt)
        - Direct string sequence
        """
        # If it's a file path
        if os.path.exists(reference_input):
            file_ext = os.path.splitext(reference_input)[1].lower()

            try:
                if file_ext in ['.fasta', '.fa', '.fna', '.fas']:
                    # Read FASTA file
                    records = list(SeqIO.parse(reference_input, "fasta"))
                    if not records:
                        raise ValueError("No sequences found in FASTA file")
                    return str(records[0].seq).upper()

                elif file_ext in ['.gb', '.gbk', '.genbank']:
                    # Read GenBank file
                    records = list(SeqIO.parse(reference_input, "genbank"))
                    if not records:
                        raise ValueError("No sequences found in GenBank file")
                    return str(records[0].seq).upper()

                elif file_ext in ['.txt', '.seq']:
                    # Read plain text file
                    with open(reference_input, 'r') as f:
                        content = f.read().strip()
                    # Remove header lines and whitespace
                    lines = [line for line in content.split('\n')
                            if not line.startswith('>') and line.strip()]
                    return ''.join(lines).upper()

                else:
                    raise ValueError(f"Unsupported file format: {file_ext}")

            except Exception as e:
                raise ValueError(f"Error reading reference file {reference_input}: {e}")

        else:
            # Treat as direct sequence string
            # Remove whitespace and convert to uppercase
            sequence = ''.join(reference_input.split()).upper()

            # Validate sequence contains only valid nucleotides
            valid_bases = set('ACGTN-')
            if not all(base in valid_bases for base in sequence):
                invalid_bases = set(sequence) - valid_bases
                raise ValueError(f"Invalid nucleotides in reference sequence: {invalid_bases}")

            return sequence

    def read_ab1_file(self, ab1_file):
        """Read AB1 file and extract sequence, quality scores, and trace data."""
        try:
            record = SeqIO.read(ab1_file, "abi")
        except Exception as e:
            raise ValueError(f"Error reading AB1 file: {e}")

        # Extract sequence and quality scores
        sequence = str(record.seq)
        quality_scores = record.letter_annotations.get("phred_quality", [])
        print(quality_scores)

        # Extract trace data for peak analysis
        trace_channels = {}
        print(record.annotations['abif_raw'].keys())
        channel_mapping = {
            'DATA9': 'A',
            'DATA10': 'C',
            'DATA11': 'G',
            'DATA12': 'T'
        }
        print()
        for channel, base in channel_mapping.items():
            if channel in record.annotations['abif_raw'].keys():
                trace_channels[base] = record.annotations['abif_raw'][channel]

        return {
            'sequence': sequence,
            'quality_scores': quality_scores,
            'trace_data': trace_channels,
            'record': record
        }

    def analyze_peaks_at_position(self, trace_data, position):
        """Analyze peaks at a given position to identify top two bases."""
        # Check if we have trace data for all bases
        required_bases = ['A', 'C', 'G', 'T']
        if not all(base in trace_data for base in required_bases):
            return None

        if position >= len(trace_data['A']):
            return None

        # Get intensity values for each base at this position
        base_intensities = {}
        for base in required_bases:
            if position < len(trace_data[base]):
                base_intensities[base] = trace_data[base][position]
            else:
                base_intensities[base] = 0

        # Sort bases by intensity
        sorted_bases = sorted(base_intensities.items(), key=lambda x: x[1], reverse=True)

        # Filter out bases with very low intensity
        max_intensity = sorted_bases[0][1]
        if max_intensity == 0:  # No signal at this position
            return None

        valid_bases = [(base, intensity) for base, intensity in sorted_bases
                      if intensity > max_intensity * self.min_peak_ratio]

        # Return top two bases if available
        if len(valid_bases) >= 2:
            print(f"doule peak detected, return {valid_bases[:2]}")
            return valid_bases[:2]
        elif len(valid_bases) == 1:
            print("single peak")
            return [valid_bases[0]]  # Single peak
        else:
            print("No peak")
            return None

    def detect_heterozygous_positions(self, ab1_data):
        """Detect positions with double peaks (heterozygous positions)."""
        trace_data = ab1_data['trace_data']
        print(trace_data)
        sequence = ab1_data['sequence']
        quality_scores = ab1_data['quality_scores']

        heterozygous_positions = {}
        base_calls = []

        for pos in range(len(sequence)):
            if pos >= len(quality_scores) or quality_scores[pos] < self.min_quality:
                base_calls.append({
                    'position': pos,
                    'primary_base': sequence[pos],
                    'secondary_base': None,
                    'is_heterozygous': False,
                    'quality': quality_scores[pos] if pos < len(quality_scores) else 0
                })
                continue

            peaks = self.analyze_peaks_at_position(trace_data, pos)

            if peaks and len(peaks) >= 2:
                primary_base, primary_intensity = peaks[0]
                secondary_base, secondary_intensity = peaks[1]

                # Only consider it heterozygous if secondary peak is significant
                ratio = secondary_intensity / primary_intensity
                if ratio > 0.3:  # Secondary peak is at least 30% of primary
                    heterozygous_positions[pos] = {
                        'primary_base': primary_base,
                        'secondary_base': secondary_base,
                        'primary_intensity': primary_intensity,
                        'secondary_intensity': secondary_intensity,
                        'quality': quality_scores[pos],
                        'peak_ratio': ratio
                    }

                    base_calls.append({
                        'position': pos,
                        'primary_base': primary_base,
                        'secondary_base': secondary_base,
                        'is_heterozygous': True,
                        'quality': quality_scores[pos],
                        'peak_ratio': ratio
                    })
                else:
                    base_calls.append({
                        'position': pos,
                        'primary_base': sequence[pos],
                        'secondary_base': None,
                        'is_heterozygous': False,
                        'quality': quality_scores[pos]
                    })
            else:
                base_calls.append({
                    'position': pos,
                    'primary_base': sequence[pos],
                    'secondary_base': None,
                    'is_heterozygous': False,
                    'quality': quality_scores[pos]
                })

        return heterozygous_positions, base_calls

    def generate_allele_sequences(self, base_calls, reference_seq):
        """Generate consensus and alternative sequences."""
        # Create primary sequence (consensus)
        primary_seq = ''.join([call['primary_base'] for call in base_calls])

        # For simple cases, we'll create alternative sequence by substituting heterozygous positions
        alternative_bases = []
        for call in base_calls:
            if call['is_heterozygous']:
                alternative_bases.append(call['secondary_base'])
            else:
                alternative_bases.append(call['primary_base'])

        alternative_seq = ''.join(alternative_bases)

        return primary_seq, alternative_seq

    def analyze(self, ab1_file, reference_input):
        """Main analysis function."""
        # Read and process reference sequence
        print(f"Reading reference sequence from: {reference_input}")
        reference_sequence = self.read_reference_sequence(reference_input)
        print(f"Reference sequence length: {len(reference_sequence)}")

        # Read AB1 file
        print(f"Reading AB1 file: {ab1_file}")
        ab1_data = self.read_ab1_file(ab1_file)
        print(f"AB1 sequence length: {len(ab1_data['sequence'])}")
        print(f"Quality scores length: {len(ab1_data['quality_scores'])}")

        # Detect heterozygous positions
        print("Analyzing peaks for heterozygous positions...")
        heterozygous_positions, base_calls = self.detect_heterozygous_positions(ab1_data)

        print(f"Found {len(heterozygous_positions)} heterozygous positions")

        # Generate allele sequences
        # 1. I need a aligner rewrite from known pairwise alignment.
        # first do ankering and then,


        print("Generating consensus and alternative sequences...")
        consensus_seq, alternative_seq = self.generate_allele_sequences(base_calls, reference_sequence)

        # Perform alignment for reporting
        alignment = self.aligner.align(reference_sequence, consensus_seq)[0]
        aligned_ref = str(alignment.target)
        aligned_consensus = str(alignment.query)

        # Prepare results
        results = {
            'input_file': ab1_file,
            'reference_source': reference_input,
            'reference_sequence': reference_sequence,
            'original_sequence': ab1_data['sequence'],
            'heterozygous_positions': heterozygous_positions,
            'base_calls': base_calls,
            'consensus_sequence': consensus_seq,
            'alternative_sequence': alternative_seq,
            'aligned_reference': aligned_ref,
            'aligned_consensus': aligned_consensus,
            'alignment_score': alignment.score,
            'alignment': alignment
        }

        return results

    def generate_report(self, results, output_file=None):
        """Generate a comprehensive analysis report."""
        report = []

        report.append("=" * 70)
        report.append("AB1 Double Peak Analysis Report")
        report.append("=" * 70)
        report.append(f"Input AB1 file: {results['input_file']}")
        report.append(f"Reference source: {results['reference_source']}")
        report.append(f"Reference length: {len(results['reference_sequence'])}")
        report.append(f"Original sequence length: {len(results['original_sequence'])}")
        report.append(f"Heterozygous positions found: {len(results['heterozygous_positions'])}")
        report.append(f"Alignment score: {results['alignment_score']:.2f}")
        report.append("")

        # Show heterozygous positions
        if results['heterozygous_positions']:
            report.append("Heterozygous Positions (Double Peaks):")
            report.append("-" * 60)
            report.append("Pos\tPrimary\tSecondary\tRatio\tQuality")
            report.append("-" * 60)
            for pos, info in sorted(results['heterozygous_positions'].items()):
                ratio = info['peak_ratio']
                report.append(f"{pos}\t{info['primary_base']}\t{info['secondary_base']}\t\t{ratio:.2f}\t{info['quality']}")
            report.append("")
        else:
            report.append("No heterozygous positions detected.")
            report.append("")

        # Show sequences
        report.append("Consensus Sequence (primary alleles):")
        report.append(results['consensus_sequence'])
        report.append("")

        report.append("Alternative Sequence (secondary alleles at heterozygous positions):")
        report.append(results['alternative_sequence'])
        report.append("")

        # Show alignment preview
        report.append("Alignment (Reference vs Consensus):")
        report.append("-" * 50)

        # Show first 100 characters of alignment
        ref_preview = results['aligned_reference'][:100]
        con_preview = results['aligned_consensus'][:100]

        report.append(f"Reference:  {ref_preview}")
        report.append(f"Consensus:  {con_preview}")
        report.append("")

        # Show differences between sequences
        report.append("Differences between consensus and alternative sequences:")
        report.append("-" * 55)

        consensus = results['consensus_sequence']
        alternative = results['alternative_sequence']

        min_len = min(len(consensus), len(alternative))
        differences = []
        for i in range(min_len):
            if consensus[i] != alternative[i]:
                differences.append((i, consensus[i], alternative[i]))

        if differences:
            report.append("Position\tConsensus\tAlternative")
            report.append("-" * 40)
            for pos, cons, alt in differences:
                report.append(f"{pos}\t\t{cons}\t\t{alt}")
        else:
            report.append("No differences between consensus and alternative sequences")

        report_text = '\n'.join(report)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_text)
            print(f"Report saved to: {output_file}")

        return report_text

def main():
    parser = argparse.ArgumentParser(
        description='Analyze double peaks in AB1 files for diploid sequencing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Reference sequence can be provided as:
- FASTA file (.fasta, .fa, .fna)
- GenBank file (.gb, .gbk)
- Text file (.txt) with sequence
- Direct sequence string

Examples:
  %(prog)s sample.ab1 reference.fasta -o report.txt
  %(prog)s sample.ab1 "ATCGATCGATCG" -o report.txt
  %(prog)s sample.ab1 reference.gb -q 25 -r 0.3
        '''
    )

    parser.add_argument('ab1_file', help='Input AB1 file path')
    parser.add_argument('reference', help='Reference sequence (file or direct sequence)')
    parser.add_argument('-o', '--output', help='Output report file')
    parser.add_argument('-q', '--min_quality', type=int, default=20,
                       help='Minimum quality score (default: 20)')
    parser.add_argument('-r', '--min_ratio', type=float, default=0.2,
                       help='Minimum peak ratio for secondary base (default: 0.2)')

    args = parser.parse_args()

    if not os.path.exists(args.ab1_file):
        print(f"Error: AB1 file '{args.ab1_file}' not found")
        sys.exit(1)

    # Initialize analyzer
    analyzer = AB1DoublePeakAnalyzer(
        min_quality=args.min_quality,
        min_peak_ratio=args.min_ratio
    )

    try:
        # Perform analysis
        results = analyzer.analyze(args.ab1_file, args.reference)

        # Generate report
        report = analyzer.generate_report(results, args.output)
        print(report)

    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
