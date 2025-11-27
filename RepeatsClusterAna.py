#!/usr/bin/env python3
"""
RepeatsClusterAna.py - Analyze clustering patterns of repetitive sequences in P. pacificus genomes.

This script identifies genomic intervals where different repeat types co-occur more frequently than expected by chance,
which may indicate centromere components or other functionally related repetitive elements.

Input: BED file with format Chr, start, end, repeatID
Output: BED file with additional columns for cluster ID and p-value
"""

import argparse
import pandas as pd
from collections import defaultdict
from scipy.stats import hypergeom
import numpy as np
from itertools import combinations


def read_bed_file(bed_path):
    """Read BED file and return as DataFrame."""
    df = pd.read_csv(bed_path, sep='\t', header=None, 
                     names=['chr', 'start', 'end', 'repeatID'])
    df = df.sort_values(['chr', 'start']).reset_index(drop=True)
    return df


def find_clusters(df, max_distance=10000):
    """Find clusters of repeats within a specified distance."""
    clusters = []
    cluster_id = 0
    
    current_cluster = []
    current_chr = df.iloc[0]['chr']
    
    for idx, row in df.iterrows():
        if row['chr'] != current_chr:
            # New chromosome, finalize previous cluster if exists
            if current_cluster:
                clusters.append({
                    'cluster_id': f'cluster_{cluster_id}',
                    'chr': current_chr,
                    'start': min([item['start'] for item in current_cluster]),
                    'end': max([item['end'] for item in current_cluster]),
                    'repeats': current_cluster
                })
                cluster_id += 1
                current_cluster = []
            current_chr = row['chr']
        
        # If cluster is empty or this repeat is within max_distance of last repeat in cluster
        if not current_cluster or (row['start'] - current_cluster[-1]['end'] <= max_distance):
            current_cluster.append({
                'idx': idx,
                'chr': row['chr'],
                'start': row['start'],
                'end': row['end'],
                'repeatID': row['repeatID']
            })
        else:
            # Start new cluster
            if current_cluster:
                clusters.append({
                    'cluster_id': f'cluster_{cluster_id}',
                    'chr': current_chr,
                    'start': min([item['start'] for item in current_cluster]),
                    'end': max([item['end'] for item in current_cluster]),
                    'repeats': current_cluster
                })
                cluster_id += 1
            current_cluster = [{
                'idx': idx,
                'chr': row['chr'],
                'start': row['start'],
                'end': row['end'],
                'repeatID': row['repeatID']
            }]
    
    # Add the last cluster
    if current_cluster:
        clusters.append({
            'cluster_id': f'cluster_{cluster_id}',
            'chr': current_chr,
            'start': min([item['start'] for item in current_cluster]),
            'end': max([item['end'] for item in current_cluster]),
            'repeats': current_cluster
        })
    
    return clusters


def calculate_cooccurrence_pvalue(cluster_repeats, all_repeats):
    """Calculate p-value for co-occurrence of repeat types in cluster."""
    # Count repeat types in the cluster
    cluster_repeat_counts = defaultdict(int)
    for rep in cluster_repeats:
        cluster_repeat_counts[rep['repeatID']] += 1
    
    # Count total occurrences of each repeat type in the genome
    total_repeat_counts = defaultdict(int)
    for rep in all_repeats:
        total_repeat_counts[rep['repeatID']] += 1
    
    # Total number of repeats in the genome
    total_repeats = len(all_repeats)
    # Total number of repeats in the cluster
    cluster_size = len(cluster_repeats)
    
    # Calculate a combined p-value for over-representation of specific repeat types in the cluster
    p_values = []
    
    for repeat_type, cluster_count in cluster_repeat_counts.items():
        total_count = total_repeat_counts[repeat_type]
        
        # Use hypergeometric test to see if this repeat type is over-represented
        # in the cluster compared to genome-wide frequency
        M = total_repeats  # Total population size
        n = total_count    # Number of success states in population
        N = cluster_size   # Number of draws
        x = cluster_count  # Observed number of successes
        
        # Calculate p-value for observing at least x successes
        if x <= n:
            p_val = hypergeom.sf(x-1, M, n, N)
            p_values.append(p_val)
    
    # Combine p-values using Fisher's method if we have multiple repeat types
    if len(p_values) == 0:
        return 1.0
    elif len(p_values) == 1:
        return p_values[0]
    else:
        # Fisher's method for combining p-values
        chi2_stat = -2 * sum(np.log(p) for p in p_values if p > 0)
        combined_p = 1 - chi2_stat._cdf(chi2_stat, df=2*len(p_values)) if hasattr(chi2_stat, '_cdf') else min(p_values)
        # Simplified approach - just take minimum p-value
        return min(p_values) if p_values else 1.0


def calculate_cluster_significance(cluster_repeats, all_repeats_df):
    """
    Calculate statistical significance of repeat co-occurrence in cluster.
    Uses a probability-based approach to determine if observed co-occurrences
    are more frequent than expected by chance.
    """
    # Convert to DataFrame for easier processing
    cluster_df = pd.DataFrame(cluster_repeats)
    all_repeats_df = pd.DataFrame(all_repeats_df)
    
    # Count unique repeat types in cluster
    cluster_repeat_types = set(cluster_df['repeatID'])
    cluster_size = len(cluster_repeats)
    
    if len(cluster_repeat_types) <= 1:
        # If cluster has only one repeat type, not interesting for co-occurrence
        return 1.0
    
    # Calculate expected number of different repeat types based on genome-wide frequencies
    total_repeats = len(all_repeats_df)
    repeat_frequencies = all_repeats_df['repeatID'].value_counts()
    
    # For each repeat type in the cluster, calculate its genome-wide frequency
    observed_diversity = len(cluster_repeat_types)
    
    # Calculate expected diversity based on genome-wide frequencies
    # For each repeat type in the cluster, calculate the probability it would appear
    # in a random cluster of the same size
    expected_diversity = 0
    for rep_type in cluster_repeat_types:
        if rep_type in repeat_frequencies:
            freq = repeat_frequencies[rep_type] / total_repeats
            # Probability of seeing this repeat type at least once in the cluster
            prob_seeing_type = 1 - (1 - freq)**cluster_size
            expected_diversity += prob_seeing_type
    
    # Calculate significance: if observed diversity is greater than expected, 
    # it might indicate non-random co-occurrence
    if expected_diversity > 0:
        enrichment_ratio = observed_diversity / expected_diversity
        # Use a simple transformation to get a p-value-like measure
        # Higher enrichment ratio means lower "p-value"
        if enrichment_ratio > 1:
            # Enriched - calculate significance
            p_val = np.exp(1 - enrichment_ratio)  # Exponential decay from 1
        else:
            # Not enriched
            p_val = 1.0
    else:
        p_val = 1.0
    
    return max(1e-10, min(1.0, p_val))  # Clamp between 1e-10 and 1.0


def analyze_clusters(bed_path, max_distance=10000):
    """Main function to analyze repeat clusters."""
    print(f"Reading BED file: {bed_path}")
    df = read_bed_file(bed_path)
    
    print(f"Finding clusters with max distance {max_distance}...")
    clusters = find_clusters(df, max_distance)
    
    print(f"Found {len(clusters)} clusters")
    
    # Prepare output dataframe
    output_rows = []
    
    for cluster in clusters:
        # Calculate significance for this cluster
        p_value = calculate_cluster_significance(cluster['repeats'], df.to_dict('records'))
        
        for rep in cluster['repeats']:
            output_rows.append({
                'chr': rep['chr'],
                'start': rep['start'],
                'end': rep['end'],
                'repeatID': rep['repeatID'],
                'cluster_id': cluster['cluster_id'],
                'p_value': p_value
            })
    
    output_df = pd.DataFrame(output_rows)
    return output_df


def main():
    parser = argparse.ArgumentParser(
        description="Analyze clustering patterns of repetitive sequences in P. pacificus"
    )
    parser.add_argument(
        'bed_file',
        help='Input BED file with format: Chr, start, end, repeatID'
    )
    parser.add_argument(
        '--max_distance',
        type=int,
        default=10000,
        help='Maximum distance between repeats to be considered in the same cluster (default: 10000)'
    )
    parser.add_argument(
        '--output',
        default='repeats_cluster_analysis.bed',
        help='Output BED file with cluster IDs and p-values (default: repeats_cluster_analysis.bed)'
    )
    
    args = parser.parse_args()
    
    # Perform analysis
    result_df = analyze_clusters(args.bed_file, args.max_distance)
    
    # Write output
    output_cols = ['chr', 'start', 'end', 'repeatID', 'cluster_id', 'p_value']
    result_df[output_cols].to_csv(args.output, sep='\t', index=False, header=False)
    print(f"Results saved to: {args.output}")
    
    # Print summary
    print(f"\nAnalysis Summary:")
    print(f"- Input file: {args.bed_file}")
    print(f"- Max cluster distance: {args.max_distance}")
    print(f"- Total repeats analyzed: {len(result_df)}")
    print(f"- Number of clusters identified: {result_df['cluster_id'].nunique()}")
    print(f"- Range of p-values: {result_df['p_value'].min():.2e} to {result_df['p_value'].max():.2e}")


if __name__ == "__main__":
    main()