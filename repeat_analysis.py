import pandas as pd
from collections import Counter
from itertools import combinations, permutations
import argparse


def parse_bed_file(bed_file):
    """Parse BED file and return a sorted DataFrame"""
    df = pd.read_csv(bed_file, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    # Sort by chromosome and start position
    df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    return df


def get_combinations_in_window(window_df, order_matters=True):
    """Get all possible combinations of repeats within a window"""
    elements = window_df['name'].tolist()
    
    all_combinations = []
    for r in range(1, len(elements) + 1):
        if order_matters:
            # Generate all possible ordered combinations (permutations)
            for combo in permutations(elements, r):
                all_combinations.append(tuple(combo))
        else:
            # Generate all possible unordered combinations
            for combo in combinations(elements, r):
                all_combinations.append(tuple(sorted(combo)))
    
    return all_combinations


def distance_filtering(window_df, max_window_span, max_adjacent_distance):
    """Apply distance filtering based on window span and adjacent distances"""
    if len(window_df) == 0:
        return True
    
    # Calculate window span (distance from smallest start to largest end)
    min_start = window_df['start'].min()
    max_end = window_df['end'].max()
    window_span = max_end - min_start
    
    if max_window_span and window_span > max_window_span:
        return False
    
    # Check adjacent distances
    if max_adjacent_distance and len(window_df) > 1:
        sorted_df = window_df.sort_values(['start'])
        for i in range(len(sorted_df) - 1):
            current_end = sorted_df.iloc[i]['end']
            next_start = sorted_df.iloc[i + 1]['start']
            distance = next_start - current_end
            if distance > max_adjacent_distance:
                return False
    
    return True


def analyze_repeats_in_windows(bed_file, window_size, order_matters=True, 
                              max_window_span=None, max_adjacent_distance=None):
    """Analyze repeats in sliding windows with distance filtering"""
    df = parse_bed_file(bed_file)
    
    all_combinations = Counter()
    
    for start_idx in range(len(df)):
        # Get the window of specified size
        end_idx = min(start_idx + window_size, len(df))
        window_df = df.iloc[start_idx:end_idx]
        
        # Apply distance filtering
        if not distance_filtering(window_df, max_window_span, max_adjacent_distance):
            continue
        
        # Get combinations in this window
        window_combinations = get_combinations_in_window(
            window_df, order_matters
        )
        
        # Update counter
        for combo in window_combinations:
            all_combinations[combo] += 1
    
    # Find most frequent combinations
    most_frequent = all_combinations.most_common(10)  # Top 10 most frequent
    
    return most_frequent, all_combinations


def main():
    parser = argparse.ArgumentParser(description='Analyze repeats in BED file with window-based approach')
    parser.add_argument('bed_file', help='Input BED file')
    parser.add_argument('--window-size', type=int, default=10, help='Window size (number of lines)')
    parser.add_argument('--window-span', type=int, help='Max window span (distance from min start to max end)')
    parser.add_argument('--adjacent-distance', type=int, help='Max distance between adjacent repeats')
    parser.add_argument('--order-matters', action='store_true', help='Consider order of repeats')
    
    args = parser.parse_args()
    
    print(f"Analyzing repeats in {args.bed_file}")
    print(f"Window size: {args.window_size}")
    print(f"Order matters: {args.order_matters}")
    if args.window_span:
        print(f"Max window span: {args.window_span}")
    if args.adjacent_distance:
        print(f"Max adjacent distance: {args.adjacent_distance}")
    
    # Analysis with order matters
    print("\n=== Analysis with Order Matters ===")
    most_freq_ordered, counter_ordered = analyze_repeats_in_windows(
        args.bed_file, args.window_size, order_matters=True,
        max_window_span=args.window_span, max_adjacent_distance=args.adjacent_distance
    )
    
    print("Most frequent combinations (order matters):")
    for combo, count in most_freq_ordered:
        print(f"  {combo}: {count} times")
    
    # Analysis with order doesn't matter
    print("\n=== Analysis with Order Doesn't Matter ===")
    most_freq_unordered, counter_unordered = analyze_repeats_in_windows(
        args.bed_file, args.window_size, order_matters=False,
        max_window_span=args.window_span, max_adjacent_distance=args.adjacent_distance
    )
    
    print("Most frequent combinations (order doesn't matter):")
    for combo, count in most_freq_unordered:
        print(f"  {combo}: {count} times")


if __name__ == "__main__":
    main()