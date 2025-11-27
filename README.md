# Repeats Cluster Analysis for P. pacificus

## Overview
`RepeatsClusterAna.py` is a Python script designed to analyze clustering patterns of repetitive sequences in P. pacificus genomes. The script identifies genomic intervals where different repeat types co-occur more frequently than expected by chance, which may indicate centromere components or other functionally related repetitive elements.

## Purpose
The script addresses the question of whether certain repeat sequences tend to co-exist in close genomic proximity, potentially indicating functional relationships such as centromere components.

## Input Format
The script takes a BED file with the following format:
```
Chr    Start    End      RepeatID
chr1   1000     1200     RepeatA
chr1   1100     1300     RepeatB
...
```

## Output Format
The output is a BED file with two additional columns:
- `cluster_id`: Identifier for the genomic cluster
- `p_value`: Statistical significance of the co-occurrence pattern

Format:
```
Chr    Start    End      RepeatID    ClusterID    P-value
```

## Methodology

### 1. Clustering
- Groups repeats that are within a specified distance (`--max_distance`, default 10,000 bp) of each other
- Each cluster represents a genomic interval with potentially co-occurring repeats

### 2. Statistical Analysis
- Calculates the expected diversity of repeat types in each cluster based on genome-wide frequencies
- Compares observed diversity to expected diversity
- Computes a p-value-like measure for enrichment of diverse repeat types in the cluster
- Uses probability calculations to determine if co-occurrence is more frequent than expected by chance

### 3. Significance Calculation
For each cluster:
- `observed_diversity` = number of unique repeat types in the cluster
- `expected_diversity` = sum of probabilities of seeing each repeat type in a random cluster of the same size
- If `observed_diversity > expected_diversity`, calculates enrichment ratio and returns a lower p-value

## Usage

```bash
python RepeatsClusterAna.py input.bed --max_distance 10000 --output output.bed
```

### Arguments:
- `bed_file`: Input BED file with repeat coordinates
- `--max_distance`: Maximum distance between repeats to be considered in the same cluster (default: 10000)
- `--output`: Output file name (default: repeats_cluster_analysis.bed)

## Interpretation
- **Low p-values** (e.g., < 0.05): Indicate clusters where different repeat types co-occur more frequently than expected by chance, suggesting potential functional relationships
- **High p-values**: Indicate clusters with repeat compositions consistent with random expectation
- **Cluster IDs**: Allow identification of specific genomic intervals with significant co-occurrence patterns

## Application to Centromere Analysis
This analysis can help identify potential centromere components by detecting clusters of repetitive elements that consistently co-occur, which may represent the complex repeat structures characteristic of centromeres.