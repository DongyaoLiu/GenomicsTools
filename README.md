# Repeats Analysis Tools for P. pacificus

## Overview
This repository contains tools for analyzing clustering patterns and repeat combinations in P. pacificus genomes. Two main tools are available:
1. `RepeatsClusterAna.py` - analyzes clustering patterns of repetitive sequences
2. `repeat_analysis.py` - analyzes repeat combinations within windows with distance filtering

## Purpose
The scripts address the question of whether certain repeat sequences tend to co-exist in close genomic proximity, potentially indicating functional relationships such as centromere components.

## Input Format
The scripts take a BED file with at least 4 columns:
```
Chr    Start    End      RepeatID    [Score]    [Strand]
chr1   1000     1200     RepeatA     0          +
chr1   1100     1300     RepeatB     0          -
...
```

## Tools

### 1. Repeats Cluster Analysis (`RepeatsClusterAna.py`)

#### Output Format
The output is a BED file with two additional columns:
- `cluster_id`: Identifier for the genomic cluster
- `p_value`: Statistical significance of the co-occurrence pattern

Format:
```
Chr    Start    End      RepeatID    ClusterID    P-value
```

#### Methodology
- **Clustering**: Groups repeats that are within a specified distance (`--max_distance`, default 10,000 bp) of each other
- **Statistical Analysis**: Calculates the expected diversity of repeat types in each cluster based on genome-wide frequencies
- **Significance Calculation**: Computes a p-value-like measure for enrichment of diverse repeat types in the cluster

#### Usage
```bash
python RepeatsClusterAna.py input.bed --max_distance 10000 --output output.bed
```

### 2. Window-based Repeat Analysis (`repeat_analysis.py`)

#### Features
- **Window-based analysis**: Analyze combinations of repeats within defined windows
- **Order consideration**: Analyze with or without considering the order of repeats
- **Distance filtering**: Apply two types of distance filtering:
  - Window span: Maximum distance between the smallest start and largest end in a window
  - Adjacent distance: Maximum distance between adjacent repeats
- **Most frequent combinations**: Identify the most common repeat combinations

#### Usage
```bash
python repeat_analysis.py [BED_FILE] [OPTIONS]
```

##### Options
- `--window-size N`: Size of the window (number of lines to consider) [default: 10]
- `--window-span N`: Maximum window span (distance from min start to max end)
- `--adjacent-distance N`: Maximum distance between adjacent repeats
- `--order-matters`: Consider order of repeats (default: order doesn't matter)

##### Examples
Basic analysis with window size of 3:
```bash
python repeat_analysis.py sample.bed --window-size 3
```

Analysis with distance filtering:
```bash
python repeat_analysis.py sample.bed --window-size 3 --window-span 2000 --adjacent-distance 1000
```

Analysis considering order of repeats:
```bash
python repeat_analysis.py sample.bed --window-size 3 --order-matters
```

#### Output
The tool outputs two analyses:
1. Analysis where order of repeats matters (permutations)
2. Analysis where order doesn't matter (combinations)

For each analysis, it shows the top 10 most frequent repeat combinations with their occurrence counts.

## Interpretation
- **Low p-values** (e.g., < 0.05): Indicate clusters where different repeat types co-occur more frequently than expected by chance, suggesting potential functional relationships
- **High p-values**: Indicate clusters with repeat compositions consistent with random expectation
- **Cluster IDs**: Allow identification of specific genomic intervals with significant co-occurrence patterns
- **Most frequent combinations**: Help identify common repeat arrangements in genomic windows

## Application to Centromere Analysis
These analyses can help identify potential centromere components by detecting clusters of repetitive elements that consistently co-occur, which may represent the complex repeat structures characteristic of centromeres.