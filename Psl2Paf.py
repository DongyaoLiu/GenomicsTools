#!/usr/bin/env python3
import argparse
import os

def psl_to_blocks(input_psl, output_tsv):
    """
    Convert PSL to block-wise 7-column format (one line per block):
    qName, qStart, qEnd, qStrand, tName, tStart, tEnd, tStrand
    - Preserves all alignment blocks
    - Correctly handles reverse-complement coordinates for '-' strand
    """
    with open(input_psl) as fin, open(output_tsv, 'w') as fout:

        # Write header

        fout.write("qName\tqStart\tqEnd\tqStrand\ttName\ttStart\ttEnd\ttStrand\n")

        for line in fin:
            if line.startswith('#'):
                continue  # Skip header lines

            fields = line.strip().split('\t')
            if len(fields) < 21:
                continue  # Skip malformed lines

            # Extract basic fields
            q_name = fields[9]
            strand = fields[8]
            q_strand = strand[0]
            t_strand = strand[1]
            q_size = int(fields[10])

            t_name = fields[13]
            t_size = int(fields[14])

            # Get all blocks
            block_sizes = [int(x) for x in fields[18].split(',') if x]
            q_starts = [int(x) for x in fields[19].split(',') if x]
            t_starts = [int(x) for x in fields[20].split(',') if x]

            # Process each block
            for i in range(len(block_sizes)):
                # Query coordinates
                if q_strand == '+':
                    q_start = q_starts[i]
                    q_end = q_starts[i] + block_sizes[i]
                else:
                    # Reverse complement coordinates for minus strand
                    q_start = q_size - (q_starts[i] + block_sizes[i])
                    q_end = q_size - q_starts[i]

                # Target coordinates (always in forward orientation)
                if t_strand == '+':
                    t_start = t_starts[i]
                    t_end = t_starts[i] + block_sizes[i]
                else:
                    t_start = t_size - (t_starts[i] + block_sizes[i])
                    t_end = t_size - t_starts[i]
                # Write output line for this block
                fout.write(f"{q_name}\t{q_start}\t{q_end}\t{q_strand}\t")
                fout.write(f"{t_name}\t{t_start}\t{t_end}\t{t_strand}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Convert PSL to block-wise 7-column format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, help="Input PSL file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()
    if os.path.getsize(args.input) != 0:
        psl_to_blocks(args.input, args.output)
        print(f"Conversion complete. Output written to {args.output}")
    else:
        print(f"{args.input} is empty!")

if __name__ == "__main__":
    main()
