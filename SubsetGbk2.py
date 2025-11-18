#!/usr/bin/env python3
import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation


class GenBankExtractor:
    """A class to handle extraction of specific regions from GenBank files."""

    def __init__(self, gene_list_path, input_folder, output_suffix):
        """
        Initialize the GenBankExtractor.

        Args:
            gene_list_path (str): Path to the gene list file.
            input_folder (str): Path to the folder containing GenBank files.
            output_suffix (str): Suffix for the output file.
        """
        self.gene_list_path = gene_list_path
        self.input_folder = input_folder
        self.output_suffix = output_suffix

    def parse_gene_list(self):
        """
        Parse the gene list file.

        Returns:
            list: A list of tuples containing (strain, start, end, chromosome).
        """
        gene_list = []
        with open(self.gene_list_path, "r") as file:
            for line in file:
                strain, start, end, chromosome = line.strip().split()
                gene_list.append((strain, int(start), int(end), chromosome))
        return gene_list

    def extract_region(self, strain, start, end, chromosome):
        """
        Extract a specific region from a GenBank file.

        Args:
            strain (str): Strain identifier.
            start (int): Start position of the region.
            end (int): End position of the region.
            chromosome (str): Chromosome identifier.

        Returns:
            SeqRecord: A SeqRecord object containing the extracted region.
        """
        input_file = f"{self.input_folder}/{strain}.gbk"
        output_file = f"{self.input_folder}/{strain}.{self.output_suffix}.gbk"

        for record in SeqIO.parse(input_file, "genbank"):
            if re.search(chromosome, record.id):
                # Adjust start to be 0-based for Python indexing
                start -= 1

                # Subset the sequence
                subseq = record.seq[start:end]

                # Create a new SeqRecord
                new_record = SeqRecord(
                    subseq,
                    id=record.id,
                    name=record.name,
                    description=record.description,
                    dbxrefs=record.dbxrefs,
                    annotations=record.annotations,
                )

                # Adjust and copy features
                for feature in record.features:
                    if start <= feature.location.start and feature.location.end <= end:
                        if isinstance(feature.location, FeatureLocation):
                            new_location = FeatureLocation(
                                feature.location.start - start,
                                feature.location.end - start,
                                feature.location.strand,
                            )
                        else:
                            new_location = CompoundLocation(
                                [
                                    FeatureLocation(loc.start - start, loc.end - start, loc.strand)
                                    for loc in feature.location.parts
                                ]
                            )
                        new_feature = SeqFeature(
                            location=new_location,
                            type=feature.type,
                            qualifiers=feature.qualifiers,
                        )
                        new_record.features.append(new_feature)

                # Save the extracted region to a new file
                SeqIO.write(new_record, output_file, "genbank")
                return new_record
        return None

    def process_all_regions(self):
        """Process all regions specified in the gene list."""
        gene_list = self.parse_gene_list()
        for strain, start, end, chromosome in gene_list:
            self.extract_region(strain, start, end, chromosome)
            print(f"Extracted region for {strain} (chromosome {chromosome})")


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract specific regions from GenBank files based on gene coordinates."
    )
    parser.add_argument(
        "-gene_list", type=str, required=True, help="Path to the gene list file (strain, start, end, chromosome)"
    )
    parser.add_argument(
        "-folder", type=str, required=True, help="Path to the folder containing GenBank files"
    )
    parser.add_argument(
        "-suffix", type=str, required=True, help="Suffix for the output file"
    )
    return parser.parse_args()


def main():
    """Main function to handle the extraction process."""
    args = parse_arguments()

    # Initialize the GenBankExtractor
    extractor = GenBankExtractor(args.gene_list, args.folder, args.suffix)

    # Process all regions in the gene list
    extractor.process_all_regions()


if __name__ == "__main__":
    main()
