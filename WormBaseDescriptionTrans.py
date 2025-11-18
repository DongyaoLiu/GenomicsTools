import argparse
import re

def transform_gene_data(input_file, output_file):
    """
    Transform gene data from the given format to a tab-separated table
    """
    # Read the input file
    with open(input_file, 'r') as f:
        content = f.readlines().strip()

    # Split the content by lines
    if re.search("WBGene", content):

        elements = content.split()
        # Extract the main gene information from first line
        gene_info = lines[0].split('\t')
        gene_id = gene_info[0]
        gene_name = gene_info[1]
        sequence_id = gene_info[2]

    # Initialize variables for other sections
    concise_desc = ""
    automated_desc = ""
    gene_class_desc = ""

    current_section = ""

    # Process each line to extract different sections
    for line in lines[1:]:
        line = line.strip()
        if line.startswith("Concise description:"):
            current_section = "concise"
            concise_desc = line.replace("Concise description:", "").strip()
        elif line.startswith("Automated description:"):
            current_section = "automated"
            automated_desc = line.replace("Automated description:", "").strip()
        elif line.startswith("Gene class description:"):
            current_section = "gene_class"
            gene_class_desc = line.replace("Gene class description:", "").strip()
        elif line and current_section == "concise":
            concise_desc += " " + line
        elif line and current_section == "automated":
            automated_desc += " " + line
        elif line and current_section == "gene_class":
            gene_class_desc += " " + line

    # Clean up the descriptions (remove extra spaces)
    concise_desc = ' '.join(concise_desc.split())
    automated_desc = ' '.join(automated_desc.split())
    gene_class_desc = ' '.join(gene_class_desc.split())

    # Write to tab-separated table
    with open(output_file, 'w') as f:
        # Write header
        f.write("Gene_ID\tGene_Name\tSequence_ID\tConcise_Description\tAutomated_Description\tGene_Class_Description\n")

        # Write data
        f.write(f"{gene_id}\t{gene_name}\t{sequence_id}\t{concise_desc}\t{automated_desc}\t{gene_class_desc}\n")

    print(f"Data successfully transformed and saved to {output_file}")

# Alternative function if you want to work with the data as a string directly
def transform_gene_data_from_string(data_string, output_file):
    """
    Transform gene data from string to tab-separated table
    """
    # Split the content by lines
    lines = data_string.strip().split('\n')

    # Extract the main gene information from first line
    gene_info = lines[0].split('\t')
    gene_id = gene_info[0]
    gene_name = gene_info[1]
    sequence_id = gene_info[2]

    # Initialize variables for other sections
    concise_desc = ""
    automated_desc = ""
    gene_class_desc = ""

    current_section = ""

    # Process each line to extract different sections
    for line in lines[1:]:
        line = line.strip()
        if line.startswith("Concise description:"):
            current_section = "concise"
            concise_desc = line.replace("Concise description:", "").strip()
        elif line.startswith("Automated description:"):
            current_section = "automated"
            automated_desc = line.replace("Automated description:", "").strip()
        elif line.startswith("Gene class description:"):
            current_section = "gene_class"
            gene_class_desc = line.replace("Gene class description:", "").strip()
        elif line and current_section == "concise":
            concise_desc += " " + line
        elif line and current_section == "automated":
            automated_desc += " " + line
        elif line and current_section == "gene_class":
            gene_class_desc += " " + line

    # Clean up the descriptions (remove extra spaces)
    concise_desc = ' '.join(concise_desc.split())
    automated_desc = ' '.join(automated_desc.split())
    gene_class_desc = ' '.join(gene_class_desc.split())

    # Write to tab-separated table
    with open(output_file, 'w') as f:
        # Write header
        f.write("Gene_ID\tGene_Name\tSequence_ID\tConcise_Description\tAutomated_Description\tGene_Class_Description\n")

        # Write data
        f.write(f"{gene_id}\t{gene_name}\t{sequence_id}\t{concise_desc}\t{automated_desc}\t{gene_class_desc}\n")

    print(f"Data successfully transformed and saved to {output_file}")

# Example usage:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='WormBase geneDescription transform',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''empty'''
     )
    parser.add_argument('data', help='Input AB1 file path')
    parser.add_argument('output', type=str,
                        help='specify the output file')
    args = parser.parse_args()

    transform_gene_data_from_string(args.data, args.output)

    # Display the result
    print("\nGenerated tab-separated table:")
    print("=" * 100)
    with open("gene_data.tsv", 'r') as f:
        print(f.read())
