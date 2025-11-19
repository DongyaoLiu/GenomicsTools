import argparse
import re

def transform_gene_data(input_file, output_file):
    """
    Transform gene data from the given format to a tab-separated table
    """
    # Read the input file
    with open(input_file, 'r') as f:
        content = f.readlines()

    for line in content:

        line = line.strip()
        if re.match("^#", line):
            continue

        # Split the content by lines
        if re.search("WBGene", line):

            elements = line.split()
            # Extract the main gene information from first line
            WBID = elements[0]
            gene_name = elements[1]
            protein_name = elements[2]
            current_section = "None"


        if line.startswith("Concise description:"):

            current_section = "concise"
            concise_desc = line.replace("Concise description: ", "")

        elif line.startswith("Automated description:"):

            current_section = "auto"
            automated_desc = line.replace("Automated description: ", "")

        elif line.startswith("Gene class description:"):


            current_section = "class"
            class_desc = line.replace("Gene class description: ", "")

        elif current_section == "concise":
            concise_desc += concise_desc + line

        elif current_section == "auto":
            automated_desc += automated_desc + line

        elif current_section == "class":
            class_desc = class_desc + line
            current_section = "None"

        if re.match("^=", line):
            concise_desc = ' '.join(concise_desc.split())
            automated_desc = ' '.join(automated_desc.split())
            class_desc = ' '.join(class_desc.split())


            with open(output_file, 'a') as fo:

                fo.write(f"{WBID}\t{gene_name}\t{protein_name}\t{class_desc}\n")


    print(f"Data successfully transformed and saved to {output_file}")


# read the data in and split by "="
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

    transform_gene_data(args.data, args.output)

    # Display the result
    print("\nGenerated tab-separated table:")
    print("=" * 20)
    with open(args.output, 'r') as f:
        print(f.read())
