import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation, CompoundLocation

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Read two files')

# Add the file arguments
parser.add_argument('-fasta', type=str, help='Path to the first file')
parser.add_argument('-gene_list', type=str, help='Path to the second file')
parser.add_argument("-folder", type=str, help="Path to the gff folder")
parser.add_argument("-extract_type", type=str, help="Path to the gff folder")
parser.add_argument("-suffix", type=str, help="suffix of your extract out put")

# Parse the command-line arguments
args = parser.parse_args()



def extend_region(mode="", gene_name="", number=[0,0]):
    ##if we use gene as a unit to do extension 
    if mode == "gene":
        left_flanking = number[0]
        right_flanking = number[1]
        gene_number = int(re.search(r"([0-9]+)", gene_name).group(1))
        left_gene = list(range(gene_number - left_flanking, gene_number))
        right_gene = list(range(gene_number, gene_number + right_flanking + 1)) 
        gene_list = left_gene + right_gene 
        new_gene_list = ['{:06d}'.format(num) for num in gene_list]
    else:
        gene_list = list(range(number[0],number[1]+1))
    return gene_list 



limit_info = dict(gff_id=["V"])
with open(args.gene_list, "r") as gene_list_file:
    for line in gene_list_file.readlines():
        line_elements = line.strip().split()
        strain = line_elements[0]
        gene_list = extend_region(mode="select",number=[int(line_elements[1]),int(line_elements[2])])
        if args.extract_type == "gff":
            with open(f"{args.folder}/caenorhabditis_{strain}.gff3") as gff, open(f"{args.folder}/caenorhabditis_{strain}.{args.suffix}.gff3","a") as out_gff:
                for line in gff.readlines(): 
                    for gene_single in gene_list:
                        if re.search(f"FUN_0+{gene_single}", line):
                            out_gff.write(line)
                            print(f"{strain}\t{line.strip()}")
            gff.close()
            out_gff.close()
        elif args.extract_type == "mRNA":
            pass
        elif args.extract_type == "genebank":
            start = int(line_elements[1]) - 1
            end = int(line_elements[2])
            locus = line_elements[3] 
            input_file = f"{args.folder}/caenorhabditis_{strain}.gbk"  # Replace with the path to your input GenBank file
            output_file = f"{args.folder}/caenorhabditis_{strain}.{args.suffix}.gbk"  # Replace with the path to the output GenBank file
            filtered_records = []

            # Iterate over each record in the input GenBank file
            for record in SeqIO.parse(input_file, "genbank"):
                print(record)
                if re.search(locus, record.id):
                    # Subset the sequence based on the desired region
                    subseq = record.seq[start:end]
                    
                    # Create a new SeqRecord with the subsequence and the same metadata as the original record
                    new_record = SeqRecord(subseq, id=record.id, name=record.name, description=record.description, dbxrefs=record.dbxrefs, annotations=record.annotations)

                    # Iterate over each feature in the record
                    for feature in record.features:
                        # Check if the feature is within the desired region
                        if start <= feature.location.start and feature.location.end <= end:
                            if isinstance(feature.location, FeatureLocation):
                                new_location = FeatureLocation(feature.location.start - start, feature.location.end - start, feature.location.strand)
                                new_feature = SeqFeature(location = new_location, type=feature.type, qualifiers=feature.qualifiers)
                            else:
                                new_location = [FeatureLocation(loc.start - start, loc.end - start, loc.strand) for loc in feature.location.parts]
                                new_location = CompoundLocation(new_location)
                                new_feature = SeqFeature(location = new_location, type=feature.type, qualifiers=feature.qualifiers)
                            # Add the new feature to the new record
                            new_record.features.append(new_feature)

                    #filtered_records.append(new_record)
                    break  # Stop searching for more features with the locus tag in the current record

            SeqIO.write(new_record, output_file, "genbank")
        elif args.extract_type == "protein":
            protein_dict = SeqIO.to_dict(SeqIO.parse(f"{args.folder}/caenorhabditis_{strain}.proteins.fa", "fasta"))
            protein_out =  open(f"{args.folder}/caenorhabditis_{strain}.{args.suffix}.proteins.fa","a")
            for gene_single in gene_list:
                for key in protein_dict.keys():
                    if re.search(str(gene_single), key):
                        protein_out.write(f">{key}_{strain}\n")
                        protein_out.write(str(protein_dict[key].seq)+"\n")
                        break
            protein_out.close()
    gene_list_file.close()


