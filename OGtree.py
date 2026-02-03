import re
from Bio import Phylo
from io import StringIO
import json
from collections import defaultdict

def find_closest_orthologs_in_tree(newick_string, og_name):
    """
    Processes a single Newick tree string to find, for each WB gene,
    the single closest gene from every other species.

    Args:
        newick_string: The Newick format tree string.
        og_name: Name of the orthogroup (for reporting).

    Returns:
        A dictionary: {WB_gene_id: [{species_name: closest_gene_id}, ...]}
    """
    # Parse the tree using Biopython
    try:
        tree = Phylo.read(StringIO(newick_string), "newick")
    except Exception as e:
        print(f"  Error parsing tree for {og_name}: {e}")
        return {}

    # Get all terminal nodes (leaves, i.e., genes)
    all_terminals = list(tree.get_terminals())

    # Extract gene info and find WB genes
    gene_data = []
    wb_terminals = []
    wb_pattern = r'WBGene\d+'

    for term in all_terminals:
        # term.name format is expected to be like "Species|GeneID"
        if term.name and '|' in term.name:
            species, gene_id = term.name.split('|', 1)
            # Clean any branch length suffix if present
            gene_id = re.split(r'[:,\s]', gene_id)[0]
        else:
            # Fallback if format differs
            species = "Unknown"
            gene_id = term.name if term.name else f"Unknown_{id(term)}"

        node_info = {
            'node': term,
            'species': species,
            'gene_id': gene_id,
            'full_name': term.name
        }
        gene_data.append(node_info)

        if re.search(wb_pattern, gene_id):
            wb_terminals.append(node_info)

    if not wb_terminals:
        return {}

    result_dict = {}
    # For each WB gene, calculate distance to all others
    for wb_info in wb_terminals:
        wb_gene_id = wb_info['gene_id']
        wb_species = wb_info['species']
        wb_node = wb_info['node']

        # Dictionary to store the closest gene found so far for each species
        # Format: {species: {'gene_id': id, 'distance': dist}}
        best_per_species = {}

        for other_info in gene_data:
            other_gene_id = other_info['gene_id']
            other_species = other_info['species']
            other_node = other_info['node']

            # Skip self-comparison and genes from the same species as WB
            if other_node is wb_node or other_species == wb_species:
                continue

            # Calculate phylogenetic distance using Biopython
            try:
                # This is the key function: calculates sum of branch lengths
                # between the two nodes (genes) on the tree.
                distance = tree.distance(wb_node, other_node)
            except Exception as e:
                # Fallback if distance calculation fails
                print(f"    Distance calc failed for {wb_gene_id} to {other_gene_id}: {e}")
                distance = float('inf')

            # Update if this gene is closer than the previous best for its species
            if (other_species not in best_per_species or
                    distance < best_per_species[other_species]['distance']):
                best_per_species[other_species] = {
                    'gene_id': other_gene_id,
                    'distance': distance
                }

        # Convert the result to the requested list-of-dictionaries format
        species_list = [{sp: data['gene_id']} for sp, data in best_per_species.items()]
        result_dict[wb_gene_id] = species_list

    return result_dict


def process_orthofinder_trees_file(input_file_path, output_file_path):
    """
    Main function to parse an OrthoFinder trees file, process each orthogroup tree,
    and save the results.

    Args:
        input_file_path: Path to the OrthoFinder 'Orthogroups_Trees.txt' file.
        output_file_path: Path where the JSON results will be saved.
    """
    with open(input_file_path, 'r') as f:
        file_content = f.read()

    # This pattern tries to capture: "OG0000001_tree" and the tree string after it.
    # OrthoFinder trees are in Newick format and end with a semicolon.
    og_tree_pattern = r'(OG\d+)_tree\s*\n([^;]+;)'
    matches = list(re.finditer(og_tree_pattern, file_content, re.MULTILINE | re.DOTALL))

    if not matches:
        print("Warning: No OG trees found with the standard pattern.")
        print("Attempting alternative parsing by splitting on 'OG' prefixes...")
        # Fallback: split the file content by 'OG' and try to parse each section
        sections = re.split(r'\n(OG\d+)', file_content)[1:]  # Skip text before first OG
        matches = []
        for i in range(0, len(sections), 2):
            if i+1 < len(sections):
                og_name = sections[i]
                # The tree part might be in the next section; join until a semicolon
                tree_section = sections[i+1]
                # Extract the first complete Newick tree (up to a semicolon)
                tree_match = re.search(r'([^;]+;)', tree_section)
                if tree_match:
                    matches.append((og_name, tree_match.group(1)))

    print(f"Found {len(matches)} orthogroup trees to process.")

    all_results = {}
    for match in matches:
        if isinstance(match, re.Match):
            og_name = match.group(1)
            tree_str = match.group(2).strip()
        else:
            og_name, tree_str = match
            tree_str = tree_str.strip()

        print(f"Processing {og_name}...")
        og_result = find_closest_orthologs_in_tree(tree_str, og_name)
        if og_result:
            all_results[og_name] = og_result

    # Save the final dictionary to a JSON file
    with open(output_file_path, 'w') as out_f:
        json.dump(all_results, out_f, indent=2)

    print(f"\nProcessing complete. Results saved to '{output_file_path}'.")

    # Print a brief summary
    total_wb_genes = sum(len(genes) for genes in all_results.values())
    total_relationships = sum(
        len(ortho_list)
        for og in all_results.values()
        for ortho_list in og.values()
    )
    print(f"  Total OG groups with results: {len(all_results)}")
    print(f"  Total WB genes processed: {total_wb_genes}")
    print(f"  Total closest-ortholog relationships found: {total_relationships}")

    return all_results


# Example usage
if __name__ == "__main__":
    # Update these paths to match your files
    input_file = "/lustre1/g/sbs_cgz/dongyao/C35_geneprediction/all_protein2/primary_transcripts/OrthoFinder/Results_Jan06/"   # Your OrthoFinder output
    output_file = "closest_orthologs_per_species.json"

    results = process_orthofinder_trees_file(input_file, output_file)

    # Optional: Print a small sample from the first OG group
    for og_name, og_data in list(results.items())[:1]:
        print(f"\nSample from {og_name}:")
        for wb_gene, ortho_list in list(og_data.items())[:2]:  # First 2 WB genes
            print(f"  {wb_gene}:")
            for ortho_dict in ortho_list[:3]:  # First 3 species
                print(f"    {ortho_dict}")
