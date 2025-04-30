import numpy as np
import ete3
import tarfile
import re
import sys
from Bio import AlignIO
from Bio.SeqIO.FastaIO import FastaWriter

data_dir = 'data'
tree_archive = f'{data_dir}/Tree.tar.gz'
msa_archive = f'{data_dir}/MSAs.tar.gz'

def phylogenetically_ordered_species(newick_string, species_list=None):
    """
    Returns species ordered by increasing evolutionary distance from queried species.
    """
    tree = ete3.Tree(newick_string, format=1)

    #if not species_list:
    query_species = "Homo_sapiens"
    # elif type(species_list) == list:
    #     query_species = str(species_list[0])
    # else:
    #     raise AttributeError("Species list is not a list.")

    # Find the Homo sapiens node
    query_node = None
    for leaf in tree.get_leaves():
        if leaf.name == query_species:
            query_node = leaf
            break
    if not query_node:
        raise KeyError(f'Leaf not found for {query_species}')
    
    # Get all nodes in the tree
    all_leaves = tree.get_leaves()
    
    # Calculate distance from each leaf to the human node
    distances = []
    for leaf in all_leaves:
        distance = tree.get_distance(query_node, leaf)
        distances.append((leaf.name, distance))
    
    # Sort by distance (ascending)
    distances.sort(key=lambda x: x[1])
    
    # Extract just the species names in order of distance
    ordered_species = [item[0].replace("_"," ") for item in distances]
    
    return ordered_species

# Get the A0A183 protein MSA
with open(tree_archive, 'rb') as f:
    with tarfile.open(fileobj=f, mode='r:gz') as tar:
        tree = tar.getmember('Tree/Q92870_pruned.tre')
        tree = tar.extractfile(tree)
        tree = tree.read().decode('utf-8')

def reorder_msa_by_human_distance(msa_file, tree_file, output_file, format="fasta"):
    """
    Reorder sequences in an MSA based on evolutionary distance from Homo sapiens.
    
    Parameters:
    - msa_file: Path to the MSA file
    - tree_file: Path to the Newick tree file
    - output_file: Path to write the reordered MSA
    - format: MSA file format (default: fasta)
    - count_species: if True, return a dict of species occurence
    """
    # Read the alignment
    alignment = AlignIO.read(msa_file)
    
    # Read the tree
    with open(tree_file, 'r') as f:
        newick_string = f.read().strip()
    
    # Extract species names from sequence IDs
    seq_to_species = {}
    species_in_msa = set()
    pattern=r"OS=([^=]+?)\sOX="
    for record in alignment:
        species = re.search(pattern, record.description).group(1)
        seq_to_species[record.id] = species
        species_in_msa.add(species)
    
    # Get species ordered by distance from humans
    ordered_species = phylogenetically_ordered_species(newick_string, list(species_in_msa))
    
    # Create a dictionary to store sequences by species
    species_to_records = {}
    for record in alignment:
        species = seq_to_species[record.id]
        if species not in species_to_records:
            species_to_records[species] = []
        species_to_records[species].append(record)
    
    # Create a new alignment with sequences ordered by distance from humans
    new_records = []
    for species in ordered_species:
        if species in species_to_records:
            new_records.extend(species_to_records[species])
    
    # Write the reordered alignment
    alignment = AlignIO.MultipleSeqAlignment(new_records)

    with open(output_file, "w") as handle:
        writer = FastaWriter(handle, wrap=None)
        writer.write_header()
        for record in alignment:
            writer.write_record(record)
        writer.write_footer()
    
    # print(f"Reordered alignment written to {output_file}")
    # print("Species ordered by distance from Homo sapiens:")
    # for species in ordered_species:
    #     print(f"- {species}")
    
    return new_records

def count_species(alignment_dir:str, proteins_list:list):
    """Count species occurences in MSA files"""
    pattern_id = r"sp\|([^|]+)\|"
    pattern=r"OS=([^=]+?)\sOX="
    species_count = {}

    # Read MSA files
    for protein in proteins_list:
        with open(f"{alignment_dir}/{protein}") as msa_file:
            alignment = AlignIO.read(msa_file, "fasta")

        # Get query protein id        
        id_protein = re.search(pattern_id, alignment[0].id).group(1)

        # Iterate over record in MSA and extract species occurences
        species_count[id_protein] = {}
        for record in alignment[1:]:
            species = re.search(pattern, record.description).group(1)
            if species in species_count[id_protein]:
                species_count[id_protein][species] += 1
            else:
                species_count[id_protein][species] = 1
    
    return species_count

def deduplicate_msa(input_file, output_file):
    """
    Deduplicate species in an MSA file by keeping first occurence of species"""

    # Read the alignment
    alignment = AlignIO.read(input_file, "fasta")
    
    # Create a set to track seen species
    seen_species = set()
    
    # Create a list to store unique records
    unique_records = []
    
    # Iterate over the records in the alignment
    for record in alignment:
        # Extract the species name from the record description
        species = re.search(r"OS=([^=]+?)\sOX=", record.description).group(1)
        
        # If the species has not been seen before, add it to the unique records
        if species not in seen_species:
            seen_species.add(species)
            unique_records.append(record)
    
    # create new alignment
    unique_alignment = AlignIO.MultipleSeqAlignment(unique_records)

    # Write the unique alignment to a new file
    with open(output_file, "w") as handle:
        writer = FastaWriter(handle, wrap=None)
        writer.write_header()
        for record in unique_alignment:
            writer.write_record(record)
        writer.write_footer()