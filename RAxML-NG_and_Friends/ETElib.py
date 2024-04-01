# Author: Jason Charamis
# Collection of functions to manipulate and visualize phylogenetic trees using the ETE3 toolkit #

#!/usr/bin/env python3

import seaborn as sns
import ete3
import re, random, os
import argparse

def parse_arguments():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating gff3 files.')
    parser.add_argument('-aln','--alignment', type=str, help='Input alignment.')
    parser.add_argument('-t','--tree', type=str, help='Nwk tree input file.')
    parser.add_argument('-phy','--phy', action="store_true", help='Convert input alignment to phy format.')
    parser.add_argument('-m','--midpoint', action="store_true",help='Midpoint root tree.')
    parser.add_argument('--collapse', action="store_true",help='Collapse nodes based on bootstrap support.')
    parser.add_argument('--cutoff', type=str,help='Bootstrap support cutoff for collapsing nodes.')
    parser.add_argument('-r','--resolve', action="store_true",help='Resolve polytomies.')
    parser.add_argument('-v','--visualize', action="store_true",help='Visualize tree.')
    parser.add_argument('-l','--layout', type=str, help='Layout for visualizing tree. Default: circular.')
    parser.add_argument('-c','--count', action="store_true",help='Count leaves')
    parser.add_argument('-n','--names', action="store_true",help='Option to replace names in nwk')
    parser.add_argument('-nf','--names_file', type = str ,help='File with names to replace. Default is second column in tab-separated format.')
    parser.add_argument('-pat','--pattern', type = str ,help='Pattern for selecting node names.')
    parser.add_argument('-al','--astral', action="store_true",help='Option to convert gene trees to ASTRAL input for species tree estimation.')
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()



def main():
    
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating multiple sequence alignments, newick files and trees.')
    args = parse_arguments()

    if args.tree and args.alignment:
        print ("Please select either --tree or --alignment input file.")
    
    elif args.tree:
        inp = re.sub (".nwk$","",args.tree)
        
        if args.midpoint:
            midpoint(args.tree)

        elif args.collapse:
            if args.cutoff:
                bootstrap_collapse(args.tree,args.cutoff)
            else:
                print ( "Please provide a bootstrap cutoff.")
            
        elif args.visualize:
            if args.layout:
                visualize_tree(args.tree, layout = args.layout)
            else:
                visualize_tree(args.tree)

        elif args.resolve:
            resolve_polytomies(args.tree)

        elif args.names:
            if args.names_file:
                if args.pattern:
                    with open ( f"{args.tree}.new_names.{args.pattern}.nwk", "w" ) as new_names:
                        print ( sub_names_nwk(tree = args.tree, file_with_names = args.names_file, pattern = args.pattern ), file = new_names )
                else:
                    with open ( f"{args.tree}.new_names.all_taxa.nwk", "w" ) as new_names:
                        print ( sub_names_nwk(tree = args.tree, file_with_names = args.names_file, pattern = False ), file = new_names )
            else:
                print ( "Please provide a list with gene names to replace.")

        elif args.count:
            print ( count_leaves(args.tree) )

        elif args.astral:
            prep_ASTRAL_input(args.tree)

    elif args.alignment:
        aln2phy(args.alignment, str(args.alignment) + ".phy")
            
    else:
        print ("Please provide a nwk or an alignment file as input.")



def identify_format(file_path):
    if re.search (".nwk$|.newick$", file_path):
        return 'newick'
    elif re.search (".fasta$|.fa$|.faa$|.fna$", file_path):
        return 'fasta'
    else:
        return 'list'


# Get support values
def extract_support_values(tree):

    """
    Extracts bootstrap support values from an ETE tree object.

    Parameters:
    - tree: Newick format string representing the ETE tree object.

    """
    
    support_values = []

    for node in tree.traverse():
        if hasattr(node, "support"):
            support_values.append(node.support)

    return support_values

    
# Main visualization function for leaf coloring and bootstrap support
def visualize_tree(tree, layout = "c", show = True):
    
    """
    Visualizes a phylogenetic tree with colored leaves representing species and bootstrap support indicated by circle size.

    Parameters:
    - tree: Newick format string representing the phylogenetic tree.
    - layout: Layout mode for tree visualization. Defaults to "c" (circular).
    - show: Boolean flag indicating whether to display the tree. Defaults to True.

    """
    
    t=ete3.Tree(tree)

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode = layout

    species_colors = {}

    for leaf in t.iter_leaves(): # Assign a unique color to each species 
        
        if re.search ( "_", leaf.name ): # If _ is present, get the prefix as species name
            species = re.sub("_.*","", leaf.name)

        else: # If no _ is found, get the species names from the first four letters of the gene ID
            species = leaf.name[0:3]
            
        if species not in species_colors:
            color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            if color != "#000000": # keep black color for special identifiers #
                species_colors[species] = color

    thresholds = {
        50: "grey",
        75: "darkgrey",
        100: "black"
    }

    for node in t.traverse():
        nstyle=NodeStyle()
        nstyle["size"] = 0
        node.set_style(nstyle)

        if node.is_leaf(): # If nodes are leaves, get color name based on species
            if re.search ( "_", node.name ): # If _ is present, get the prefix as species name
                species_n = re.sub("_.*","", node.name)
                
            else: # If no _ is found, get the species names from the first four letters of the gene ID
                species_n = node.name[0:3]
                
            color_n = species_colors[species_n]               
            species_face = TextFace(node.name,fgcolor=color_n, fsize=500,ftype="Arial")
            node.add_face(species_face, column=1, position='branch-right')

        for threshold, col in thresholds.items(): # Bootstrap values in internal nodes represented as circles
            if max(extract_support_values(t)) <= 10:
                print ( "Bootstrap support values are in the 0-1 scale and will be converted to the 1-100 scale.")
                node.support = node.support * 100
                
            if node.support <= 50:
                color_b = "lightgrey"
            elif node.support >= threshold:
                color_b = col
                
        if color_b:
            node_style=NodeStyle()
            node_style["fgcolor"] = color_b
            node_style["size"] = 500
            node.set_style(node_style)


    t.render(tree+".svg", w=500, units="mm", tree_style=ts)

    if show == True:
        t.show(tree_style=ts)

    return

        
## TREE MANIPULATION FUNCTIONS ##
def midpoint(input):
    
    """
    Sets the midpoint root as the outgroup for the phylogenetic tree.

    Parameters:
    - input: Newick format string representing the phylogenetic tree.

    """
    
    tree = ete3.Tree(input, format = 2)   
    midpoint = tree.get_midpoint_outgroup()

    # set midpoint root as outgroup #
    tree.set_outgroup(midpoint)
    tree.write(format=2, outfile=input+".tree")
    return


def bootstrap_collapse(tree, threshold=50):

    """
    Collapses tree branches with bootstrap support below a specified threshold.

    Parameters:
    - tree: Newick format string representing the phylogenetic tree.
    - threshold: Minimum bootstrap support value for branch collapse. Defaults to 50.

    """
    
    t=ete3.Tree(tree)
    for node in t.traverse():
        if node.support < threshold:
            return node.delete()
        else:
            return node

def resolve_polytomies(input):

    """
    Resolves polytomies in the phylogenetic tree.

    Parameters:
    - input: Newick format string representing the phylogenetic tree.

    """

    tree = ete3.Tree(input, format = 2)   
    tree.resolve_polytomy(recursive=True) # resolve polytomies in tree #
    tree.write(format = 2, outfile=input+".resolved_polytomies")
    return


# Leaf counting functions
def count_leaves ( tree ):

    """
    Counts the total number of leaves in the phylogenetic tree.

    Parameters:
    - tree: Newick format string representing the phylogenetic tree.

    Returns:
    - counts: Total number of leaves in the tree.

    """
    
    nleaves = []
    t = ete3.Tree(tree)
    
    for leaf in t.iter_leaves():
        nleaves.append(leaf)

    counts = len(nleaves)
    return counts


def count_descendant_leaves ( tree, node ):

    """
    Counts the number of descendant leaves for a given node in the phylogenetic tree.

    Parameters:
    - tree: Newick format string representing the phylogenetic tree.
    - node: Node name for which descendant leaves are counted.

    """
    
    t=ete3.Tree(tree)
    descendant_leaves = []

    for node in t.traverse ("preorder"):
        descendant_leaves.append ( node.get_leaf_names() )
        
    print ("You have", len(descendant_leaves), "leaves" )
    return


def count_leaves_by_taxon ( tree, taxon_ID ):

    """
    Counts the number of leaves associated with a specific taxon in the phylogenetic tree.

    Parameters:
    - tree: Newick format string representing the phylogenetic tree.
    - taxon_ID: Identifier for the taxon of interest.

    """
    
    t=ete3.Tree(tree)
    descendant_leaves = []

    for node in t.traverse ("preorder"):
        if node.is_leaf():
            if re.search ( taxon_ID, node.name ):
                descendant_leaves.append ( node.get_leaf_names() )
                
    print ("You have", len(descendant_leaves), "leaves for", taxon_ID )
    return


# Substitute taxon names in newick
def sub_names_nwk(tree, file_with_names, pattern=False):

    """
    Substitutes taxon names in a Newick format tree with names from a provided file.

    Parameters:
    - tree: Newick format string representing the phylogenetic tree.
    - file_with_names: Path to a file containing taxon names.
    - pattern: Regular expression pattern for matching taxon names in the tree. Defaults to False.

    Returns:
    - newick_string: Newick format tree with substituted taxon names.

    """
    
    n = ete3.Tree(tree)
    name = {}

    file_format = identify_format(file_with_names)

    if file_format == 'newick':
        if pattern == False:
            t = ete3.Tree(file_with_names)

            for node in t.traverse():
                if node.is_leaf():
                    if re.search(pattern, node.name):
                        name[node.name] = line.strip()[1:]

    elif file_format == 'fasta' or file_format == 'list':
        with open(file_with_names, "r") as file:
            lines = file.readlines()

        for line in lines:
            if file_format == 'fasta' and line.startswith('>'):
                name[line.strip()[1:]] = line.strip()[1:]
                
            elif file_format == 'list':
                name[line.strip()[1:]] = line.strip()[1:]
                
    for node in n.traverse():
        if node.is_leaf():
            if node.name in name.keys():
                print ( name[node.name] )
                node.name = str(name[node.name])

    return n.write(format = 2)

# Concatenate multiple gene trees together for ASTRAL input
def prep_ASTRAL_input (tree):

    """
    Prepares multiple gene trees for input into ASTRAL by concatenating them into a single file.

    Parameters:
    - tree: Newick format string representing the phylogenetic tree.

    """
    
    output_file = re.sub(".nwk|.tree|.tre",".astral.nwk",tree)
    
    with open(output_file, "a") as out:    
        with open ( tree, "r" ) as treefile:
            tl = treefile.readlines()

            for tre in tl:        
                t = ete3.Tree ( tre, quoted_node_names = True )

                for node in t.traverse():        
                    if node.is_leaf:
                        if re.search ( "_" , node.name ):
                            node.name = re.sub("_.*","",node.name)
                        elif re.search ( "." , node.name ):
                            node.name = re.sub("..*","",node.name)

                out.write (t.write(format=2) + "\n")


# Multiple Sequence Alignment (MSA) manipulation
def aln2phy ( input_file, output_file ):
    
    """
    Converts a multiple sequence alignment (MSA) file to PHYLIP format.

    Parameters:
    - input_file: Path to the input MSA file.
    - output_file: Path to the output PHYLIP format file.

    """
    
    sequences = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
        current_sequence = None
        for line in lines:
            if line.startswith('>'):
                if current_sequence is not None:
                    sequences.append(current_sequence)
                current_sequence = {'header': line.strip()[1:], 'sequence': ''}
            else:
                current_sequence['sequence'] += line.strip()
        if current_sequence is not None:
            sequences.append(current_sequence)

    num_sequences = len(sequences)
    seq_length = len(sequences[0]['sequence'])

    with open(output_file, 'w') as f:
        f.write(f"{num_sequences} {seq_length}\n")
        for sequence in sequences:
            f.write(f"{sequence['header']} {sequence['sequence']}\n")
            
               
if __name__ == "__main__":
    main()
