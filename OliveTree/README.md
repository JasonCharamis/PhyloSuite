# OliveTree
Library of functions for efficient tree manipulation and beautiful visualization in R, based on ggtree, treeio, ape and phytools.

Contains functions for:

1. Reading phylogenetic trees as treedata objects
2. Writing phylogenetic tree object as newick files while preserving bootstrap values
3. Collapsing branches based on bootstrap support
4. Producing the node labels of a tree with an option to print leaves matching a specific pattern
5. Flipping nodes based on user-provided descendant leaves
6. Grouping descendant leaves
7. Extracting subtrees based on the IDs of two user-provided anchor leaves, while preserving bootstrap values and branch lengths
8. Visualizing trees with highlighted descendant leaves of user-specified nodes
9. Visualizing trees with advanced features, including producing color and shape mappings for user-specified tip labels, isolating specific clades of tree and representing bootstrap support as colored circles at parent nodes of branches

To load it run:
source("OliveTree.R")
