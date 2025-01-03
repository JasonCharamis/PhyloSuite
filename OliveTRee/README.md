# OliveTree: R package for efficient tree manipulation and beautiful visualization based on ggtree, treeio and phytools.

Advanced phylogenetic tree manipulation and visualization in R

## Overview

OliveTree extends and unifies functionality from `ggtree`, `treeio`, `phytools`, and other phylogenetic tools to provide:

- Streamlined tree manipulation
- Enhanced visualization capabilities
- Robust bootstrap support handling
- Flexible clade management

## Installation

```R
# Install dependencies first
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ggtree", "treeio"))
install.packages(c("ape", "phytools"))

# Then source OliveTree
source("OliveTree.R")
```

## Key Features

### Tree Input/Output

- Read trees from Newick, Nexus and IQ-TREE formats with automatic bootstrap detection
- Export trees with preserved bootstrap values
- Unified tree object handling via `treedata` format

### Tree Manipulation

- Extract subtrees while preserving branch lengths and support values
- Collapse branches below specified bootstrap thresholds
- Identify and flip nodes based on descendant tips
- Group and manipulate descendant tips
- Print node labels with pattern-based tip highlighting

### Visualization

- Advanced tree layouts (rectangular, circular, fan etc.)
- Bootstrap value display options:
  - Numeric labels
  - Color-coded circles
  - Custom thresholds
- Tip label customization:
  - Color mapping
  - Shape mapping (using `ggstar` shapes)
  - Pattern-based highlighting
- Clade manipulation:
  - Highlight specific clades
  - Add clade labels and bars
  - Group descendants
- Legend customization options

## Usage Examples

### Basic Tree Reading and Visualization

```R
# Read tree
tree <- read_tree("path/to/tree.nwk")

# Basic visualization
visualize_tree(tree, form = "circular")

# With bootstrap circles
visualize_tree(tree, 
              bootstrap_circles = TRUE,
              bootstrap_circle_size = 2)
```

### Advanced Visualization

```R
# Color and shape mapping
visualize_tree(tree,
    color = c("SpeciesA" = "red", "SpeciesB" = "blue"),
    shape = c("SpeciesA" = "star", "SpeciesB" = "circle"),
    legend_name = "Species")

# Highlight clades
visualize_tree(tree,
    highlight_clades = c("CladeA" = 15, "CladeB" = 20),
    highlight_colors = c("CladeA" = "pink", "CladeB" = "lightblue"))
```

### Tree Manipulation

```R
# Extract subtree
subtree <- extract_subtree(tree, tip1 = "taxonA", tip2 = "taxonB")

# Collapse low support branches
collapsed_tree <- bootstrap_collapse(tree, cutoff = 70)

# Print node IDs
print_internal_nodes(tree, pattern = "taxonA")
```

## Dependencies

- ape
- phytools
- treeio
- tidytree
- TreeTools
- ggstar
- ggtree
- dplyr
- ggplot2

## Citation

If you use OliveTree in your research, please cite:

```
@software{olivetree2024,
  author = {Original Author},
  title = {OliveTree: Advanced Phylogenetic Tree Manipulation and Visualization},
  year = {2024},
  url = {https://github.com/username/olivetree}
}
```

## Contributing

Contributions are welcome! Please submit issues and pull requests.
