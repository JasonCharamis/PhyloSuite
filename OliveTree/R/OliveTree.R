#' Library of functions for advanced tree manipulation and visualization using ggtree, ape, phytools, and other related tools.
#' @import treeio
#' @import ape
#' @import phytools
#' @import tidytree
#' @import ggstar
#' @import ggtree
#' @import dplyr
#' @import BiocManager
#' 
#' @title .load_packages
#' @description This function checks if a package is installed. If not, it installs the package using BiocManager if available, otherwise using install.packages.
#' @param tools A character vector of package names to be checked and installed.
#' @return NULL
#' @export

# Function to check if a package is installed, and if not, install it. 
# Then load it in memory.

.load_packages <- function(tools) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built, 1, 1)))
  tmp <- tmp[as.numeric(substr(tmp$Built, 1, 1)) == max_version,]
  
  for (pkg in tools) {
    if (pkg %in% tmp$Package) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    } else {
      print(sprintf("%s %s", pkg, "is not installed. Installing it!"))
      
      if (pkg %in% BiocManager::available(pkg)) {
        BiocManager::install(pkg, dependencies = TRUE, update = TRUE)
      } else {
        install.packages(pkg, dependencies = TRUE, ask = FALSE)
      }
    }
  }
}

# Load required packages or install them if necessary
dependencies <- c(
  "ape",
  "phytools",
  "treeio",
  "tidytree",
  "TreeTools",
  "ggstar",
  "ggtree",
  "dplyr",
  "ggplot2"
)

.load_packages(dependencies)

#==================================== TREE MANIPULATION FUNCTIONS ====================================#

#' read_tree
#' Read a phylogenetic tree from a file with auto-detection of bootstrap values.
#'
#' @param input_file Path to the file containing the phylogenetic tree. Input tree can be in newick or nexus format.
#' 
#' @return A tree data object representing the phylogenetic tree.
#' 
#' @export

# Function to read and preprocess a tree
read_tree <- function(input_file) {
  if (grepl("\\.(newick|nwk|tre|tree)$", tolower(input_file))) {
    t <- phytools::read.newick(input_file)
  } else if (grepl("\\.(nxs|nex)$", tolower(input_file))) {
      t <- ape::read.nexus(input_file)
  } else {
      stop("Unsupported file format.")
  }
  
  # Check if the tree contains bootstrap values, and if it contains load them
  if ( any(!is.null(t$node.label) && any(as.numeric(t$node.label) > 0)) == TRUE) {
    to <- treeio::as.treedata(t, node.label = "support")
  } else {
    to <- treeio::as.treedata(t)
  }
  return(to)
}


#' .load_tree_object
#' Load phylogenetic tree as treedata object
#' 
#' @param tree An object representing the phylogenetic tree, providing a single tree entrypoint used in the whole package. 
#'             Makes use of the read_tree function.
#'
#' If no file path is provided, the user will be prompted to provide an input phylogenetic tree through Rstudio's graphical user interface (GUI).
#' 
#' @return A tree data object representing the phylogenetic tree.
#' 
#' @export


.load_tree_object <- function(tree = NULL) {
  
  if ( is.null(tree) ) {
    tree <- read_tree(file.choose())
  } else if (typeof(tree) == "character") {
    if (file.exists(tree)) {
      if (any(grepl(".newick|.nwk|.tre|.support|.nxs|.nex", tree))) {
        tree_obj <- read_tree(tree)
      }
    } else {
      stop("Provided file is not a newick or nexus file.")
    }
  } else if (is(tree, "treedata")) {
      tree_obj <- tree
  } else if (is(tree, "phylo")) {
      tree_obj <- treeio::as.treedata(tree)
  } else {
      stop("Provided file is not a treedata, a phylo or a tibble_df object.")
  }
  return(tree_obj)
}


#' write_tree
#' This function writes a tree in newick format, while preserving the bootstrap values.
#' 
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo', or a file.
#' @param output Name of the output file. If no name, is provided it will use the tree object name with a '.nwk' suffix.
#' 
#' @return A ggplot object representing the phylogenetic tree with node IDs and optional highlighted tip labels.
#' 
#' @examples
#' \dontrun{
#' # Write a sample tree
#' write_tree( tree = tree_obj, output = 'tree_obj.nwk' )
#' }
#' @export

write_tree <- function(tree, output = "output_tree.nwk") {
  tree_obj <- .load_tree_object(tree)
  
  # Copy bootstrap values in the labels of the internal (non-Tip) nodes
  phylo_data <- treeio::as.treedata(
    left_join(
      as_tibble(tree_obj@phylo),
      as_tibble(tree_obj@data)
    ) %>%
      mutate(
        isTip = ifelse(!is.na(label), TRUE, FALSE)
      ) %>%
      mutate(
        label = ifelse(
          isTip == FALSE,
          as.character(support),
          label
        )
      ) %>%
      select(-support)
  )
  
  ape::write.tree(phy = phylo_data@phylo, output)
}

#' export_plot: Function to export plot with provided options
#' 
#' 
#' @param plot A ggtree plot object.
#' @param save Option to save plot to output svg file. Default: True
#' @param output Name of the output file. Default: 'tree_plot_visualized.svg'
#' 
#' @return A ggplot object representing the phylogenetic tree and optionally an output svg file with this phylogenetic tree.
#'         Main plot exporting function of the package. 
#' 
#' @examples
#' dontrun{
#' # Export plot with option to save
#' export_plot (plot, save = TRUE, output = "tree_plot_visualized.svg")
#' }
#' @export

export_plot <- function(
    plot,
    save = TRUE,
    output = "tree_plot_visualized.svg"
) {
  
  if (!(is.null(output))) {
    save == TRUE
  }
  
  # Export plot with provided options
  if (save == TRUE) {
    if (is.null(output)) {
      if (typeof(tree) == "character") {
        if (file.exists(tree)) {
          if (any(grepl(".newick|.nwk|.tre|.support|.nxs|.nex", tree))) {
              
            ggsave(
                plot = plot, 
                sprintf( "%s_visualized.svg", sub(".newick|.nwk|.tre|.support|.nxs|.nex", "", tree)), 
                dpi = 600
              )
              
              message <- sprintf("Tree plotted and saved as %s_visualized.svg", sub(".newick|.nwk|.tre|.support|.nxs|.nex", "", tree))
              print(message)
            }
          }
      } else {
            ggsave(plot = plot, "tree_plot_visualized.svg", dpi = 600)
            print("Tree plotted and saved as tree_plot_visualized.svg!")
      }
    } else {
          ggsave(plot = plot, output, dpi = 600)
          print(paste("Tree plotted and saved as", output))
    }
  } else {
        print("Plot will not be saved! Use the options save = TRUE and output = <OUTPUT_NAME> for saving the output plot.")
  }
    
  return(plot)
}


#' print_internal_nodes: This function generates a plot of the provided phylogenetic tree with node IDs displayed. 
#'                       Additionally, it offers the option to highlight specific tip labels that match a user-provided pattern.
#'
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo'.
#' @param form Layout of the tree based on ggtree options. Layout can be rectangular, circular, roundrect, slanted, ellipse, fan, equal_angle, daylight (Default: "rectangular").
#' @param node_id_color Color of the printed node ids. Default is "darkred".
#' @param ... Pattern(s) for printing tip labels that match them.
#' 
#' @return A ggplot object representing the phylogenetic tree with node IDs and optional highlighted tip labels.
#'         Used for identifying nodes for interest for extracting subtrees, flipping or labeling them in the visualize_tree function.
#' 
#' @examples
#' \dontrun{
#' # Load a sample tree
#' tree <- read_tree("path/to/tree.nwk")
#'
#' # Print the tree with node IDs and print tip labels matching "taxon_A"
#' print_internal_nodes(tree, color = "darkred", "taxon_A")
#' }
#' 
#' @export

print_internal_nodes <- function(
    tree, 
    form = "circular", 
    node_id_color = "red", 
    tip_label_size = 2, 
    references = NULL,
    save = TRUE,
    output = "highlighted_tree.svg"
) {
  
  tree_obj <- .load_tree_object(tree)
  
  # Create the base tree plot with node labels
  plot <- ggtree(tree_obj, layout = form) +
          geom_nodelab(aes(label = node), hjust = -0.1, color = node_id_color, size = 3)
  
  if ( is.null(references) ) {
    print("Only node IDs will be printed.")
    return(plot)
  } else {
    matching_labels <- unlist(
      sapply(references, function(reference) {
        grep(reference, tree_obj@phylo$tip.label, value = TRUE, fixed = TRUE)
      }))
    
    matching_dict <- data.frame(
      label = matching_labels,
      matching_flag = TRUE
    )
    
    plot <- plot %<+% matching_dict
    
    plot <- plot + 
      geom_tiplab(aes(label = ifelse(matching_flag == TRUE, label, NA), 
                  fill =  sample(colors(), 1)), 
                  size = tip_label_size, align.tip.label = TRUE)
  }
  
  export_plot(
    plot, 
    save = save,
    output = output
  )
  
}

#' bootstrap_collapse
#' This function collapses nodes in a phylogenetic tree where the bootstrap support is below a specified cutoff.
#'
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo'.
#' @param cutoff The threshold value for collapsing nodes based on bootstrap support. Defaults to 0.5.
#'
#' @return A 'phylo' object with collapsed branches with bootstrap less than the provided cutoff.
#'
#' @examples
#' \dontrun{
#' # Load a sample tree
#' tree <- read_tree("path/to/tree.nwk")
#'
#' # Collapse nodes with bootstrap support less than 0.7
#' collapsed_tree <- bootstrap_collapse(tree, 0.7)
#' }
#'
#' @export

bootstrap_collapse <- function(tree, cutoff = 0.5) {
  tree_obj <- .load_tree_object(tree)
  return(as.polytomy(tree_obj, feature = 'support', fun = function(x) as.numeric(x) < cutoff))
}

#' group_descendants
#'
#' Group all descendant branches of specified node(s) in a phylogenetic tree.
#'
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo'.
#' @param node1 The first node for grouping descendants.
#' @param node2 The second node for grouping descendants.
#' @param node3 The third node for grouping descendants.
#' @param node4 The fourth node for grouping descendants.
#'
#' @return A 'phylo' object with grouped descendants.
#'
#' @examples
#' \dontrun{
#' # Load a sample tree
#' tree <- read_tree("path/to/tree.nwk")
#'
#' # Group all descendants of nodes 5, 8, and 10 in the tree
#' grouped_tree <- group_descendants(tree, 5, 8, 10)
#' }
#'
#' @export

# Function to group all descendant branches of node(s)
group_descendants <- function(tree, ...) {
  tree_obj <- .load_tree_object(tree)
  nodes <- list(...)
  return(tidytree::groupClade(tree_obj, .node = nodes))
}

#' extract_subtree
#' Function to extract a subtree by finding the MRCA (Most Recent Common Ancestor) of two anchor leaves or 
#' the subtree that contains a user-provided list of tip labels,  while preserving branch lengths and bootstrap values.
#'
#' @param tree A phylogenetic tree object of class 'phylo' or 'treedata'.
#' @param tip1,tip2 The labels or pattern of the two anchor tips.
#' @param taxon_list User-provided list of tip labels, the subtree which contains them will be extracted.
#'
#' @return A treedata object representing the extracted subtree with preserved branch lengths and bootstrap values.
#'
#' @examples
#' # Extract a subtree preserving branch lengths and bootstrap values
#' tree <- read.tree("path/to/your/treefile.newick")
#' subtree <- extract_subtree(tree, "species_A", "species_B", branch_length = TRUE)
#' visualize_tree(subtree)
#'
#' # Extract a cladogram with equal branch lengths
#' tree <- read.tree("path/to/your/treefile.newick")
#' subtree <- extract_subtree(tree, "tip1", "tip")
#' visualize_tree(subtree)
#'
#' @seealso
#' \code{\link{visualize_tree}} for visualizing phylogenetic trees.
#' \code{\link{findMRCA}} for finding the Most Recent Common Ancestor (MRCA) of multiple nodes.
#'
#' @importFrom dplyr left_join select mutate as_tibble
#' @importFrom ape extract.clade
#' @importFrom phytools findMRCA 
#' @importFrom treeio as.treedata read.newick
#'
#' @export

# Function to extract a subtree by finding the MRCA of two anchor nodes while preserving branch lengths and bootstrap values 

extract_subtree <- function(
    tree, 
    tip1 = NULL, 
    tip2 = NULL, 
    taxon_list = NULL
) {
  
  tree_obj <- .load_tree_object(tree)
  
  if ("support" %in% colnames(tree_obj@data)) {
    
    # Append bootstrap values in the label section, in the nodes which are NOT tips
    # The bootstrap values will be kept there and transferred along with the labels in the new node numbering of the subtree
    phylo_data <- treeio::as.treedata(
      left_join(
        as_tibble(tree_obj@phylo),
        as_tibble(tree_obj@data)
      ) %>%
        mutate(
          isTip = ifelse(!is.na(label), TRUE, FALSE)
        ) %>%
        mutate(
          label = ifelse(isTip == FALSE, as.character(support), label)
        ) %>%
        select(-support)
    )
  } else {
      phylo_data <- tree_obj
  }
  
  t <- TreeTools::Preorder(phylo_data@phylo)
  
  if (is.null(taxon_list) && !is.null(tip1) && !is.null(tip2)) {
    # Extract subtree based on anchor tip labels and/or patterns
    tip1m <- t$tip.label[grepl(tip1, t$tip.label)]
    tip2m <- t$tip.label[grepl(tip2, t$tip.label)]
    
    if (any(sapply(list(tip1m, tip2m), function(x) is.null(x)))) {
      stop("Provided tip labels are NULL or were not found among tree tip labels.")
    } else if (any(sapply(list(tip1m, tip2m), function(x) length(x) > 1))) {
        stop("Provided tip names are not unique.")
    }
    
    subtree <- ape::extract.clade(t, node = phytools::findMRCA(t, c(tip1m, tip2m)))
    
  } else if ( is.null(tip1) && is.null(tip2) && !is.null(taxon_list)) {
    taxon_tips <- t$tip.label %in% taxon_list
    
    if (any(taxon_tips)) {
      mrca <- getMRCA(t, which(taxon_tips))
      subtree <- ape::extract.clade(t, mrca)
    } else {
        stop("No matching species found in the tree.")
        return(NULL)
    }
  } else {
      stop("Please choose either tip1,tip2 OR taxon_list as identifiers for extracting subtree.")
  }
  
  if (!("support" %in% colnames(tree_obj@data))) {
    return(subtree)
  } else {
      # If node is NOT tip, add the labels in the bs_support column and make the labels NA - mapping according to new node numbering
      subtree_f <- treeio::as.treedata(subtree) %>%
        mutate(
          bs_support = ifelse(isTip == FALSE, label, NA)
        ) %>%
        mutate(
          label = ifelse(isTip == FALSE, NA, label)
        )
      
      # To create an object compatible with visualize_tree, initialize the 'node' and 'support' columns in the data slice of the treedata class
      subtree_f@data <- tibble(
        node = rep(NA, length = length(subtree_f@extraInfo$node)),
        support = rep(NA, length = length(subtree_f@extraInfo$bs_support))
      )
      
      # ... and then add the node number and bs_support values
      subtree_f@data <- tibble(
        node = subtree_f@extraInfo$node,
        support = as.numeric(subtree_f@extraInfo$bs_support)
      )
      return(subtree_f)
  }
}

# ==================================== TREE VISUALIZATION FUNCTIONS ====================================#

#' highlight_tree: Highlight Nodes and Descendant Leaves on a Phylogenetic Tree
#'
#' This function highlights specified nodes and their descendant leaves on a phylogenetic tree,
#' allowing for visualization emphasis on selected parts of the tree. Custom colors can be assigned
#' to different groups of nodes, or random colors can be automatically generated.
#'
#' @param tree A phylogenetic tree object, typically of class 'phylo' from the ape package,
#'        or a file path to a tree file in Newick format that can be read into such an object.
#' @param highlight_nodes An associative vector where names are group labels and values are
#'        node IDs indicating the nodes and their descendants to be highlighted on the tree.
#' @param colors An associative vector where names are group labels and values are colors
#'        for the corresponding groups defined in `highlight_nodes`. If not provided,
#'        random colors will be assigned to each group.
#' @param form A character string specifying the layout of the tree. Common values include
#'        "rectangular", "circular", and others supported by ggtree.
#' @param tip_fill_color Fill color for the tip label symbols. Default: lightgrey
#' @param tip_border_color Border color for the tip label symbols. Default: black
#' @param legend_name A character string specifying the legend name for the highlighted plot. This parameter
#'        is optional and primarily used for labeling within complex workflows.
#' @param legend_position Position of the legend in the plot ("right", "left", "top", "bottom").
#' @param legend_text_position Horizontal alignment of the text within the legend ("left" or "right").
#' @param legend_orientation Orientation of the legend ("horizontal" or "vertical").
#' @param legend_key_size Size of the legend keys, specified in cm.
#' @param legend_font Font family for the text in the legend.
#' @param legend_fontsize Size of the text in the legend, specified in points.
#' @param legend_font_face Boolean indicating whether the text in the legend should be bold.
#' @param legend_spacing_x Horizontal spacing between legend items, specified in cm.
#' @param legend_spacing_y Vertical spacing between legend items, specified in cm.
#' @param legend_key_width Width of the keys in the legend, specified in cm.
#' @param legend_title_hjust Horizontal justification of the legend title (0 to 1).
#' @param save Logical; if TRUE, the plot will be saved to the file specified in `output`.
#' @param output The filename where the plot will be saved if `save` is TRUE.
#'
#' @return A ggplot object representing the phylogenetic tree with highlighted nodes.
#'
#' @examples
#' \dontrun{
#' library(ggtree)
#' library(treeio)
#' # Load a tree from a Newick file
#' tree <- read.tree("path/to/your/treefile.newick")
#' # Highlight specific nodes with specified colors
#' highlight_tree(tree,
#'                highlight_nodes = c("Node1" = 5, "Node2" = 10),
#'                colors = c("Node1" = "red", "Node2" = "blue"),
#'                form = "circular")
#' }
#'
#' @seealso \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{ggsave}},
#'          \code{\link[ggtree]{ggtree}}, \code{\link[ggtree]{geom_hilight}}
#' @importFrom ggplot2 aes geom_point scale_fill_manual theme ggsave
#' @importFrom ggtree ggtree geom_hilight
#' @importFrom treeio read.newick
#' @export

# Function to highlight nodes on a tree
highlight_tree <- function(
    tree,
    highlight_nodes,
    legend_position = "right",
    legend_text_position = "right",
    legend_orientation = "horizontal",
    legend_key_size = 0.1,
    legend_font = "Arial",
    legend_fontsize = 11,
    legend_font_face = "bold",
    legend_spacing_x = 0.5,
    legend_spacing_y = 0.5,
    legend_key_width = 1,
    legend_title_hjust = 0.5,
    tip_fill_color = "lightgrey",
    tip_border_color = "white",
    colors = NULL,
    form = "circular",
    legend_name = "",
    name = NULL,
    save = TRUE,
    output = "highlighted_tree.svg"
) {
  
  tree_obj <- .load_tree_object(tree)
  
  plot <- ggtree(tree_obj, layout = form) +
    geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) +
    scale_starshape_identity() 
  
  # Create a data frame for highlighting
  if (is.vector(highlight_nodes)) {
    if (is.null(colors)) {
      colors <- setNames(sample(rainbow(length(highlight_nodes)), length(highlight_nodes)), names(highlight_nodes))
    }
    
    highlight <- data.frame(
      Group = names(highlight_nodes),
      Label = highlight_nodes,
      Color = colors
    )
  } else {
      stop("highlight_nodes should be a vector.")
  }
  
  plot <- plot +
    geom_hilight(
      data = highlight,
      aes(node = Label, fill = Color),
      alpha = 0.3,
      extend = 0.10,
      linetype = 1,
      linewidth = 0.9
    ) +
      scale_fill_manual(
        values = highlight$Color,
        labels = gsub("_", " ", highlight$Group),
        name = legend_name
    ) + 
    theme(
      legend.direction = legend_orientation,
      legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
      legend.key.size = unit(legend_key_size, "cm"),
      legend.spacing.x = unit(legend_spacing_x, "cm"),
      legend.spacing.y = unit(legend_spacing_y, "cm"),
      legend.key.width = unit(legend_key_width, "cm"), 
      legend.title = element_text(hjust = legend_title_hjust)
    )
  
  export_plot(
    plot, 
    save = save,
    output = output
  )
  
}


#' visualize_tree: Visualize Phylogenetic Tree with Various Customizable Features
#'
#' Main tree visualization function of the package, which provides comprehensive options 
#' for visualizing a phylogenetic tree, allowing for customizations such as 
#' displaying bootstrap values, tip labels as text and shape with different colors based on user-provided mappings, 
#' flipping nodes and label user-specified clades. It supports various tree layouts.
#'
#' @param tree A phylogenetic tree in Newick format, treedata, or phylo class object.
#' @param form The layout of the tree, with options including "rectangular", "circular",
#'        "roundrect", "slanted", "ellipse", "fan", "equal_angle", and "daylight". Default is "rectangular".
#' @param flip_nodes Logical; whether to flip the tree at specified nodes.
#' @param node1,node2 The ID of the two nodes to flip based on their most common recent ancestor (MRCA)
#' @param tiplabels Logical; whether to print tip labels on the tree. Default: FALSE.
#' @param pattern_id A pattern used to match and selectively print tip labels as text.
#' @param tip_label_size Numeric specifying the size of the tip labels.
#' @param tip_label_color Color for the tip labels, default is black.
#' @param taxon_group_separator A separator string for splitting taxon IDs to generate color
#'        and shape mappings. If NULL, mappings are generated based on the whole taxon ID string. 
#' @param color A vector specifying color mappings for tip label shapes based on the taxon ID. 
#'              Is different to tip_label_color, to allow for a user to provide colors for both text and shape mappings.
#' @param shape A vector specifying ggstar starshapes for tip labels based on the taxon ID. 
#' @param tip_shape_size Numeric specifying the size of the ggstar starshapes representing tip labels.
#' @param legend_position The position of the legend on the plot ("right" by default).
#' @param legend_text_position The position of the text in the legend ("right" by default).
#' @param legend_orientation The orientation of the legend, either "horizontal" or "vertical" (default).
#' @param legend_key_size Size of the keys in the legend, in cm.
#' @param legend_font The font family for the legend text.
#' @param legend_fontsize The size of the font in the legend, in points.
#' @param legend_font_face Whether the font in the legend should be bold ("bold" or normal).
#' @param legend_spacing_x Horizontal spacing in the legend, in cm.
#' @param legend_spacing_y Vertical spacing in the legend, in cm.
#' @param legend_key_width The width of the legend keys, in cm.
#' @param legend_title_hjust Horizontal justification of the legend title (0 to 1).
#' @param legend_name The name to be displayed as the legend title.
#' @param bootstrap_numbers Logical; whether to display bootstrap values on branches.
#' @param bootstrap_number_nudge_y Numeric value controlling the vertical offset of bootstrap numbers.
#' @param bootstrap_circles Logical; whether to display bootstrap values as circles at nodes.
#' @param bootstrap_legend Logical; whether to display a legend for bootstrap values if displayed as circles.
#' @param bootstrap_circle_size Numeric specifying the size of the bootstrap circles.
#' @param branch_length Logical indicating whether branch lengths should be shown or omitted. Default: TRUE
#' @param clades A vector or list specifying clades and their labels for highlighting.
#' @param labeldist The distance for clade labels from their respective nodes.
#' @param bardist The distance for clade label bars from their respective nodes.
#' @param clade_label_size The size of the clade labels.
#' @param save Logical; if TRUE, the plot will be saved using the specified output file name.
#' @param output The file name or path where the plot will be saved if `save` is TRUE.
#'
#' @return A ggplot object representing the visualized tree.
#'
#' @examples
#' \dontrun{
#' library(ggtree)
#' library(ggplot2)
#' tree <- rtree(50)  # Generate a random tree
#' visualize_tree(tree, form = "circular", bootstrap_circles = TRUE,
#'                color = c("Species1" = "red", "Species2" = "blue", "Species3" = "orange"),
#'                shape = c("Species1" = "circle", "Species2" = "star", "Species3" = "triangle"),
#'                tiplabels = TRUE, pattern_id = "Species")
#' }
#'
#' @seealso \link[ggtree]{ggtree}, \link[ggplot2]{geom_point}
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual theme
#' @importFrom ggtree ggtree geom_tiplab geom_nodelab
#' @export


visualize_tree <- function(
    tree = NULL,
    form = "rectangular",
    flip_nodes = FALSE,
    node1, node2 = NULL,
    tiplabels = FALSE,
    pattern_id = NULL,
    color = NULL,
    shape = NULL,
    taxon_group_separator = "NULL",
    legend_position = "right",
    legend_text_position = "right",
    legend_orientation = "horizontal",
    legend_key_size = 0.1,
    legend_font = "Arial",
    legend_fontsize = 11,
    legend_font_face = "bold",
    legend_spacing_x = 0.5,
    legend_spacing_y = 0.5,
    legend_key_width = 1,
    legend_title_hjust = 0.5,
    legend_name = "",
    tip_shape_size = 3,
    tip_label_size = 2,
    tip_label_color = "black",
    bootstrap_numbers = FALSE,
    bootstrap_number_nudge_x = 0,
    bootstrap_number_nudge_y = 0.2,
    node_label_size = 3, 
    bootstrap_circles = FALSE,
    bootstrap_legend = FALSE,
    bootstrap_circle_size = 1.7,
    branch_length = TRUE,
    clades = NULL,
    labeldist = 1,
    bardist = 1,
    clade_label_size = 3,
    save = TRUE,
    output = NULL
) {
  
  # Load tree with option to print branch length and/or flip nodes
  tree_obj <- .load_tree_object(tree)
  
  # If option for branch length is not TRUE, make branch lengths NULL
  if (!(branch_length == TRUE | branch_length == T)) {
    print("Option branch length is deactivated. As a result, the subtree will be extracted as a cladogram with equal branch lengths.")
    tree_obj@phylo$edge.length <- NULL
  } 
  
  if ( flip_nodes == TRUE ) {
    if ( !is.null(node1) || !is.null(node2) ) {
      plot <- ggtree::flip( ggtree(tree_obj, layout = form), node1, node2 )
    } else {
        stop ("Please provide node IDs to flip. Find node IDs of interest using the print_internal_nodes( tree, references = c(taxon_pattern1, taxon_pattern1) ) function.")
    }
  } else {
      plot <- ggtree(tree_obj, layout = form)
  }

  # Manipulate bootstrap values and control if and how to print them
  if (any(!is.null(tree_obj@data)) && any(!is.na(tree_obj@data))) {
    if ( all(tree_obj@data$support <= 10, na.rm = TRUE) == TRUE) { # Bootstrap values are in scale of 1. 
      warning("Provided bootstrap values were in scale of 1 and were multiplied by 100 to change scales.")
      
      tree_obj@data$support <- sapply(
        tree_obj@data$support, function(x) {
          x <- x * 100
        }
      )
    } 
    
    tree_obj <- treeio::as.treedata (
      as_tibble(tree_obj) %>%
        mutate (
          bs_color = case_when(
            support < 50 ~ "snow2",
            support >=   50 & support < 75 ~ "grey",
            support >=   75 ~ "black",
            TRUE ~ NA 
          )
        )
    )
    
    if (bootstrap_numbers == TRUE) {
      plot <- plot + geom_nodelab(
        aes(label = round(support, 2)), 
        nudge_y = bootstrap_number_nudge_y,
        size = node_label_size
      )
      print("Bootstrap values displayed as labels.")
      
     } else if (bootstrap_circles == TRUE) {
        plot <- plot + geom_point(
          aes(size = support),
          shape = 21,
          fill = 'gray',
          alpha = 0.6
        ) + scale_size_continuous(range = c(1, bootstrap_circle_size))
        print("Bootstrap values displayed as circles.")
        
     } else if ( bootstrap_numbers == TRUE && bootstrap_circles == TRUE ) {
        stop("Please select either bootstrap_numbers or bootstrap_circles.")
       
     } else {
        message("Bootstrap values will not be printed.")
       
     } 
    
  } else {
      message ("Bootstrap support was not found in tree.")
  }
  
  # Options for printing tip labels as text
  if (tiplabels == TRUE) {
    if (is.null(pattern_id)) {
      plot <- plot + 
        geom_tiplab(
          size = tip_label_size, 
          color = ifelse(is.null(color), "black", color), 
          align.tip.label = TRUE, 
        )
    } else {
        matching_labels <- unlist(
          sapply(pattern_id, function(pattern) {
            grep(pattern, tree_obj@phylo$tip.label, fixed = TRUE, value = TRUE)
          }), use.names = FALSE
        )
        
        if (length(matching_labels) > 0) {
          matching_dict <- data.frame(
            label = matching_labels,
            color = tip_label_color[matching_labels], 
            matching_flag = TRUE
          )
          
          plot <- plot %<+% matching_dict 
          
          plot <- plot + 
            geom_tiplab(
              aes(
                subset = isTip,
                label = label,
                color = ifelse(matching_flag == TRUE, color, NA)
              ),
              size = tip_label_size, 
              align.tip.label = TRUE
            ) + scale_color_identity()
          
        if (length(pattern_id) < 10 ) {
          cat("Printing phylogenetic tree plot with tip labels matching", pattern_id, "!")
        } else {
            print("Printing phylogenetic tree plot with >10 tip label patterns!")
        }
        
      } else {
          message( cat("Provided", pattern_id, "was not found among tip labels! Printing phylogenetic tree plot without tip labels!") )
      }
    }
  } else { 
      if (!is.null(pattern_id) || !is.null(tip_label_color)) {
        stop ("Please enable tip label printing with tiplabels = TRUE.")
      } else {
          print ("Printing phylogenetic tree plot without tip labels!")
    }
  }
  
  # Visualize phylogenetic tree with color and shape mappings, and an option to print selected taxa as tip labels
  
  # If color and/or shape mappings are present for tip labels
  if (!is.null(color) || !is.null(shape)) {
    
    # Use the taxon group separator to generate mappings
    if ( !is.null(taxon_group_separator) ) {
      taxa_dict <- lapply(tree_obj@phylo$tip.label, function(label) {
        list( tip_label = label, group = sub(paste0(taxon_group_separator,"_.*"), "", label)) 
      })
    } else {
      taxa_dict <- lapply(tree_obj@phylo$tip.label, function(label) {
        list(tip_label = label, group = label)
        }
      )
    }
  
    # If tip is also defined in the color and/or shape mappings, the mappings will override the text tip label
    if (length(taxa_dict) > 0) {
      tip_labels <- sapply(taxa_dict, function(entry) entry$tip_label)
      taxa_names <- sapply(taxa_dict, function(entry) entry$group)
      
      if (!is.null(color)) {
        missing_species_color <- setdiff(unique(sapply(taxa_dict, function(entry) entry$group)), names(color))
        
        if (length(missing_species_color) > 0) {
          warning(paste("Color mapping not found for the following species: ", paste(missing_species_color, collapse = ", ")))
        }
        
        tip_colors_df <- data.frame(
          label = tip_labels,
          taxa_colors = taxa_names,
          s_color = color[match(sapply(taxa_dict, function(entry) entry$group), names(color))]
        )
        
        if (!is.null(shape)) {
          if (tiplabels == TRUE) {
            warning("If both tip labels and shapes are enabled, shape mapping will override the text tip labels provided." )
          }
          
          missing_taxa_shape <- setdiff(unique(sapply(taxa_dict, function(entry) entry$group)), names(shape))
          
          if (length(missing_taxa_shape) > 0) {
            warning(paste("Shape mapping not found for the following species: ", paste(missing_taxa_shape, collapse = ", ")))
          }
          
          tip_shapes_df <- data.frame(
            label = tip_labels,
            taxa_shapes = taxa_names,
            s_shape = shape[match(sapply(taxa_dict, function(entry) entry$group), names(shape))]
          )
        }
        
        # Add ggtree aesthetics to plots based on the user-provided arguments
        if ( !is.null(color) && !is.null(shape)) {
          
          # Merge color and shape mappings into a single dataframe to allow for mappings legend
          tip_mapping_df <- inner_join(
            tip_colors_df, 
            tip_shapes_df, 
            by = "label"
          ) %>%
          distinct()
          
          # Join data frames to create a unique mapping tibble
          plot <- plot %<+% tip_mapping_df 
  
          plot <- plot + geom_star(
            mapping = aes( 
              subset = isTip,
              fill = ifelse(!is.na(taxa_colors), s_color, NA),
              starshape = ifelse(!is.na(taxa_shapes), s_shape, NA),
            ), size = tip_shape_size,
            show.legend = mappings_legend 
          ) + 
            scale_fill_identity( 
              guide = guide_legend(theme = theme(
                legend.direction = legend_orientation,
                legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
                legend.key.size = unit(legend_key_size, "cm"),
                legend.spacing.x = unit(legend_spacing_x, "cm"),
                legend.spacing.y = unit(legend_spacing_y, "cm"),
                legend.key.width = unit(legend_key_width, "cm"), 
                legend.title = element_text(hjust = legend_title_hjust)
              )),
              label = gsub("_", " ", tip_colors_df$taxa_colors),
              name = legend_name
            ) + 
            scale_starshape_identity( 
              guide = guide_legend(theme = theme(
                legend.direction = legend_orientation,
                legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
                legend.key.size = unit(legend_key_size, "cm"),
                legend.spacing.x = unit(legend_spacing_x, "cm"),
                legend.spacing.y = unit(legend_spacing_y, "cm"),
                legend.key.width = unit(legend_key_width, "cm"), 
                legend.title = element_text(hjust = legend_title_hjust)
              )),
              label = gsub("_", " ", tip_colors_df$taxa_colors),
              name = legend_name
            )
          
        } else if ( !is.null(color) && is.null(shape)) {
          
            plot <- plot %<+% tip_colors_df
            
            plot <- plot + geom_star(
              mapping = aes(
                subset = isTip,
                fill = ifelse(!is.na(taxa_colors), s_color, NA)
              ), 
              starshape = "circle",
              size = tip_shape_size
            ) + 
              scale_fill_identity( 
                guide = guide_legend(theme = theme(
                  legend.direction = legend_orientation,
                  legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
                  legend.key.size = unit(legend_key_size, "cm"),
                  legend.spacing.x = unit(legend_spacing_x, "cm"),
                  legend.spacing.y = unit(legend_spacing_y, "cm"),
                  legend.key.width = unit(legend_key_width, "cm"), 
                  legend.title = element_text(hjust = legend_title_hjust)
                )),
                label = gsub("_", " ", tip_colors_df$taxa_colors),
                name = legend_name
              )
            
          } else if ( is.null(color) && !is.null(shape)) {
          
              plot <- plot %<+% tip_shapes_df
            
              plot <- plot + geom_star(
                mapping = aes( 
                  subset = isTip & !is.na(s_shape),
                  starshape = ifelse(!is.na(taxa_shapes), s_shape, NA)
                ), 
                fill = "black", 
                size = tip_shape_size
              ) + 
                scale_starshape_identity(
                  guide = guide_legend(theme = theme(
                    legend.direction = legend_orientation,
                    legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
                    legend.key.size = unit(legend_key_size, "cm"),
                    legend.spacing.x = unit(legend_spacing_x, "cm"),
                    legend.spacing.y = unit(legend_spacing_y, "cm"),
                    legend.key.width = unit(legend_key_width, "cm"), 
                    legend.title = element_text(hjust = legend_title_hjust)
                  )),
                  label = gsub("_", " ", tip_colors_df$taxa_colors),
                  name = legend_name
              )
        }
      }
    } else {
        stop("None of the defined mapping tips was found in the tree.")
    }
  }
  
  # Option to label specific clades of interest 
  if (!is.null(clades)) {
    for (i in names(clades)) {
      plot <- plot + geom_cladelab(
        node = clades[i], 
        label = i, 
        align = TRUE, 
        fill = 'black',
        offset.text = labeldist, 
        barsize = 0.9, 
        offset.bar = bardist, 
        clade_label_size = clade_label_size
      )
    }
  }
  
  export_plot(
    plot = plot, 
    save = save,
    output = output
  )
}