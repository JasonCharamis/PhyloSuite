#' Library of functions for advanced tree manipulation and visualization using ggtree, ape, phytools, and other related tools.
#' @import treeio
#' @import ape
#' @import phytools
#' @import tidytree
#' @import TreeTools
#' @import ggstar
#' @import ggtree
#' @import dplyr
#' @import plotly
#' @import BiocManager
#' 
#' @title .load_packages
#' @description This function checks if a package is installed. If not, it installs the package using BiocManager if available, otherwise using install.packages.
#' @param tools A character vector of package names to be checked and installed.
#' @return NULL
#' @export

# Function to check if a package is installed, and if not, install it.
# Function to install or load a package

.load_packages <- function(tools) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built,
                                      1,
                                      1)))
  tmp <- tmp[as.numeric(substr(tmp$Built,
                                1,
                                1)) == max_version,
             ]

  for (pkg in tools) {
    if (pkg %in% tmp$Package) {
      suppressPackageStartupMessages(library(pkg,
                                             character.only = TRUE))
    } else {
      print(sprintf("%s %s",
                    pkg,
                    "is not installed. Installing it!"))
      
      if (pkg %in% BiocManager::available(pkg)) {
        BiocManager::install(pkg,
                             dependencies = TRUE,
                             update = TRUE)
      } else {
        install.packages(pkg,
                         dependencies = TRUE,
                         ask = FALSE)
      }
    }
  }
}

# Load required packages or install them if necessary
dependencies <- c(
  "optparse",
  "ape",
  "phytools",
  "treeio",
  "tidytree",
  "TreeTools",
  "ggstar",
  "ggtree",
  "dplyr",
  "plotly"
)


.load_packages(dependencies)

#==================================== TREE MANIPULATION FUNCTIONS ====================================#

#' read_tree
#' Read a phylogenetic tree from a file.
#'
#' @param input_file Path to the file containing the phylogenetic tree.
#' 
#' @return An object representing the phylogenetic tree.
#' 
#' @export

# Function to read and preprocess a tree
read_tree <- function(input_file, bootstrap_support = TRUE) {
  if (bootstrap_support == TRUE) {
    t <- treeio::read.newick(input_file, node.label = 'support')
  } else {
    t <- treeio::read.newick(input_file)
  }
  
  return(t)
}

#' .load_tree_object
#' Load phylogenetic tree file and/or object as treedata object
#' 
#' @param tree An object representing the phylogenetic tree. Should be a newick or nexus file, or an object of class 'treedata' or 'phylo'.
 
.load_tree_object <- function(tree) {
  if (typeof(tree) == "character") {
    if (file.exists(tree)) {
      if (any(grepl(".newick|.nwk|.tre|.support|.nxs|.nex", tree))) {
        tree_obj <- read_tree(tree)
      } else {
        stop("Provided file is not a newick or nexus file.")
      }
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
#' 
#' @export

write_tree <- function(tree, output = "output_tree.nwk") {
  tree_obj <- .load_tree_object(tree)

  # Copy bootstrap values in the labels of the internal (non-Tip) nodes
  phylo_data <- as.treedata(
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


#' print_internal_nodes
#' This function generates a plot of the provided phylogenetic tree with node IDs displayed. Additionally, it offers the option
#' to highlight specific tip labels that match a user-provided pattern.
#'
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo'.
#' @param form Layout of the tree based on ggtree options. Layout can be rectangular, circular, roundrect, slanted, ellipse, fan, equal_angle, daylight (Default: "rectangular").
#' @param node_id_color Color of the printed node ids. Default is "darkred".
#' @param ... Pattern(s) for printing tip labels that match them.
#' 
#' @return A ggplot object representing the phylogenetic tree with node IDs and optional highlighted tip labels.
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
 
print_internal_nodes <- function(tree, form = "circular", node_id_color = "darkred", tip_label_size = 2, ...) {
  tree_obj <- .load_tree_object(tree)

  # Create the base tree plot with node labels
  plot <- ggtree(tree_obj, layout = form) +
    geom_nodelab(aes(label = node), hjust = -0.1, color = node_id_color, size = 3)

  references <- list(...)

  if (length(references) == 0) {
    print("Only node IDs will be printed.")
    return(plot)
  } else {
    reference <- as.character(unlist(references))
    matching_labels <- tree_obj@phylo$tip.label[grepl(reference, tree_obj@phylo$tip.label)]

    matching_dict <- data.frame(
      label = matching_labels,
      matching_flag = TRUE
    )

    plot <- plot %<+% matching_dict
    plot <- plot + geom_tiplab(aes(label = ifelse(matching_flag, label, ""), fill = "black"), size = tip_label_size, align.tip.label = TRUE)
    return(plot)
  }
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

#' flip_nodes
#'
#' Flip nodes on a phylogenetic tree based on the specified descendant nodes. It can flip internal nodes or leaves, if the node is terminal.
#'
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo'.
#' @param node1 The first node for flipping.
#' @param node2 The second node for flipping.
#'
#' @return A 'phylo' object with flipped nodes.
#'
#' @examples
#' \dontrun{
#' # Load a sample tree
#' tree <- read_tree("path/to/tree.nwk")
#'
#' # Use print_internal_nodes to identify the nodes you want to flip, ideally with a pattern matching a leaf of interest.
#' print_internal_nodes(tree)
#'
#' # Flip nodes 5 and 8 in the tree
#' flipped_tree <- flip_node(tree, 5, 8)
#' }
#'
#' @export

flip_nodes <- function(tree, node1,  node2) {
  tree_obj <- .load_tree_object(tree)
  return( ggtree::flip(ggtree(tree_obj), node1, node2) )
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
#' Function to extract a subtree by finding the MRCA (Most Recent Common Ancestor) of two anchor nodes while preserving branch lengths and bootstrap values.
#'
#' @param tree A phylogenetic tree object of class 'phylo' or 'treedata'.
#' @param tip1 The label or pattern of the first tip node.
#' @param tip2 The label or pattern of the second tip node.
#' @param branch_length Logical. If TRUE, preserve branch lengths in the extracted subtree; if FALSE, create a cladogram with equal branch lengths.
#'
#' @return A treedata object representing the extracted subtree with preserved branch lengths and bootstrap values.
#'
#' @examples
#' # Extract a subtree preserving branch lengths and bootstrap values
#' tree <- read.tree("path/to/your/treefile.newick")
#' subtree <- extract_subtree(tree, "taxon_A", "taxon_B")
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

extract_subtree <- function(tree, tip1, tip2) {
  tree_obj <- .load_tree_object(tree)

  if ("support" %in% colnames(tree_obj@data)) {
    # Append bootstrap values in the label section, in the nodes which are NOT tips
    # The bootstrap values will be kept there and transferred along with the labels in the new node numbering of the subtree
    phylo_data <- as.treedata(
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

  t <- Preorder(phylo_data@phylo)

  # Extract subtree based on anchor tip labels and/or patterns
  tip1m <- t$tip.label[grepl(tip1, t$tip.label)]
  tip2m <- t$tip.label[grepl(tip2, t$tip.label)]

  if (any(sapply(list(tip1m, tip2m), function(x) is.null(x)))) {
    stop("Provided tip labels are NULL or were not found among tree tip labels.")
  } else if (any(sapply(list(tip1m, tip2m), function(x) length(x) > 1))) {
    stop("Provided tip names are not unique.")
  }

  subtree <- ape::extract.clade(t, node = phytools::findMRCA(t, c(tip1m, tip2m)))

  if (!("support" %in% colnames(tree_obj@data))) {
    return(subtree)
  } else {
    # If node is NOT tip, add the labels in the bs_support column and make the labels NA - mapping according to new node numbering
    subtree_f <- as.treedata(subtree) %>%
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

#==================================== TREE VISUALIZATION FUNCTIONS ====================================#

#' highlight_tree
#'
#' Function to highlight nodes and descendant leaves on a phylogenetic tree.
#'
#' @param tree A phylogenetic tree object of class 'phylo' or a file path to a tree file (e.g., in Newick format).
#' @param highlight_nodes A associative vector of node labels with associated group names to highlight on the tree.
#' @param colors An associative vector of colors with optional group names to connect with the highlight nodes. If not provided, random colors will be assigned.
#' @param layout A character string specifying the layout of the tree. Options include "rectangular", "circular", etc.
#' @param name A character string specifying the name for the highlighted plot.
#'
#' @return A ggtree object with highlighted nodes.
#'
#' @examples
#' # Highlight specific nodes with random colors
#' tree <- read.tree("path/to/your/treefile.newick")
#' highlight_tree(tree, c("taxon_A", "taxon_B"))
#'
#' # Highlight specific nodes with specified colors
#' tree <- read.tree("path/to/your/treefile.newick")
#' highlight_tree(tree, 
#'                highlight_nodes = c("Category_A" = node_number_A, 
#'                                    "Category_B" = node_number_B), 
#'                colors = c("Category_A" = color_A, 
#'                           "Category_B" = color_B)
#'                )
#'
#' @seealso
#' \code{\link{visualize_tree}} for visualizing phylogenetic trees.
#'
#' @importFrom ggtree ggtree geom_hilight  
#' @importFrom ggstar geom_star
#' @importFrom ggplot2 scale_fill_manual ggsave
#' @importFrom treeio read.newick
#'
#' @export

# Function to highlight nodes on a tree
highlight_tree <- function(tree, form = "circular", highlight_nodes, colors = NULL, name = NULL, ...) {
  tree_obj <- .load_tree_object(tree)

  plot <- ggtree(tree_obj, layout = form) +
    geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) +
    scale_starshape_identity() + 
    scale_fill_identity()

  # Check if highlight_nodes is a list or a vector
  if (is.vector(highlight_nodes)) {
    # It's a vector, create a data frame with provided colors or generate random colors
    highlight <- data.frame(
      Group = names(highlight_nodes),
      Label = highlight_nodes,
      Color = sapply(names(highlight_nodes), function(x) {
        if (is.null(colors)) {
          sample(colors(), 1)
        } else {
          colors[x]
        }
      })
    )
  } else {
    stop("highlight_nodes should a vector.")
  }

  # Highlight the nodes
  plot <- plot + geom_hilight(
    data = highlight,
    aes(node = Label, fill = Color),
    color = "black",
    alpha = 0.3,
    extend = 0.10,
    linetype = 1,
    linewidth = 0.9,
    show.legend = TRUE
  )
  return(plot)
}


#' visualize_tree
#' Function to visualize a tree with a wide variety of options and customizable features
#'
#' @param tree Phylogenetic tree in Newick format, treedata, or phylo class object.
#' @param form Layout of the tree based on ggtree options. Layout can be rectangular, circular, roundrect, slanted, ellipse, fan, equal_angle, daylight (Default: "rectangular").
#' @param tiplabels Logical. Print tip labels on the tree.
#' @param pattern_id Pattern to match tip labels for printing. If tip_label_colors is provided, the tip labels matching that pattern will be printed with specified color. Default: black.
#' @param bootstrap_numbers Logical. Display bootstrap values on branches. Default: TRUE
#' @param bootstrap_number_nudge_y Numeric. Controls the relative height of the bootstrap number on top of the branch.
#' @param bootstrap_circles Logical. Display bootstrap values as circles on parent nodes of branches.
#' @param bootstrap_legend Logical. Display legend for bootstrap circles.
#' @param color Vector specifying tip color mappings.
#' @param shape Vector specifying tip shape mappings.
#' @param mappings_legend Logical. Display legend for color and shape mappings.
#' @param tip_label_size Numeric. Size of tip labels.
#' @param tip_shape_size Numeric. Size of tip shapes.
#' @param tip_label_colors Associative vector ("taxon_name" = "color") or single string ("color") combined with pattern_id.
#' @param reference Vector or Dataframe with patterns or geneids to be printed as tip labels, when color and/or shape mappings are provided.
#' @param clades Associative vector of node IDs and related labels.
#' @param save Logical. Save the plot.
#' @param output Character. Output file name if saving the plot.
#' @param interactive Logical. Generate an interactive plot using plotly.
#' @param ... Add pattern(s) for reference taxon IDs. Works only if color and/or shape mappings are provided
#'
#' @return A ggplot or ggplotly object representing the visualized tree.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(ggplot2)
#' library(ggtree)
#'
#' # Example: Visualize a tree with bootstrap values as circles and custom colors
#' visualize_tree(tree, 
#'                form = "rectangular", 
#'                bootstrap_circles = TRUE,
#'                color = c( taxon1 = "red", 
#'                           taxon2 = "blue",
#'                           taxon3 = "orange" ),
#'                shape = c( taxon1 = "circle", 
#'                           taxon2 = "star",
#'                           taxon3 = "circle" ),
#'                references = c( "Taxon1",
#'                                "Taxon2" ) 
#'                )
#' }
#'
#' @seealso \code{\link[ggtree]{ggtree}}, \code{\link[ggtree]{geom_tiplab}},
#' \code{\link[ggtree]{geom_nodelab}}, \code{\link[ggtree]{geom_cladelab}},
#' \code{\link[ggstar]{geom_star}}, \code{\link[ggplot2]{theme_tree}}, \code{\link[ggplot2]{ggsave}},
#' \code{\link[plotly]{ggplotly}}
#' 
#' @importFrom ggplot2 aes theme ggsave
#' @importFrom ggtree ggtree geom_tiplab geom_nodelab geom_cladelab theme_tree 
#' @importFrom ggstar geom_star
#' @importFrom plotly ggplotly
#'
#' @export 


visualize_tree <- function(tree = tree,
                           form = "rectangular",
                           tiplabels = FALSE,
                           pattern_id = NULL,
                           color = NULL,
                           shape = NULL,
                           mappings_legend = FALSE,
                           tip_shape_size = 3,
                           references = NULL,
                           tip_label_size = 2,
                           tip_label_colors = NULL, 
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
                           tree_width_limits = NULL,
                           save = TRUE,
                           output = NULL,
                           interactive = FALSE) {


  tree_obj <- .load_tree_object(tree)
  
  # If option for branch length is not TRUE, make branch lengths NULL
  if (!(branch_length == TRUE | branch_length == T)) {
    print("Option branch length is deactivated. As a result, the tree will be visualized as a cladogram with equal branch lengths.")
    tree_obj@phylo$edge.length <- NULL
  } 
  
  # Check if bootstrap support is present in the tree object
  if (any(!is.null(tree_obj@data$support)) && any(!is.na(tree_obj@data$support))) {
      if ( all(tree_obj@data$support < 10, na.rm = TRUE) == TRUE) { # Bootstrap values are in scale of 1. 
        warning("Provided bootstrap values were in scale of 1 and were multiplied by 100 to change scales.")
        
        tree_obj@data$support <- sapply(tree_obj@data$support, function(x) {
            x <- x * 100
        })
      }
    
   tree_obj <- as.treedata ( as_tibble(tree_obj) %>%
                             mutate (bs_color = case_when(support < 50 ~ "snow2",
                                                          support >= 50 & support < 75 ~ "grey",
                                                          support >= 75 ~ "black",
                                                          TRUE ~ NA  # Default case if none of the conditions are met
                                                         )))
      
   plot <- ggtree(tree_obj, layout = form)
   
   # Adjust tree line width based on the number of tips -- larger trees require smaller width
   if ( !is.null(tree_width_limits) ) {
     if ( length(tree_obj@phylo$tip.label) > min(tree_width_limits) ) {
       plot <- plot + geom_treescale(width = 1)
     } else if ( length(tree_obj@phylo$tip.label) > min(tree_width_limits) && max(tree_width_limits) > length(tree_obj@phylo$tip.label) ) {
       plot <- plot + geom_treescale(width = 0.5)
     } else if (  length(tree_obj@phylo$tip.label) > max(tree_width_limits) ) {
       plot <- plot + geom_treescale(width = 0.3)
     } 
   } 
   
   # Print bootstrap values either as numbers on top of branches or as circles in the parent nodes
   if ( bootstrap_numbers == TRUE && bootstrap_circles != TRUE ) {
    if (any(!is.null(tree_obj@data$support)) && any(!is.na(tree_obj@data$support))) {
      # bootstrap_number_nudge_y controls the relative height of the number on top of branch - use LARGER values as the SIZE of the tree INCREASES
      plot <- plot + geom_nodelab(  mapping = aes(x = branch, label = support), nudge_x = bootstrap_number_nudge_x, nudge_y = bootstrap_number_nudge_y, size = node_label_size) 
      } else {
          print ("Bootstrap values do not exist.")
      }
    } else if ( bootstrap_numbers != TRUE && bootstrap_circles == TRUE ) {
        plot <- plot + geom_point2(fill = ifelse(!is.na(tree_obj@data$bs_color), tree_obj@data$bs_color, NA),
                                   color = ifelse(!is.na(tree_obj@data$bs_color), "black", NA),
                                   shape = 21, size = bootstrap_circle_size)
          
        if (bootstrap_legend == TRUE) {
          plot <- plot + theme_tree(legend.position = c(0.2, 0.2))
        }
            
    } else if ( bootstrap_numbers == TRUE && bootstrap_circles == TRUE ) {
          stop ("bootstrap_numbers = TRUE (default) and bootstrap_circles = TRUE. Please select either bootstrap circles or bootstrap numbers option.")
    } else {
        print ("Both bootstrap_circles and bootstrap_numbers options are FALSE. Bootstrap values will not be printed.")
    }
  } else {
      tree_obj <- tree_obj
      plot <- ggtree(tree_obj, layout = form)
      print ("Bootstrap values are not present in the tree object.")
  }
  
   # Visualize phylogenetic tree with an option to select the printed tips, based on pattern
   if (is.null(color) && is.null(shape) ) {
     if (length(references) > 0) { 
       stop("The 'references' option can be used only if tip mappings are provided. If you want to print tip labels with a specific pattern use the tiplabels = TRUE and pattern = 'desired pattern' options")
     }
     
   # If you want to print some of the taxon labels as text, do not include them into the color and shape mapping vectors, and provide them as 'references'
   if (tiplabels == TRUE) { 
     if (is.null(pattern_id)) {
         plot <- plot + geom_tiplab(size = tip_label_size, color = ifelse(is.null(tip_label_colors), "black", tip_label_colors), 
                                    show.legend = mappings_legend, align.tip.label = TRUE)
         
         print("Printing phylogenetic tree with all tip labels!")
       } else {
          if ( !is.null(tip_label_colors) ) {
            if ( length(tip_label_colors) != length(pattern_id) ) {
             stop("tip_label_colors and pattern_ids should have the same length!")
            } 
          }
         
         if (any(sapply(pattern_id, function(pattern) any(grepl(pattern, tree_obj@phylo$tip.label))))) {
           matching_labels <- unlist(sapply(pattern_id, function(pattern) {
             grep(pattern, tree_obj@phylo$tip.label, fixed = TRUE, value = TRUE)
           }), use.names = FALSE)
           
             matching_dict <- data.frame(
               label = matching_labels,
               matching_flag = TRUE,
               color_l = if (is.null(tip_label_colors)) {
                            "black"
                         } else if (length(tip_label_colors) > 0) {
                            unlist( lapply(names(group_colors), function(label) { 
                              ifelse( grepl(label, matching_labels), tip_label_colors[label], NA )
                            }))
                          } else {
                              print(paste("No color was found in the associative tip_label_colors vector for", x, "!\n"))
                              print(paste("Will use", tip_label_colors, "for all", x, "instances!"))
                              tip_label_colors
                          }
              )
      
         matching_dict <- matching_dict[complete.cases(matching_dict), ] # Remove lines with NA values
         
         plot <- plot %<+% matching_dict
         plot <- plot + geom_tiplab(aes(label = ifelse(matching_flag == TRUE, label, NA), 
                                         color = ifelse(matching_flag == TRUE, color_l, NA) ),
                                         size = tip_label_size, show.legend = mappings_legend, align.tip.label = TRUE) +
                             scale_color_identity()
              
              cat("Printing phylogenetic tree plot with tip labels matching", pattern_id, "!\n")
          } else {
              cat ("Provided", pattern_id, "was not found among tip labels! Printing phylogenetic tree plot without tip labels!")
              plot <- plot
          }
       }
   } else { 
      if (!is.null(pattern_id) || !is.null(tip_label_colors)) {
            print ("Please enable tip label printing with tiplabels = TRUE.")
      } else {
            print ("Printing phylogenetic tree plot without tip labels!")
      }
   }
  } else { # Visualize phylogenetic tree with color and shape mappings, and an option to print reference taxon as text

         taxon_dict <- list()
    
         if ( is.null(references) ) { 
            taxon_dict <- lapply(seq_along(tree_obj@phylo$tip.label), function(i) {
             list(tip_label = tree_obj@phylo$tip.label[i], taxon = sub("_.*", "", tree_obj@phylo$tip.label[i]))
           })
         } else if (any(sapply(references, function(pattern) any(grepl(pattern, tree_obj@phylo$tip.label))))) { # Manipulate tip labels to generate taxon_names and isolate tip labels which match reference
                  
                   matching_labels <- c()
                
                   matching_labels <- unlist(lapply(references, function(reference) {
                       grep(reference, tree_obj@phylo$tip.label, fixed = TRUE, value = TRUE)
                   }))
                   
                   ref_taxons <- data.frame( label = matching_labels, reference_flag = rep(TRUE, length(matching_labels)) )
                   
                   plot_ref <- plot %<+% ref_taxons 
                   plot <- plot_ref + geom_tiplab(aes(label = ifelse(reference_flag, label, "")), 
                                                  size = tip_label_size, color = "black", align.tip.label = TRUE)
                     
                   non_ref_taxons <- setdiff(tree_obj@phylo$tip.label, ref_taxons)
                   
                   taxon_dict <- lapply(non_ref_taxons, function(label) {
                     list(tip_label = label, taxon = sub("_.*", "", label))
                   })
                } else {
                  print ("Provided 'references' pattern was not found.")
                  
                  taxon_dict <- lapply(tree_obj@phylo$tip.label, function(label) {
                    list(tip_label = label, taxon = sub("_.*", "", label))
                  })
          }
  
         if (length(taxon_dict) > 0) {
           
           tip_labels <- sapply(taxon_dict, function(entry) entry$tip_label)
           taxon_names <- sapply(taxon_dict, function(entry) entry$taxon)
           
           mappings_legend = TRUE
           
           if (!is.null(color) && !is.null(shape)) {
             
             sapply(references, function(reference) { 
               if ( any(grepl(reference, names(color))) || any(grepl (reference, names(shape))) ) {
                 cat(reference, "was found among the colors and/or shapes names, but it will be printed as a 'reference' tip label!\n")
               }
             })
             
             missing_taxon_color <- setdiff(unique(sapply(taxon_dict, function(entry) entry$taxon)), names(color))
             
             if (length(missing_taxon_color) > 0) {
               warning(paste("Color mapping not found for the following taxon: ", paste(missing_taxon_color, collapse = ", ")))
             }
              
             # Create data frames for tip colors and shapes
             tip_colors_df <- data.frame(
               label = tip_labels,
               taxon_colors = taxon_names,
               s_color = color[match(sapply(taxon_dict, function(entry) entry$taxon), names(color))]
             )
              
             # Check for missing taxon in color dataframe
             missing_taxon_shape <- setdiff(unique(sapply(taxon_dict, function(entry) entry$taxon)), names(shape))
             if (length(missing_taxon_shape) > 0) {
                warning(paste("Shape mapping not found for the following taxon: ", paste(missing_taxon_shape, collapse = ", ")))
             }
               
              tip_shapes_df <- data.frame(
                label = tip_labels,
                taxon_shapes = taxon_names,
                s_shape = shape[match(sapply(taxon_dict, function(entry) entry$taxon), names(shape))]
              )
              
              # Integrate color and shape mappings into the ggtree object
              plot <- plot %<+% tip_colors_df 
              plot <- plot %<+% tip_shapes_df  
  
              plot <- plot + geom_star(mapping = aes( subset = isTip,
                                                      fill = ifelse(!is.na(taxon_colors), s_color, NA),
                                                      starshape = ifelse(!is.na(taxon_shapes), s_shape, NA)), 
                                                      show.legend = mappings_legend, size = tip_shape_size) + 
                             scale_fill_identity() +
                             scale_starshape_identity()
              
              if (tiplabels == TRUE) {
                plot <- plot + geom_tiplab(mapping = aes( subset = isTip & !is.na(s_color), color = ifelse(!is.na(s_color), s_color, NA)),
                                           show.legend = FALSE, size = tip_shape_size, align.tip.label = TRUE) +
                  scale_fill_identity()
              }
            } else if ( !is.null(color) && is.null(shape)) {
              
                tip_colors_df <- data.frame(
                  label = tip_labels,
                  taxon_colors = taxon_names,
                  s_color = color[match(sapply(taxon_dict, function(entry) entry$taxon), names(color))]
                  )
                
                plot <- plot %<+% tip_colors_df
                
                if (tiplabels == TRUE) {
                  plot <- plot + geom_tiplab(mapping = aes( subset = isTip & !is.na(s_color), color = ifelse(!is.na(s_color), s_color, NA)),
                                             show.legend = FALSE, size = tip_shape_size, align.tip.label = TRUE) +
                                 scale_fill_identity()
                } else {
                    plot <- plot + geom_star(mapping = aes( subset = isTip & !is.na(s_color),
                                                            fill = ifelse(!is.na(taxon_colors), s_color, NA)),
                                                            show.legend = mappings_legend, starshape = "circle", size = tip_shape_size) +
                                   scale_fill_identity()
                }
                
            } else if ( is.null(color) && !is.null(shape)) {
	    
                tip_shapes_df <- data.frame(
                  label = tip_labels,
                  taxon_shapes = taxon_names,
                  s_shape = shape[match(sapply(taxon_dict, function(entry) entry$taxon), names(shape))]
                  )
                
                plot <- plot %<+% tip_shapes_df
                
                if (tiplabels == TRUE) {
                  plot <- plot + geom_tiplab(mapping = aes( subset = isTip & !is.na(s_color), color = ifelse(!is.na(s_color), s_color, NA)),
                                             show.legend = FALSE, size = tip_shape_size, align.tip.label = TRUE) +
                    scale_fill_identity()
                } else {
                  plot <- plot + geom_star(mapping = aes( subset = isTip & !is.na(s_shape),
                                                        starshape = ifelse(!is.na(taxon_shapes), s_shape, NA)),
                                         fill = "black", size = tip_shape_size, show.legend = mappings_legend) +
                                 scale_starshape_identity()
                }
            }
         } else {
            stop("Provided tip mappings were not found.")
        }
      }
    
   # Option to denote specific clades in tree and include associated labels
   if (!is.null(clades)) {
       for (i in names(clades)) {
         plot <- plot + geom_cladelab(node = clades[i], label = i, align = TRUE, fill = 'black',
                                      offset.text = labeldist, barsize = 0.9, offset.bar = bardist, clade_label_size = clade_label_size)
       }
   }
  
  if (!(is.null(output))) {
    save == TRUE
  }
  
   # Export plot with provided options
   if (exists("plot")) {
     if (save == TRUE) {
       if (is.null(output)) {
         if (typeof(tree) == "character") {
           if (file.exists(tree)) {
             if (any(grepl(".newick|.nwk|.tre|.support|.nxs|.nex", tree))) {
                ggsave(plot = plot, sprintf("%s_visualized.svg", sub(".newick|.nwk|.tre|.support|.nxs|.nex", "", tree)), dpi = 600)
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
    
    if (interactive == FALSE) {
        return (plot)
    } else {
        return (plotly::ggplotly(plot))
    }
  }
}