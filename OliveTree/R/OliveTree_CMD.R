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
.load_packages <- function( tools ) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built, 1, 1)))
  tmp <- tmp[as.numeric(substr(tmp$Built, 1, 1)) == max_version, ]

  for ( pkg in tools ) {
    if ( pkg %in% tmp$Package ) {
      suppressPackageStartupMessages( library (pkg, character.only = TRUE) )
    } else {
         print(sprintf("%s %s", pkg, "is not installed. Installing it!"))
      
         if ( pkg %in% BiocManager::available(pkg) ) {
            BiocManager::install(pkg, dependencies = TRUE, update = TRUE)
         } else {
            install.packages(pkg, dependencies = TRUE, ask = FALSE)
      }
    }
  }
}

# Load required packages or install them if necessary
dependencies <- c( "optparse", "ape", "phytools", "treeio", "tidytree", 
                   "TreeTools", "ggstar", "ggtree","dplyr", "plotly" )

.load_packages(dependencies)

#==================================== TREE MANIPULATION FUNCTIONS ====================================#

#' read_tree
#' Read a phylogenetic tree from a file.
#'
#' @param input_file Path to the file containing the phylogenetic tree.
#' 
#' @return An object representing the phylogenetic tree.
#' 
#' @examples
#' \dontrun{
#' tree <- read_tree("path/to/tree.nwk")
#' }
#' 
#' @export

# Function to read and preprocess a tree
read_tree <- function(input_file, bootstrap_support = TRUE) {
  if (bootstrap_support == TRUE) {
    t <- treeio::read.newick(input_file, node.label = 'support' )
  } else {
      t <- treeio::read.newick(input_file)
  } 
  
  return(t)
}

#' .load_tree_object
#' Load phylogenetic tree file and/or object as treedata object
#' 
#' @param tree An object representing the phylogenetic tree. Should be a newick or nexus file, or an object of class 'treedata' or 'phylo'.
 
.load_tree_object <- function (tree) {
  if (typeof(tree) == "character") {
    if (file.exists(tree)) {
      if (any(grepl(".newick|.nwk|.tre|.support|.nxs|.nex", tree))) {
        tree_obj <- read_tree(tree)
      } else {
          stop ("Provided file is not a newick or nexus file.")
      }
    } 
  } else if (is(tree,"treedata")) {
      tree_obj <- tree
  } else if (is(tree,"phylo")) {
      tree_obj <- treeio::as.treedata(tree)
  } else {
      stop ("Provided file is not a treedata, a phylo or a tibble_df object.")
  }
  return(tree_obj)
}


#' node_ids
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
#' node_ids(tree, color = "darkred", "taxon_A")
#' }
#' 
#' @export
 
node_ids <- function(tree, form = "circular", node_id_color = "darkred", tip_label_size = 2, ...) {
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
      plot <- plot + geom_tiplab(aes(label = ifelse(matching_flag, label, ""), fill = "black"), size = tip_label_size)
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
#' # Use node_ids to identify the nodes you want to flip, ideally with a pattern matching a leaf of interest.
#' node_ids(tree)
#'
#' # Flip nodes 5 and 8 in the tree
#' flipped_tree <- flip_node(tree, 5, 8)
#' }
#'
#' @export

flip_nodes <- function(tree, node1, node2) {
  tree_obj <- .load_tree_object(tree)
  return(as.phylo(ggtree::flip(ggtree(tree_obj), node1, node2)))
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
#' subtree <- extract_subtree(tree, "species_A", "species_B", branch_length = TRUE)
#' visualize_tree(subtree)
#'
#' # Extract a cladogram with equal branch lengths
#' tree <- read.tree("path/to/your/treefile.newick")
#' cladogram <- extract_subtree(tree, "species_A", "species_B", branch_length = FALSE)
#' visualize_tree(cladogram)
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
  
  if ( "support" %in% colnames(tree_obj@data) ) { 
    
    # Append bootstrap values in the label section, in the nodes which are NOT tips
    # The bootstrap values will be kept there and transferred along with the labels in the new node numbering of the subtree

    phylo_data <- as.treedata(left_join(as_tibble(tree_obj@phylo), as_tibble(tree_obj@data)) %>%
                              mutate(isTip = ifelse(!is.na(label), TRUE, FALSE)) %>% 
                              mutate(label = ifelse(isTip == FALSE, as.character(support), label)) %>% 
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
  } else if (any(sapply(list(tip1m, tip2m), function(x) length(x) > 1 ))) {
    stop("Provided tip names are not unique.")
  }
    
  subtree <- ape::extract.clade(t, node = phytools::findMRCA(t, c(tip1m, tip2m)))
  
  if ( !( "support" %in% colnames(tree_obj@data)) ) { 
    return ( subtree )
  } else {
      # If node is NOT tip, add the labels in the bs_support column and make the labels NA - mapping according to new node numbering 
      subtree_f <- as.treedata(subtree) %>%
                   mutate(bs_support = ifelse(isTip == FALSE, label, NA)) %>%
                   mutate(label = ifelse(isTip == FALSE, NA, label))
      
      # To create an object compatible with visualize_tree, initialize the 'node' and 'support' columns in the data slice of the treedata class
      subtree_f@data <- tibble(
          node = rep(NA, length = length(subtree_f@extraInfo$node)),
          support = rep(NA, length = length(subtree_f@extraInfo$bs_support))
      )
      
      # ... and the add the node number and bs_support values
      subtree_f@data <- tibble(
          node = subtree_f@extraInfo$node,
          support = as.numeric(subtree_f@extraInfo$bs_support)
      ) 
      return ( subtree_f )
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
#' highlight_tree(tree, c("species_A", "species_B"))
#'
#' # Highlight specific nodes with specified colors
#' tree <- read.tree("path/to/your/treefile.newick")
#' highlight_tree(tree, c("species_A", "species_B"), colors = c("red", "blue"))
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
highlight_tree <- function(tree, highlight_nodes, colors = NULL, layout = "circular", name = NULL, ...) {
  tree_obj <- .load_tree_object(tree)
  
  plot <- ggtree(tree_obj, layout = layout) + 
          geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) + 
          scale_starshape_identity() + scale_fill_identity()

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
  plot <- plot + geom_hilight(data = highlight,
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
#' @param pattern_id Pattern to match tip labels for printing.
#' @param bootstrap_numbers Logical. Display bootstrap values on branches. Default: TRUE
#' @param bootstrap_number_nudge_y Numeric. Controls the relative height of the bootstrap number on top of the branch.
#' @param bootstrap_circles Logical. Display bootstrap values as circles on parent nodes of branches.
#' @param bootstrap_legend Logical. Display legend for bootstrap circles.
#' @param color Vector specifying tip color mappings.
#' @param shape Vector specifying tip shape mappings.
#' @param mappings_legend Logical. Display legend for color and shape mappings.
#' @param tip_label_size Numeric. Size of tip labels.
#' @param tip_shape_size Numeric. Size of tip shapes.
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
#' visualize_tree(tree, form = "rectangular", bootstrap_circles = TRUE,
#'                color = c(species1 = "red", species2 = "blue"))
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

visualize_tree <- function(tree, form = "rectangular", tiplabels = FALSE, pattern_id = NULL,
                           color = NULL, shape = NULL, mappings_legend = FALSE, tip_shape_size = 3,  
                           tip_label_size = 2, fontsize = 3,
                           bootstrap_numbers = TRUE, bootstrap_number_nudge_x = 1, bootstrap_number_nudge_y = 3.2, node_label_size = 3, bootstrap_circles = FALSE, bootstrap_legend = FALSE, 
                           branch_length = TRUE, clades = NULL, labeldist = 1, bardist = 1,
                           save = TRUE, output = NULL, interactive = FALSE, ...) {
  
  tree_obj <- .load_tree_object(tree)
  
  # If option for branch length is not TRUE, make branch lengths NULL
  if (!(branch_length == TRUE | branch_length == T)) {
    print("Option branch length is deactivated. As a result, the subtree will be extracted as a cladogram with equal branch lengths.")
    tree_obj@phylo$edge.length <- NULL
  } 
  
  # Check if bootstrap support is present in the tree object
  if (any(!is.null(tree_obj@data$support)) && any(!is.na(tree_obj@data$support))) {
      tree_obj <- as.treedata (as_tibble(tree_obj) %>%
                  mutate(bs_color = case_when(support < 50 ~ "snow2",
                                              support >= 50 & support < 75 ~ "grey",
                                              support >= 75 ~ "black",
                                              TRUE ~ NA  # Default case if none of the conditions are met
                                              )))
      
      plot <- ggtree(tree_obj, layout = form)
  
      # Manipulate bootstrap values
      if ( bootstrap_numbers == TRUE && bootstrap_circles != TRUE ) {
        if (any(!is.null(tree_obj@data$support)) && any(!is.na(tree_obj@data$support))) {
          # bootstrap_number_nudge_y controls the relative height of the number on top of branch
          # Use larger values as the size of the tree increases
          plot <- plot + geom_nodelab(  mapping = aes(x = branch, label = support), nudge_x = bootstrap_number_nudge_x, nudge_y = bootstrap_number_nudge_y, size = node_label_size) 
          } else {
              print ("Bootstrap values do not exist.")
          }
      } else if ( bootstrap_numbers != TRUE && bootstrap_circles == TRUE ) {
            plot <- plot + geom_point2(fill = ifelse(!is.na(tree_obj@data$bs_color), tree_obj@data$bs_color, NA),
                                         color = ifelse(!is.na(tree_obj@data$bs_color), "black", NA),
                                         shape = 21, size = 1.7)
          
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
     if (length(list(...)) > 0) { 
       stop("The 'reference' option can be used only if tip mappings are provided.")
     }
     
     # Tip labels == TRUE and tip mappings & reference options are incongruent. 
     # If you want to print some of the species labels as text, do not include them into the color and shape mapping vectors, and provide them as 'reference'
     if (tiplabels == TRUE) { 
       if (is.null(pattern_id)) {
         plot <- plot + geom_tiplab(size = tip_label_size) # Here default size is 1, because it will print all tip labels
         print("Printing phylogenetic tree with all tip labels!")
       } else {
           if (any(grepl(pattern_id, tree@phylo$tip.label)) == TRUE) {
              matching_labels <- tree@phylo$tip.label[grepl(pattern_id, tree@phylo$tip.label)]
               
              matching_dict <- data.frame(
                label = matching_labels,
                matching_flag = TRUE
              )
               
             plot <- plot %<+% matching_dict
             plot <- plot + geom_tiplab(aes(label = ifelse(matching_flag, label, "")), color = "black", size = tip_label_size)
              
             cat("Printing phylogenetic tree plot with tip labels matching", pattern_id, "!\n")
                
           } else {
              if (is.null(pattern_id)) {
                print ("Please enable tip label printing with tiplabels = TRUE.")
              }
             
              cat ("Provided", pattern_id, "was not found among tip labels! Printing phylogenetic tree plot without tip labels!")
              plot <- plot
            }
          } 
     } else {
        print ("Printing phylogenetic tree plot without tip labels!")
     }
   } else { # Visualize phylogenetic tree with color and shape mappings, and an option to print reference species as text
       # Manipulate tip labels to generate species_names and isolate tip labels which match reference
     
      if (tiplabels == TRUE) { 
        stop("The tip labels option can be used only if tip mappings are NOT provided.\n  Add pattern you want print as text using the 'reference' option.")
      }
       
       references <- list(...)
       reference_args <- names(references)
       
       if (length(reference_args) == 0) {
         species_dict <- lapply(seq_along(tree_obj@phylo$tip.label), function(i) {
           list(tip_label = tree_obj@phylo$tip.label[i], species = sub("_.*", "", tree_obj@phylo$tip.label[i]))
         })
       } else {
           reference <- as.character(unlist(references))
           
           if (any(grepl(reference, tree_obj@phylo$tip.label))) {
             
               ref_species <- lapply(references, function(reference) {
                 matching_labels <- tree_obj@phylo$tip.label[which(sapply(tree_obj@phylo$tip.label, function(x) any(grepl(reference, x))))]
                 data.frame(label = matching_labels, reference_flag = TRUE)
               })
               
               ref_dict <- do.call(rbind, ref_species)

               plot_ref <- plot %<+% ref_dict 
               plot <- plot_ref + geom_tiplab(aes(label = ifelse(reference_flag, label, "")), size = tip_label_size, color = "black")
               
               non_ref_species <- setdiff(tree_obj@phylo$tip.label, ref_species)
               species_dict <- lapply(non_ref_species, function(label) {
                 list(tip_label = label, species = sub("_.*", "", label))
               })
           } else {
             print ("Provided 'reference' pattern was not found.")
             
             species_dict <- lapply(tree_obj@phylo$tip.label, function(label) {
               list(tip_label = label, species = sub("_.*", "", label))
             })
           }
       }

       if (length(species_dict) > 0) {
         
         tip_labels <- sapply(species_dict, function(entry) entry$tip_label)
         species_names <- sapply(species_dict, function(entry) entry$species)
         
         mappings_legend = TRUE

         if (!is.null(color) && !is.null(shape)) {
            
           missing_species_color <- setdiff(unique(sapply(species_dict, function(entry) entry$species)), names(color))
           if (length(missing_species_color) > 0) {
             warning(paste("Color mapping not found for the following species: ", paste(missing_species_color, collapse = ", ")))
           }
            
           # Create data frames for tip colors and shapes
           tip_colors_df <- data.frame(
             label = tip_labels,
             species_colors = species_names,
             s_color = color[match(sapply(species_dict, function(entry) entry$species), names(color))]
           )
            
           # Check for missing species in color dataframe
           missing_species_shape <- setdiff(unique(sapply(species_dict, function(entry) entry$species)), names(shape))
            if (length(missing_species_shape) > 0) {
              warning(paste("Shape mapping not found for the following species: ", paste(missing_species_shape, collapse = ", ")))
            }
             
            tip_shapes_df <- data.frame(
              label = tip_labels,
              species_shapes = species_names,
              s_shape = shape[match(sapply(species_dict, function(entry) entry$species), names(shape))]
            )
             
            # Join data frames to create a unique mapping tibble
            plot <- plot %<+% tip_colors_df 
            plot <- plot %<+% tip_shapes_df  
            
            plot <- plot + geom_star(mapping = aes( subset = isTip,
                                                    fill = ifelse(!is.na(species_colors), s_color, NA),
                                                    starshape = ifelse(!is.na(species_shapes), s_shape, NA)),
                                     size = tip_shape_size, show.legend = mappings_legend) + 
                           scale_fill_identity() +
                           scale_starshape_identity()
                           
          } else if ( !is.null(color) && is.null(shape)) {
              tip_colors_df <- data.frame(
                label = tip_labels,
                species_colors = species_names,
                s_color = color[match(sapply(species_dict, function(entry) entry$species), names(color))]
                )
              
              plot <- plot %<+% tip_colors_df
              plot <- plot + geom_star(mapping = aes( subset = isTip & !is.na(s_color),
                                                      fill = ifelse(!is.na(species_colors), s_color, NA)),
                                       starshape = "circle", 
                                       size = tip_shape_size, show.legend = mappings_legend) +
                             scale_fill_identity()
              
          } else if ( is.null(color) && !is.null(shape)) {
              tip_shapes_df <- data.frame(
                label = tip_labels,
                species_shapes = species_names,
                s_shape = shape[match(sapply(species_dict, function(entry) entry$species), names(shape))]
                )
              
              plot <- plot %<+% tip_shapes_df
              plot <- plot + geom_star(mapping = aes( subset = isTip & !is.na(s_shape),
                                                      starshape = ifelse(!is.na(species_shapes), s_shape, NA)),
                                       fill = "black",
                                       size = tip_shape_size, show.legend = mappings_legend) +
                             scale_starshape_identity()
          }
       } else {
          stop("Provided tip mappings were not found.")
       }
   }
        
   if (!is.null(clades)) {
       for (i in names(clades)) {
         plot <- plot + geom_cladelab(node = clades[i], label = i, align = TRUE, fill = 'black',
                                      offset.text = labeldist, barsize = 0.9, offset.bar = bardist, fontsize = fontsize)
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


#=================================== PARSE ARGUMENTS ====================================#

parse_arguments <- function() {
  option_list <- list(
    make_option(c("--input_tree"), type = "character", required = TRUE),
    make_option(c("--form"), type = "character", default = "rectangular"),
    make_option(c("--node_id_color"), type = "character", default = "darkred"),
    make_option(c("--tip_label_size"), type = "numeric", default = 2),
    make_option(c("--highlight_nodes"), type = "character", nargs = "*"),
    make_option(c("--bootstrap_collapse_cutoff"), type = "numeric", default = 0.5),
    make_option(c("--flip_nodes"), type = "character", nargs = "*"),
    make_option(c("--group_descendants"), type = "character", nargs = "*"),
    make_option(c("--extract_subtree_tips"), type = "character", nargs = 2),
    make_option(c("--highlight_tree_colors"), type = "character", nargs = "*"),
    make_option(c("--highlight_tree_layout"), type = "character", default = "circular"),
    make_option(c("--visualize_tree_tiplabels"), action = "store_true"),
    make_option(c("--visualize_tree_color"), type = "character", nargs = "*"),
    make_option(c("--visualize_tree_shape"), type = "character", nargs = "*"),
    make_option(c("--visualize_tree_interactive"), action = "store_true"),
    make_option(c("--output"), type = "character"),
    make_option(c("--branch_color"), type = "character"),
    make_option(c("--branch_width"), type = "numeric"),
    make_option(c("--node_size"), type = "numeric"),
    make_option(c("--legend"), action = "store_true"),
    make_option(c("--legend_title"), type = "character"),
    make_option(c("--legend_position"), type = "character"),
    make_option(c("--help"), action = "store_true")
  )
  
  parser <- OptionParser(option_list = option_list)
  args <- parse_args(parser)
  
  if (args$help) {
    print_help(parser)
    q("no", status = 0)
  }
  
  return(args)
}


# Function to execute the main workflow
main <- function() {
  # Parse command line arguments
  args <- parse_arguments()
  
  # Read the input tree file
  input_tree <- read.tree(args$input_tree)
  
  #==================================== TREE MANIPULATION FUNCTIONS ====================================#
  
  # Bootstrap Collapse
  if (!is.null(args$bootstrap_collapse_cutoff)) {
    input_tree <- bootstrap_collapse(input_tree, args$bootstrap_collapse_cutoff)
  }
  
  # Flip Nodes
  if (!is.null(args$flip_nodes)) {
    input_tree <- flip_nodes(input_tree, args$flip_nodes)
  }
  
  # Group Descendants
  if (!is.null(args$group_descendants)) {
    input_tree <- group_descendants(input_tree, args$group_descendants)
  }
  
  # Extract Subtree
  if (!is.null(args$extract_subtree_tips)) {
    input_tree <- extract_subtree(input_tree, args$extract_subtree_tips)
  }
  
  #==================================== TREE VISUALIZATION FUNCTIONS ====================================#
  
  # Highlight Tree
  if (!is.null(args$highlight_nodes)) {
    input_tree <- highlight_tree(input_tree, args$highlight_nodes, args$highlight_tree_colors, args$highlight_tree_layout)
  }
  
  # Visualize Tree
  visualize_tree(input_tree, args$form, args$visualize_tree_tiplabels, color = args$visualize_tree_color,
                 shape = args$visualize_tree_shape, output = args$output, interactive = args$visualize_tree_interactive)
}

