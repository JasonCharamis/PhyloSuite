
# OliveTRee: R package for advanced tree manipulation and visualization by combining and extending the functionality of ggtree, ape, phytools and other related tools."

#' @name .load_packages
#' @description 
#' This function checks if a package is installed. If not, it installs the package using BiocManager if available, otherwise using install.packages. Then it loads it.
#' 
#' @param tools A vector or list of package names to be checked and installed.
#' @return NULL
#' 
#' @export

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

bioconductor_version <- ifelse(getRversion() >= "4.4", "3.19", ifelse(getRversion() >= "4.3", "3.18", "3.17"))

#==================================== Reading, Loading and Exporting Trees ====================================#

#' @name read_tree 
#' @description 
#' Read a phylogenetic tree from a file with auto-detection of bootstrap values.
#'
#' @param input_file Path to the file containing the phylogenetic tree. Input tree can be in newick or nexus format.
#' 
#' @return A tree data object representing the phylogenetic tree.
#' @examples
#' Read a sample tree from newick file
#'  tree <- read_tree("path/to/tree.nwk")
#'
#' Read a sample tree from nexus file
#'  tree <- read_tree("path/to/tree.nexus")
#'
#' @export


read_tree <- function(input_file) {
  if (grepl("\\.(newick|nwk|tre|tree)$", tolower(input_file))) {
    t <- phytools::read.newick(input_file)
    print (t)
  } else if (grepl("\\.(nexus|nxs|nex)$", tolower(input_file))) {
      t <- ape::read.nexus(input_file)
  } else if (grepl("\\.(treefile)$", tolower(input_file))) {
      t <- treeio::read.iqtree(input_file)
      colnames(t@data)[colnames(t@data) == "UFboot"] <- "support"
      return(t)
  } else {
      stop("Unsupported file format.")
  }
  
  # Check if the tree contains bootstrap values, and if it contains load them
  if ( any(is.numeric(t$node.label) || is.integer(t$node.label))) {
    if (any(as.numeric(t$node.label) > 0) ) {
      to <- treeio::as.treedata(t, node.label = "support")
    } 
  } else {
      to <- treeio::as.treedata(t)
  }
  return(to)
}


#' @name .load_tree_object
#' @description 
#' Load phylogenetic tree as treedata object.
#' 
#' @param tree An object representing the phylogenetic tree, providing a single tree entrypoint used in the whole package. 
#'             Makes use of the read_tree function. 
#'             If no file path is provided, the user will be prompted to provide an input phylogenetic tree through Rstudio's graphical user interface (GUI).
#' 
#' @return A tree data object representing the phylogenetic tree.
#' 
#' @examples
#' Load a tree data object from newick file
#'   tree_obj <- .load_tree_obj("path/to/tree.nwk")
#'  
#' tree: 'phylo' or 'treedata' object or 'path/to/tree.nwk' or 'path/to/tree.nexus' files
#'   tree_obj <- .load_tree_obj(tree)
#' 
#' @export


.load_tree_object <- function(tree = NULL) {
  
  if ( is.null(tree) ) {
    tree <- read_tree(file.choose())
  } else if (typeof(tree) == "character") {
    if (file.exists(tree)) {
      if (any(grepl(".newick|.nwk|.tre|.support|.nxs|.nex|treefile", tree))) {
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


#' @name write_newick
#' @description 
#' Writes a tree in newick format, while preserving bootstrap support values (if present).
#' 
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo', or a file.
#' @param output Name of the output file. If no output is specified, is provided it will use the tree object name with a '.nwk' suffix.
#' 
#' @return NULL
#' 
#' @examples
#' Read a tree in nexus format
#'   nx_tree <- read_tree(input_file = 'tree.nxs')
#'   
#' Write it in newick format
#'   write_newick( tree = nx_tree, output = 'nx_tree.nwk' )
#' 
#' @export


write_newick <- function(
    tree, 
    output = NULL
) {
  
  tree_obj <- .load_tree_object(tree)
  
  if (!is.null(tree_obj@data)) {
    
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
  } else {
     phylo_data <- tree_obj
  }
  
  if (!is.null(output)) {
    ape::write.tree(phy = phylo_data@phylo, output)
  } else {
      ape::write.tree(phy = phylo_data@phylo, paste0(tree, ".nwk"))
  }
}


#' @name export_plot
#' @description Function to export plot with provided options
#' 
#' @param plot A ggtree plot object.
#' @param output Name of the output file. 
#' @param format Format of the output file. Default: 'svg'
#' @param dpi Number of DPIs in output file. Default: 600
#' 
#' @return A ggplot object representing the phylogenetic tree and optionally an output svg file with this phylogenetic tree.
#' 
#' @examples
#'   Export a phylogenetic tree object
#'   export_plot (plot, output = "tree_plot_visualized.svg")
#' 
#' @export


export_plot <- function(
    plot,
    output = NULL,
    format = 'svg',
    dpi = 600
) {
  
  if (!(is.null(output))) {
    ggsave(
      plot = plot, 
      sprintf(paste0("%s.",format), output), 
      width = 10, 
      height = 10,
      dpi = dpi
    )
    
    message <- sprintf("Tree plotted and saved as %s.svg", output)
    message(message)
    
  } else {
      message("Plot will not be saved! Provide an output = <OUTPUT_NAME> for saving the output plot.")
  }
  
  return(plot)
}


#==================================== Manipulating Trees ====================================#


#' @name bootstrap_collapse
#' @description Collapse nodes in a phylogenetic tree where the bootstrap support is below a user-specified cutoff.
#'              If bootstrap support is in 0-1 scale, is converted at the 1-100 scale.
#'
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo'.
#' @param cutoff The threshold value for collapsing nodes based on bootstrap support. Default: 50
#'
#' @return A treedata object with collapsed branches where bootstrap support is less than the user-specified cutoff.
#'
#' @examples
#' Load a sample tree
#'   tree <- read_tree("path/to/tree.nwk")
#'
#' Collapse nodes with bootstrap support less than 0.7
#'   collapsed_tree <- bootstrap_collapse(tree, cutoff = 0.7)
#' 
#' 
#' @export

bootstrap_collapse <- function(tree, cutoff = 50) {
  tree_obj <- .load_tree_object(tree)
  
  if (all(tree_obj@data$support <= 1)) {
    message("Bootstrap values are in scale 0-1. Will be multiplied by 100 to be transferred to a scale of 1-100.")
    
    tree_obj@data$support <- sapply(
      tree_obj@data$support, function (bootstrap) {
        bootstrap <- bootstrap*100
        return(data.frame(bootstrap))
    })
  } 
  
  collapsed_tree <- as.treedata(
    as.polytomy(
      tree_obj, 
      feature = 'support', 
      fun = function(x) as.numeric(x) < cutoff
    )
  )
  
  return(collapsed_tree)
}



#' @name print_internal_nodes
#' @description This function generates a plot of the provided phylogenetic tree with node labels displayed. 
#'              Optional printing of tip labels matching a user-provided pattern.
#'
#' @param tree An object representing the phylogenetic tree. Should be of class 'treedata' or 'phylo'.
#' @param form Layout of the tree based on ggtree options. Layout can be rectangular, circular, roundrect, slanted, ellipse, fan, equal_angle, daylight. Default: "rectangular"
#' @param node_id_color Color of the printed node labels. Default: 'darkred'
#' @param pattern Pattern or list of patterns for matching tip labels of interest. If pattern is NULL, then only node labels will be printed. Default: NULL
#' @param tip_label_size Size of printed tip labels matching the user-specified pattern. Default: 2
#' 
#' @return A ggplot object representing the phylogenetic tree with node labels and optional highlighted tip labels.
#'         Used for identifying nodes for interest for extracting subtrees, flipping or labeling them in the visualize_tree function.
#' 
#' @examples
#' \dontrun{
#'   # Load a sample tree
#'   tree <- read_tree('tree.nwk')
#'
#'   # Print the tree with node labels and print tip labels matching "taxon_A"
#'   print_internal_nodes(tree, color = 'darkred', pattern = 'taxon_A' )
#' }
#' 
#' @export


print_internal_nodes <- function(
    tree, 
    form = "circular", 
    node_id_color = "darkred", 
    pattern = NULL,
    tip_label_size = 2, 
    output = "node_with.svg"
) {
  
  tree_obj <- .load_tree_object(tree)
  
  # Create the base tree plot with node labels
  plot <- ggtree(tree_obj, layout = form) +
    geom_nodelab(aes(label = node), hjust = -0.1, color = node_id_color, size = 3)
  
  if ( is.null(pattern) ) {
    print("Only node labels will be printed.")
    return(plot)
  } else {
    matching_labels <- unlist(
      sapply(pattern, function(patterns) {
        grep(patterns, tree_obj@phylo$tip.label, value = TRUE, fixed = TRUE)
      }))
    
    matching_dict <- data.frame(
      label = matching_labels,
      matching_flag = TRUE
    )
    
    plot <- plot %<+% matching_dict
    
    plot <- plot + 
      geom_tiplab(
        aes(
          label = ifelse(matching_flag == TRUE, label, NA), 
          fill = sample(colors(), 1)
        ), 
        size = tip_label_size, align.tip.label = TRUE
      )
  }
  
  export_plot(
    plot = plot, 
    output = output
  )
}


#' Group and Highlight Tree Descendants
#' 
#' @description
#' Creates a phylogenetic tree visualization with highlighted clades and their descendants
#' 
#' @param tree Phylogenetic tree object
#' @param highlight_clades Named vector of node IDs to highlight
#' @param colors Named vector of colors for each clade. If NULL, random colors are generated
#' @param form Tree layout format (default: "circular")
#' @param tiplabels Logical, whether to show tip labels (default: FALSE)
#' @param tip_label_size Numeric, size of tip labels. Ignored if tiplabels=FALSE
#' @param tip_fill_color Color for tip points when tiplabels=FALSE (default: "lightgrey")
#' @param legend_position Position of legend (default: "right")
#' @param legend_name Title for the legend (default: "")
#' @param output Output file path. If NULL, returns plot object
#' @param format Output file format (default: 'svg')
#' @param dpi Resolution for output file (default: 600)
#' 
#' @return Plot object or exported file if output is specified
#' 
#' @import ggtree
#' @import tidytree
#' 
#' @examples
#' highlight_nodes <- c("CYP2" = 1063, "CYP3" = 1426)
#' colors <- c("CYP2" = "gold", "CYP3" = "green")
#' group_descendants(tree = tree, highlight_clades = highlight_nodes, colors = colors)
#' 
#' @export

group_descendants <- function(
  tree,
  highlight_clades = NULL,
  colors = NULL,
  form = "circular",
  tiplabels = FALSE,
  tip_label_size = NULL,
  tip_fill_color = "lightgrey",
  legend_position = "right",
  legend_name = "",
  output = NULL,
  format = 'svg',
  dpi = 600
) {

  if (is.null(highlight_clades)) {
    stop("No grouping information was provided.")
  }
    
  if (is.null(colors)) {
    colors <- setNames(sample(rainbow(length(highlight_clades)), length(highlight_clades)), 
                      names(highlight_clades))
  }
  
  tree_obj <- .load_tree_object(tree)

  # For each highlighted node, find all its descendants and group them
  groups <- list()
  
  for(group_name in names(highlight_nodes)) {
    node_id <- highlight_nodes[[group_name]]
    descendants <- c(node_id, tidytree::offspring(tree_obj, node_id))
    groups[[group_name]] <- descendants
  }

  # Create basic phylogenetic tree with tip labels or not
  if (!tiplabels) {
    plot <- ggtree(tree_obj, layout = form) +
      geom_star(aes(subset = isTip, starshape = "circle"), 
                fill = tip_fill_color, size = 0.8) +
      scale_starshape_identity()

  } else {
    tip_label_size <- ifelse(is.null(tip_label_size), 1, tip_label_size)

    plot <- ggtree(tree_obj, layout = form) +
      geom_tiplab(
        size = tip_label_size, 
        color = "black",
        align.tip.label = TRUE
      )
  }

  # Create plot and apply groupings
  plot <- ggtree::groupOTU(plot, groups)
  
  highlight <- data.frame(
    Group = names(highlight_nodes),
    Label = highlight_nodes,
    Color = colors
  )

  # Add stylistic information
  plot <- plot + 
    aes(color=group) +
    scale_fill_manual(
      values = highlight$Color,
      labels = gsub("_", " ", highlight$Group),
      name = legend_name
    ) + 
    theme(legend.position = legend_position)

  export_plot(
      plot = plot, 
      output = output,
      format = format,
      dpi = dpi 
  )
}


#' @name extract_subtree
#' @description
#' Function to extract a subtree from a larger phylogenetic tree based on:
#'    1. the MRCA (Most Recent Common Ancestor) of two user-provided anchor tips
#'    2. a user-provided list of tip labels, which are contained within this subtree
#' while preserving branch lengths and bootstrap values.
#'
#' @param tree A phylogenetic tree object of class 'phylo' or 'treedata', a newick or a nexus file.
#' @param tip1,tip2 The labels or pattern ids of the two anchor tips.
#' @param taxon_list User-provided list of tip labels, the subtree which contains them will be extracted.
#'
#' @return A treedata object representing the extracted subtree with preserved branch lengths and bootstrap values.
#'
#' @examples
#' 
#' Extract a subtree by providing labels of anchor tips
#'   subtree <- extract_subtree(tree, tip1 = "tip1", tip2 = "tip2")
#'   
#' Extract a subtree by providing a list of taxons contained within this subtree
#'   subtree <- extract_subtree(tree, taxon_list = list_of_tip_labels)
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


# ==================================== Visualizing Trees ====================================#

#' @name highlight_tree
#' @description
#' Highlight Nodes and Descendant Leaves on a Phylogenetic Tree
#'
#' This function highlights specified nodes and their descendant leaves on a phylogenetic tree,
#' allowing for visualization emphasis on selected parts of the tree. Custom colors can be assigned
#' to different groups of nodes, or random colors can be automatically generated.
#'
#' @param tree A phylogenetic tree object, typically of class 'phylo' from the ape package,
#'        or a file path to a tree file in Newick format that can be read into such an object.
#' @param highlight_nodes An associative vector where names are group labels and values are
#'        node labels indicating the nodes and their descendants to be highlighted on the tree.
#' @param colors An associative vector where names are group labels and values are colors
#'        for the corresponding groups defined in `highlight_nodes`. If not provided,
#'        random colors will be assigned to each group.
#' @param form Layout of the tree, with options including "rectangular", "circular",
#'        "roundrect", "slanted", "ellipse", "fan", "equal_angle", and "daylight". Default: "rectangular"
#' @param tip_fill_color Fill color for the tip label symbols. Default: lightgrey
#' @param tip_border_color Border color for the tip label symbols. Default: black
#' @param legend_name A character string specifying the legend name for the highlighted plot. Default: ''
#' @param legend_position Position of the legend in the plot ("right", "left", "top", "bottom"). Default: 'right'
#' @param legend_text_position Horizontal alignment of the text within the legend ("left" or "right"). Default: 'right'
#' @param legend_orientation Orientation of the legend ("horizontal" or "vertical"). Default: 'vertical'
#' @param legend_key_size Size of the legend keys, specified in cm. Default: 0.1 because typically highlights are not that many in number
#' @param legend_font Font family for the text in the legend. Default: 'Arial'
#' @param legend_fontsize Size of the text in the legend, specified in points.
#' @param legend_font_face Fontface of legend text, same like for geom_text. Default: 'bold' 
#' @param legend_spacing_x Horizontal spacing between legend items, specified in cm. Default: 0.5 
#' @param legend_spacing_y Vertical spacing between legend items, specified in cm. Default: 0.5
#' @param legend_key_width Width of the keys in the legend, specified in cm. Default: 1
#' @param legend_title_hjust Horizontal justification of the legend title. Default: 0.5
#' @param output The filename where the plot will be saved.
#' @param format Format of the output file. Default: 'svg'
#' @param dpi Number of DPIs in output file. Default: 600
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


highlight_tree <- function(
    tree,
    highlight_nodes,
    colors = NULL,
    form = "circular",
    tiplabels = FALSE,
    tip_label_size = NULL,
    tip_fill_color = "lightgrey",
    tip_border_color = "white",
    legend_name = "",
    legend_position = "right",
    legend_text_position = "right",
    legend_orientation = "vertical",
    legend_key_size = 0.1,
    legend_font = "Arial",
    legend_fontsize = 11,
    legend_font_face = "bold",
    legend_spacing_x = 0.5,
    legend_spacing_y = 0.5,
    legend_key_width = 1,
    legend_title_hjust = 0.5,
    output = NULL,
    format = 'svg',
    dpi = 600
) {
  
  tree_obj <- .load_tree_object(tree)
  
  if (tiplabels == FALSE || tiplabels == F) {
    if ( !is.null(tip_label_size) ) {
      stop("This action is not permitted. Tip labels need to be enabled (tiplabels = TRUE)")
    }

    plot <- ggtree(tree_obj, layout = form) +
      geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) +
      scale_starshape_identity() 
  } else {
      if (is.null(tip_label_size)) {
        tip_label_size = 1
      }

      plot <- ggtree(tree_obj, layout = form) + 
        geom_tiplab(
          size = tip_label_size,
          color =  "black",
          align.tip.label = TRUE
        )
  }

  # Create a data frame for highlighting
  if (is.vector(highlight_nodes)) {
    if (is.null(colors)) {
      colors <- setNames(sample(rainbow(length(highlight_nodes)), length(highlight_nodes)), names(highlight_nodes))
    } else if (length(colors) != length(highlight_nodes)) {
        stop("Length of colors and highlight_nodes is not equal.")
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
      legend.position = legend_position,
      legend.direction = legend_orientation,
      legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
      legend.key.size = unit(legend_key_size, "cm"),
      legend.spacing.x = unit(legend_spacing_x, "cm"),
      legend.spacing.y = unit(legend_spacing_y, "cm"),
      legend.key.width = unit(legend_key_width, "cm"), 
      legend.title = element_text(hjust = legend_title_hjust)
    )
  
  export_plot(
    plot = plot, 
    output = output,
    format = format,
    dpi = dpi
  )
}


#' visualize_tree: Visualize Phylogenetic Tree with many customizable Features
#'
#' Main tree visualization function of the package, which provides comprehensive options 
#' for visualizing a phylogenetic tree, allowing for customizations such as 
#' displaying bootstrap values, tip labels as text and shape with different colors based on user-provided mappings, 
#' flipping nodes and label user-specified clades. It supports various tree layouts.
#'
#' @param tree A phylogenetic tree in Newick or Nexus format, or an object of treedata or phylo class.
#' @param form Layout of the tree, with options including "rectangular", "circular",
#'        "roundrect", "slanted", "ellipse", "fan", "equal_angle", and "daylight". Default: "rectangular"
#' @param flip_nodes Logical; whether to flip the tree at specified nodes.
#' @param node1,node2 The ID of the two nodes to flip based on their most common recent ancestor (MRCA)
#' @param tiplabels Logical; whether to print tip labels on the tree. Default: FALSE.
#' @param pattern_id A pattern used to match and selectively print tip labels as text.
#' @param tip_label_size Numeric specifying the size of the tip labels.
#' @param taxon_group_separator A separator string for splitting tip labels to generate color
#'        and shape mappings. If NULL, mappings are generated based on the whole tip label string. 
#' @param color A vector specifying color mappings for tip label text and/or shapes.
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
#' @param output The file name or path where the plot will be saved.
#' @param format Format of the output file. Default: 'svg'
#' @param dpi Number of DPIs in output file. Default: 600
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
    legend_position = "left",
    legend_text_position = "right",
    legend_orientation = "horizontal",
    legend_key_size = 0.1,
    legend_font = "Arial",
    legend_title_fontsize = 12,
    legend_fontsize = 11,
    legend_title_font_face = "bold",
    legend_font_face = "bold",
    legend_spacing_x = 0.5,
    legend_spacing_y = 0.5,
    legend_key_width = 1,
    legend_title_hjust = 0.5,
    legend_name = "",
    tip_shape_size = 3,
    tip_label_size = 2,
    bootstrap_numbers = FALSE,
    bootstrap_number_nudge_x = 0,
    bootstrap_number_nudge_y = 0.2,
    node_label_size = 3,
    nudge_x_label = 0.01,
    bootstrap_circles = FALSE,
    bootstrap_legend = FALSE,
    bootstrap_circle_size = 1.7,
    branch_length = TRUE,
    clades = NULL,
    labeldist = 1,
    bardist = 1,
    clade_label_size = 3,
    output = NULL,
    format = 'svg',
    dpi = 600
) {

  # Load tree with option to print branch length and/or flip nodes
  tree_obj <- .load_tree_object(tree)

  # If option for branch length is not TRUE, make branch lengths NULL
  if (!(branch_length == TRUE | branch_length == T)) {
    print(
      "Option branch length is deactivated. 
       As a result, the subtree will be extracted as 
       a cladogram with equal branch lengths."
    )
    tree_obj@phylo$edge.length <- NULL
  }

  if ( flip_nodes == TRUE ) {
    if ( !is.null(node1) || !is.null(node2) ) {
      plot <- ggtree::flip( ggtree(tree_obj, layout = form), node1, node2 )
    } else {
        stop(
          "Please provide node labels to flip. 
           Find node labels of interest using the 
           print_internal_nodes( tree, references = c(taxon_pattern1, taxon_pattern1) ) 
           function."
          )
    }
  } else {
      plot <- ggtree(
        tree_obj, layout = form
      )
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
          color =  "black",
          align.tip.label = TRUE,
          nudge_x = nudge_x_label
        )
    } else {
        matching_labels <- unlist(
          sapply(pattern_id, function(pattern) {
            matches <- grep(pattern, tree_obj@phylo$tip.label, fixed = TRUE, value = TRUE)
            return (matches)
          })
        )
        
        if (!is.null(names(color))) {
          color <- color
        } else {
            warning(
            "color has not been defined correctly.\nIt is an associative vector where names are taxon_names and values are colors.\n
             color <- c(
                  'taxonA' = 'red',
                  'taxonB' = 'green'
             )"
            )
        }
        
        if (length(matching_labels) > 0) {
          matching_dict <- data.frame(
            label = matching_labels,
            color = ifelse(
              !is.na(color[matching_labels]), 
              color[matching_labels], 
              "black" 
            )
          , matching_flag = TRUE
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
              align.tip.label = TRUE,
              nudge_x = nudge_x_label
            ) + scale_color_identity()
          
          if (length(pattern_id) < 10 ) {
            cat("Printing phylogenetic tree plot with tip labels matching", pattern_id, "!\n")
          } else {
              print("Printing phylogenetic tree plot with >10 tip label patterns!")
          }
        
      } else {
          message( cat("Provided", pattern_id, "was not found among tip labels! Printing phylogenetic tree plot without tip labels!") )
      }
    }
  } else { 
      if (!is.null(pattern_id)) {
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
        list( tip_label = label, group = sub(paste0(taxon_group_separator,".*"), "", label)) 
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
          s_color = color[
            match(sapply(
              taxa_dict, function(entry) 
                entry$group
              ), names(color)
            )
          ]
        )
      }
        
       if (!is.null(shape)) {
         
         # Associative vector of starshapes and numbers, got from ggstar's official GitHub page (https://github.com/xiangpin/ggstar/blob/master/R/geom-star.R, lines 109-141)
         starshape_table <- c(
           "pentagram"                = 1,
           "magen david"              = 2,
           "seven pointed star"       = 3,
           "anise star"               = 4,
           "regular pentagon"         = 5,
           "hexagon"                  = 6,
           "regular heptagon"         = 7,
           "regular octagon"          = 8,
           "anise star2"              = 9,
           "anise star3"              = 10,
           "regular triangle"         = 11,
           "rhombus"                  = 12,
           "square"                   = 13,
           "four-pointed star"        = 14,
           "circle"                   = 15,
           "heart"                    = 16,
           "left-triangle1"           = 17,
           "right-triangle1"          = 18,
           "left-triangle2"           = 19,
           "right-triangle2"          = 20,
           "rectangle"                = 21,
           "triangle star"            = 22,
           "regular triangle down"    = 23,
           "hexagonal star"           = 24,
           "ellipse"                  = 25,
           "thin triangle"            = 26,
           "anise star4"              = 27,
           "square diamond"           = 28,
           "plus filled"              = 29,
           "antiparallelogram"        = 30,
           "semicircle"               = 31
         )
         
         lapply(names(shape), function (name) {
           if (!(shape[name] %in% names(starshape_table))){
             print(show_starshapes()) # Show starshape to assist the user select the desired starshape
             stop(cat(shape[name], "is not found among ggstar's starshapes."))
           }
         })
         
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
      # Add ggtree aesthetics to plots based on the user-provided arguments
      if (!is.null(color) && !is.null(shape)) {
        
        # Function to extract group from tip label
        extract_group <- function(label) {
          if (!is.null(taxon_group_separator)) {
            parts <- strsplit(label, taxon_group_separator, fixed = TRUE)[[1]]
            return(parts[1])
          } else {
            return(label)
          }
        }

        # Merge color and shape mappings into a single dataframe to allow for mappings legend
        tip_mapping_df <- data.frame(
          label = tree_obj@phylo$tip.label,
          taxa_group = sapply(tree_obj@phylo$tip.label, extract_group),
          stringsAsFactors = FALSE
        )
        
        tip_mapping_df$s_color <- color[tip_mapping_df$taxa_group]
        tip_mapping_df$s_shape <- shape[tip_mapping_df$taxa_group]
        
        # Join data frames to create a unique mapping tibble
        plot <- plot %<+% tip_mapping_df 
        
        plot <- plot + geom_star(
          mapping = aes( 
            subset = isTip,
            fill = s_color,
            starshape = s_shape
          ), 
          size = tip_shape_size,
          color = "black"
        ) + 
          scale_fill_identity(
            breaks = unique(tip_mapping_df$s_color[!is.na(tip_mapping_df$s_color)]),
            labels = unique(tip_mapping_df$taxa_group[!is.na(tip_mapping_df$s_color)]),
            name = legend_name,
            guide = guide_legend(
              override.aes = list(
                starshape = unique(tip_mapping_df$s_shape[!is.na(tip_mapping_df$s_shape)])
              ),
              theme = theme(
                legend.position = legend_position,
                legend.direction = legend_orientation,
                legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
                legend.key.size = unit(legend_key_size, "cm"),
                legend.spacing.x = unit(legend_spacing_x, "cm"),
                legend.spacing.y = unit(legend_spacing_y, "cm"),
                legend.key.width = unit(legend_key_width, "cm"), 
                legend.title = element_text(face = legend_title_font_face, family = legend_font, size = legend_title_fontsize,  hjust = legend_title_hjust)
              )
            ) 
          ) + 
          scale_starshape_identity(guide = "none")
      
      } else if ( !is.null(color) && is.null(shape)) {
              plot <- plot %<+% tip_colors_df
             
              plot <- plot + geom_star(
                mapping = aes(
                  subset = isTip,
                  fill = ifelse(!is.na(taxa_colors), s_color, NA)
                ), color = "black",
                starshape = "circle",
                size = tip_shape_size,
              ) + 
                scale_fill_identity( 
                  breaks = tip_colors_df$s_color,
                  labels = gsub("_", " ", tip_colors_df$taxa_colors),
                  name = legend_name,
                  
                  guide = guide_legend(
                      theme = theme(
                        legend.position = 'none',
                        legend.direction = legend_orientation,
                        legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
                        legend.key.size = unit(legend_key_size, "cm"),
                        legend.spacing.x = unit(legend_spacing_x, "cm"),
                        legend.spacing.y = unit(legend_spacing_y, "cm"),
                        legend.key.width = unit(legend_key_width, "cm"), 
                        legend.title = element_text(face = legend_title_font_face, family = legend_font, size = legend_title_fontsize,  hjust = legend_title_hjust)
                    )
                  ),
                )
            
          } else if ( is.null(color) && !is.null(shape)) {
            
            # Translate shapes from words to numbers using the starshape table to get the number associated with each starshape
            tip_shapes_df$s_shape <- starshape_table[shape[match(sapply(taxa_dict, function(entry) entry$group), names(shape))]]
            
            plot <- plot %<+% tip_shapes_df
            
            plot <- plot + geom_star(
              mapping = aes( 
                subset = isTip & !is.na(s_shape),
                starshape = ifelse(!is.na(taxa_shapes), s_shape, NA)
              ), 
              fill = "black", 
              color = "black",
              size = tip_shape_size
            ) + 
              scale_starshape_identity(
                labels = gsub("_", " ", tip_shapes_df$taxa_shapes),
                breaks = tip_shapes_df$s_shape,
                name = legend_name,
                guide = guide_legend(
                    theme = theme(
                      legend.position = legend_position,
                      legend.direction = legend_orientation,
                      legend.text = element_text(face = legend_font_face, family = legend_font, size = legend_fontsize, color = "black"),
                      legend.key.size = unit(legend_key_size, "cm"),
                      legend.spacing.x = unit(legend_spacing_x, "cm"),
                      legend.spacing.y = unit(legend_spacing_y, "cm"),
                      legend.key.width = unit(legend_key_width, "cm"), 
                      legend.title = element_text(face = legend_title_font_face, family = legend_font, size = legend_title_fontsize,  hjust = legend_title_hjust)
                  )
                ),
              )
          }
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
      output = output,
      format = format,
      dpi = dpi 
    )
}