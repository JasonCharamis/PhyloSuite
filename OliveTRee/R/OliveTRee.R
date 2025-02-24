#' @title OliveTRee: Advanced Tree Manipulation and Visualization
#'
#' @description
#' A comprehensive R package that combines and extends functionality from ggtree, ape,
#' phytools and other phylogenetic tools for advanced tree manipulation and visualization.
#' Provides unified interfaces for common tree operations, enhanced visualization options,
#' and robust data handling capabilities.
#'
#' @docType package
#' @name OliveTRee
#' @importFrom ape read.tree write.tree read.nexus extract.clade
#' @importFrom phytools findMRCA read.newick
#' @importFrom treeio as.treedata
#' @importFrom dplyr left_join select mutate as_tibble
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual theme
#' @importFrom ggtree geom_tiplab geom_nodelab
#' @export NULL

#' @name .load_packages
#'
#' @description
#' Checks for required package dependencies and installs missing ones using BiocManager
#' when available, otherwise using install.packages(). Then loads the packages.
#'
#' @param tools Character vector of package names to check and install
#' @return NULL invisibly
#' @examples
#' \dontrun{
#' # Load multiple packages
#' .load_packages(c("ape", "phytools", "ggtree"))
#' }
#' @keywords internal
#' @export

.load_packages <- function(tools) {
  tmp <- as.data.frame(installed.packages())
  max_version <- max(as.numeric(substr(tmp$Built, 1, 1)))
  tmp <- tmp[as.numeric(substr(tmp$Built, 1, 1)) == max_version, ]

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

bioconductor_version <- ifelse(
  getRversion() >= "4.4",
  "3.19",
  ifelse(getRversion() >= "4.3", "3.18", "3.17")
)

#==================================== Reading, Loading and Exporting Trees ====================================#

#' @name read_tree
#'
#' @description
#' Reads a phylogenetic tree from various file formats with automatic bootstrap value detection.
#' Supports Newick, Nexus, and IQ-TREE formats.
#'
#' @param input_file Character string specifying path to tree file
#' @return A treedata object containing the phylogenetic tree and associated data
#' @examples
#' \dontrun{
#' # Read from Newick file
#' tree <- read_tree("path/to/tree.nwk")
#'
#' # Read from Nexus file
#' tree <- read_tree("path/to/tree.nexus")
#' }
#' @export

read_tree <- function(input_file) {
  if (grepl("\\.(newick|nwk|tre|tree)$", tolower(input_file))) {
    t <- phytools::read.newick(input_file)
    print(t)
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
  if (any(is.numeric(t$node.label) || is.integer(t$node.label))) {
    if (any(as.numeric(t$node.label) > 0)) {
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
  if (is.null(tree)) {
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
  # Enforce import as 'treedata' object
  tree_obj <- .load_tree_object(tree)

  if (!is.null(tree_obj@phylo)) {
    # Copy bootstrap values in the labels of the internal (non-Tip) nodes
    phylo_data <- treeio::as.treedata(
      cross_join(
        as_tibble(tree_obj@phylo),
        as_tibble(tree_obj@data)
      ) %>%
        mutate(
          isTip = ifelse(!is.na(label), TRUE, FALSE)
        ) %>%
        mutate(
          label = ifelse(
            isTip == FALSE,
            as.character(node.label),
            label
          )
        )
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
  plot_width = 20,
  plot_height = 15,
  dpi = 600,
  limitsize = FALSE
) {
  if (!(is.null(output))) {
    output_f = ifelse(
      grepl(format, output, fixed = TRUE),
      output,
      sprintf(paste0("%s.", format), output)
    )

    ggsave(
      plot = plot,
      output_f,
      width = plot_width,
      height = plot_width,
      dpi = dpi,
      limitsize = limitsize
    )

    message <- sprintf("Tree plotted and saved as %s", output_f)
    message("\033[34m\033[1m\033[30m", message)
  } else {
    message(
      "\033[34m\033[1m\033[30m",
      "Plot will not be saved! Provide an output = <OUTPUT_NAME> for saving the output plot."
    )
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
    message(
      "Bootstrap values are in scale 0-1. Will be multiplied by 100 to be transferred to a scale of 1-100."
    )

    tree_obj@data$support <- sapply(
      tree_obj@data$support,
      function(bootstrap) {
        bootstrap <- bootstrap * 100
        return(data.frame(bootstrap))
      }
    )
  } else if (all(tree_obj@data$support <= 10)) {
    message(
      "Bootstrap values are in scale 0-1. Will be multiplied by 100 to be transferred to a scale of 1-100."
    )

    tree_obj@data$support <- sapply(
      tree_obj@data$support,
      function(bootstrap) {
        bootstrap <- bootstrap * 10
        return(data.frame(bootstrap))
      }
    )
  }

  collapsed_tree <- as.treedata(
    as.polytomy(
      as.phylo(tree_obj),
      feature = 'node.label',
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
  output = "node_with.svg",
  format = "svg",
  dpi = 600
) {
  tree_obj <- .load_tree_object(tree)

  # Create the base tree plot with node labels
  plot <- ggtree(tree_obj, layout = form) +
    geom_nodelab(
      aes(label = node),
      hjust = -0.1,
      color = node_id_color,
      size = 3
    )

  if (is.null(pattern)) {
    print("Only node labels will be printed.")
    plot <- plot
  } else {
    matching_labels <- unlist(
      sapply(pattern, function(patterns) {
        grep(patterns, tree_obj@phylo$tip.label, value = TRUE, fixed = TRUE)
      })
    )

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
        size = tip_label_size,
        align.tip.label = TRUE
      )
  }

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
      stop(
        "Provided tip labels are NULL or were not found among tree tip labels."
      )
    } else if (any(sapply(list(tip1m, tip2m), function(x) length(x) > 1))) {
      stop("Provided tip names are not unique.")
    }

    subtree <- ape::extract.clade(
      t,
      node = phytools::findMRCA(t, c(tip1m, tip2m))
    )
  } else if (is.null(tip1) && is.null(tip2) && !is.null(taxon_list)) {
    taxon_tips <- t$tip.label %in% taxon_list

    if (any(taxon_tips)) {
      mrca <- getMRCA(t, which(taxon_tips))
      subtree <- ape::extract.clade(t, mrca)
    } else {
      stop("No matching species found in the tree.")
      return(NULL)
    }
  } else {
    stop(
      "Please choose either tip1,tip2 OR taxon_list as identifiers for extracting subtree."
    )
  }

  if (!any(c("support", "node.label") %in% names(tree_obj@phylo))) {
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


#' @name visualize_tree
#'
#' @description
#' Visualize Phylogenetic Tree with Advanced Customization
#' Main visualization function providing comprehensive options for displaying phylogenetic trees.
#' Supports multiple layouts, bootstrap value display methods, tip label customization, clade highlighting,
#' and various aesthetic adjustments. Integrates with ggplot2 and ggtree for high-quality visual output.
#'
#' @param tree Phylogenetic tree (treedata/phylo object or file path)
#' @param form Layout type: "rectangular", "circular", "roundrect", "slanted", "ellipse",
#'        "fan", "equal_angle", or "daylight" (default: "rectangular")
#' @param flip_nodes Logical; whether to flip nodes (default: FALSE)
#' @param node1,node2 Node IDs for flipping based on MRCA
#' @param tiplabels Logical; display tip labels (default: FALSE)
#' @param pattern_id Pattern to match for selective tip label display
#' @param color Named vector mapping taxa to colors
#' @param shape Named vector mapping taxa to ggstar shapes
#' @param taxon_group_separator Character separator in tip labels (default: "_")
#' @param taxon_group_field Position of taxon name in label (default: 1)
#' @param font Font family (default: "Arial")
#' @param legend_position Legend position (default: "left")
#' @param legend_text_position Text alignment in legend (default: "right")
#' @param legend_orientation "horizontal" or "vertical" (default)
#' @param legend_key_size Size of legend keys in cm (default: 0.1)
#' @param legend_title_fontsize Title font size (default: 12)
#' @param legend_fontsize Text font size (default: 11)
#' @param legend_title_font_face Title font style (default: "bold")
#' @param legend_font_face Text font style (default: "bold")
#' @param legend_spacing_x,legend_spacing_y Legend spacing in cm (default: 0.5)
#' @param legend_key_width Key width in cm (default: 1)
#' @param legend_title_hjust Title horizontal justification (default: 0.5)
#' @param legend_title Legend title (default: "")
#' @param tip_shape_size Size of tip shapes (default: 3)
#' @param tip_label_size Size of tip labels (default: 2)
#' @param bootstrap_numbers Logical; show bootstrap as numbers (default: FALSE)
#' @param bootstrap_number_nudge_x,bootstrap_number_nudge_y Offset for bootstrap numbers
#' @param node_label_size Size of node labels (default: 3)
#' @param nudge_x_label Horizontal label offset (default: 0.01)
#' @param bootstrap_circles Logical; show bootstrap as circles (default: FALSE)
#' @param bootstrap_legend Logical; show bootstrap legend (default: FALSE)
#' @param bootstrap_circle_size Size of bootstrap circles (default: 1.7)
#' @param branch_length Logical; show true branch lengths (default: TRUE)
#' @param clades Named vector mapping clade nodes to labels
#' @param labeldist,bardist Distance of clade labels/bars from nodes (default: 1)
#' @param clade_label_size Size of clade labels (default: 3)
#' @param tree_linewidth Width of tree branches (default: 0.5)
#' @param highlight_clades,highlight_nodes Nodes to highlight
#' @param highlight_colors Colors for highlighted regions
#' @param tip_fill_color Fill color for tip shapes (default: "lightgrey")
#' @param tip_border_color Border color for tip shapes (default: "white")
#' @param highlight_alpha Transparency of highlights (default: 0.3)
#' @param highlight_extend Extension of highlights (default: 0.10)
#' @param highlight_linetype Line type for highlights (default: 1)
#' @param highlight_linewidth Width of highlight lines (default: 0.9)
#' @param output Output file name
#' @param format Output file format (default: 'svg')
#' @param dpi Resolution in dots per inch (default: 600)
#'
#' @return A ggplot object representing the visualized phylogenetic tree
#'
#' @examples
#' \dontrun{
#' # Basic rectangular tree
#' visualize_tree(tree)
#'
#' # Circular tree with bootstrap circles
#' visualize_tree(tree, form = "circular", bootstrap_circles = TRUE)
#'
#' # Tree with colored tip shapes and custom legend
#' visualize_tree(tree,
#'   color = c("Species1" = "red", "Species2" = "blue"),
#'   shape = c("Species1" = "circle", "Species2" = "star"),
#'   legend_title = "Species"
#' )
#'
#' # Tree with highlighted clades
#' visualize_tree(tree,
#'   highlight_clades = c("CladeA" = 15, "CladeB" = 20),
#'   highlight_colors = c("CladeA" = "pink", "CladeB" = "lightblue")
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual theme
#' @importFrom ggtree ggtree geom_tiplab geom_nodelab geom_cladelab
#' @export
#'

visualize_tree <- function(
  tree = NULL,
  form = "rectangular",
  flip_nodes = FALSE,
  node1,
  node2 = NULL,
  tiplabels = FALSE,
  pattern_id = NULL,
  color = NULL,
  shape = NULL,
  taxon_group_separator = "_",
  taxon_group_field = 1,
  font = "Arial",
  legend_position = "left",
  legend_text_position = "right",
  legend_orientation = "horizontal",
  legend_key_size = 0.1,
  legend_title_fontsize = 12,
  legend_fontsize = 11,
  legend_title_font_face = "bold",
  legend_font_face = "bold",
  legend_spacing_x = 0.5,
  legend_spacing_y = 0.5,
  legend_key_width = 1,
  legend_title_hjust = 0.5,
  legend_title = "",
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
  tree_linewidth = 0.5,
  highlight_clades = NULL,
  highlight_nodes = NULL,
  highlight_colors = NULL,
  highlight_legend_title = "Group",
  highlight_border_color = "black",
  tip_fill_color = "lightgrey",
  tip_border_color = "white",
  highlight_alpha = 0.3,
  highlight_extend = 0.10,
  highlight_linetype = 1,
  highlight_linewidth = 0.9,
  output = NULL,
  format = 'svg',
  dpi = 600,
  plot_width = 20,
  plot_height = 15
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

  # Initialize base plot
  if (flip_nodes == TRUE) {
    if (!is.null(node1) || !is.null(node2)) {
      plot <- ggtree::flip(
        ggtree(tree_obj, width = tree_linewidth, layout = form),
        node1,
        node2
      )
    } else {
      stop(
        "Please provide node labels to flip.
         Find node labels of interest using the
         print_internal_nodes(tree, references = c(taxon_pattern1, taxon_pattern1))
         function."
      )
    }
  } else {
    plot <- ggtree(tree_obj, width = tree_linewidth, layout = form)
  }

  # Enforce graphic parameters
  plot <- plot +
    geom_tree(size = 0.5)

  # Handle highlight_nodes if provided (highlight_tree functionality)
  if (!is.null(highlight_nodes)) {
    if (!is.vector(highlight_nodes) || is.null(names(highlight_nodes))) {
      stop("highlight_nodes should be a named vector.")
    } else if (!is.null(highlight_clades)) {
      stop(
        "You should choose either highlight_nodes or highlight_clades. Both cannot be enabled at the same time."
      )
    }

    # Sort named vectors based on alphabetical order of names
    highlight_nodes <- highlight_nodes[order(names(highlight_nodes))]

    # Ensure highlight_colors is correctly sorted
    if (!is.null(highlight_colors) && is.vector(highlight_colors)) {
      highlight_colors <- highlight_colors[order(names(highlight_colors))]
    }

    # Handle missing or inconsistent colors for highlighting
    if (
      is.null(highlight_colors) ||
        length(highlight_colors) != length(highlight_nodes)
    ) {
      highlight_colors <- setNames(
        sample(rainbow(length(highlight_nodes))),
        names(highlight_nodes)
      )
    }

    # Create highlight data frame for nodes
    highlight <- data.frame(
      Group = names(highlight_nodes),
      Label = unname(highlight_nodes),
      Color = sapply(names(highlight_nodes), function(group) {
        if (
          !is.null(highlight_colors) &&
            group %in% names(highlight_colors) &&
            !is.na(highlight_colors[[group]])
        ) {
          highlight_colors[[group]]
        } else {
          "black" # Default color if no mapping exists
        }
      }),
      stringsAsFactors = FALSE
    )

    # Add highlighting for nodes
    plot <- plot +
      geom_hilight(
        data = highlight,
        aes(node = Label, color = Group),
        color = highlight_border_color,
        alpha = highlight_alpha,
        extend = highlight_extend,
        linetype = highlight_linetype,
        linewidth = highlight_linewidth
      ) +
      scale_color_manual(
        name = highlight_legend_title,
        values = setNames(highlight$Color, highlight$Group), # Ensure node color mapping
        guide = "legend"
      )
  } else if (!is.null(highlight_clades)) {
    # Handle highlight_clades if provided (group_descendants functionality)

    if (is.null(highlight_colors)) {
      warning(
        "No colors provided for highlighting. Using black as default color."
      )
    }

    # Sort highlight_clades by group names
    highlight_clades <- highlight_clades[order(names(highlight_clades))]

    # Create highlight data frame for clades with default black for any missing colors
    highlight <- data.frame(
      Group = names(highlight_clades),
      Label = highlight_clades,
      Color = sapply(names(highlight_clades), function(group) {
        if (
          is.null(highlight_colors) ||
            !group %in% names(highlight_colors) ||
            is.na(highlight_colors[[group]])
        ) {
          "black" # Default color if no mapping exists
        } else {
          highlight_colors[[group]]
        }
      }),
      stringsAsFactors = FALSE
    )

    # Create groups for clades highlighting
    groups <- lapply(names(highlight_clades), function(group_name) {
      node_id <- highlight_clades[[group_name]]
      descendants <- c(node_id, unlist(tidytree::offspring(tree_obj, node_id)))
      return(descendants)
    })

    names(groups) <- names(highlight_clades)

    # Apply grouping for clades
    plot <- ggtree::groupOTU(plot, groups, group_name = highlight_legend_title)

    # Apply group coloring for clades
    plot <- plot +
      aes(color = !!sym(highlight_legend_title)) + # Ensure correct column mapping
      scale_color_manual(
        name = highlight_legend_title,
        values = setNames(highlight$Color, highlight$Group), # Ensure proper color mapping for clades
        guide = "legend"
      )
  }

  # Manipulate bootstrap values and control if and how to print them
  if (
    any(!is.null(tree_obj@data$support)) ||
      any(!is.na(tree_obj@data$support)) ||
      any(!is.null(tree_obj@phylo$node.label)) ||
      any(!is.na(tree_obj@phylo$node.label))
  ) {
    if (
      all(
        tree_obj@data$support <= 10,
        tree_obj@data$node.label <= 10,
        na.rm = TRUE
      ) ==
        TRUE
    ) {
      warning(
        "Provided bootstrap values were in scale of 1 and were multiplied by 100 to change scales."
      )
      if (!is.na(tree_obj@data$support)) {
        tree_obj@data$support <- sapply(
          tree_obj@data$support,
          function(x) x * 100
        )
      } else if (!is.na(tree_obj@phylo$node.label)) {
        tree_obj@data$support <- sapply(
          tree_obj@phylo$node.label,
          function(x) x * 100
        )
      }
    } else {
      tree_obj@data$support <- tree_obj@data$node.label
    }

    if (bootstrap_numbers == TRUE) {
      plot <- plot +
        geom_nodelab(
          aes(label = round(support, 2)),
          nudge_x = bootstrap_number_nudge_x,
          nudge_y = bootstrap_number_nudge_y,
          size = node_label_size,
          family = font
        )
      message(
        "\033[34m\033[1m\033[30m",
        "Bootstrap values displayed as labels."
      )
    } else if (bootstrap_circles == TRUE) {
      tree_obj <- treeio::as.treedata(
        as_tibble(tree_obj) %>%
          mutate(
            bs_fill = case_when(
              support < 50 | is.na(support) ~ "snow2",
              support >= 50 & support < 75 ~ "grey",
              support >= 75 ~ "black"
            ),
            bs_group = factor(
              case_when(
                # Add grouping for legend
                support < 50 | is.na(support) ~ "< 50%",
                support >= 50 & support < 75 ~ "50-75%",
                support >= 75 ~ "≥ 75%"
              ),
              levels = c("< 50%", "50-75%", "≥ 75%") # Specify order
            )
          )
      )

      plot <- plot +
        geom_point2(
          aes(
            subset = !isTip,
            fill = tree_obj@data$bs_group
          ),
          colour = "black",
          shape = 21,
          alpha = 2,
          stroke = 0.30,
          size = bootstrap_circle_size
        ) +
        scale_fill_manual(
          values = c(
            "< 50%" = "snow2",
            "50-75%" = "grey",
            "≥ 75%" = "black"
          ),
          name = "Bootstrap Support",
          na.value = NA
        )

      print("Bootstrap values displayed as circles.")
    } else if (bootstrap_numbers == TRUE && bootstrap_circles == TRUE) {
      stop("Please select either bootstrap_numbers or bootstrap_circles.")
    } else {
      message("Bootstrap values will not be printed.")
    }
  } else {
    message("Bootstrap support was not found in tree.")
  }

  # Options for printing tip labels as text
  if (tiplabels == TRUE) {
    if (is.null(pattern_id)) {
      plot <- plot +
        geom_tiplab(
          size = tip_label_size,
          color = "black",
          align.tip.label = TRUE,
          nudge_x = nudge_x_label,
          family = font
        )
    } else {
      matching_labels <- unlist(
        sapply(pattern_id, function(pattern) {
          matches <- grep(
            pattern,
            tree_obj@phylo$tip.label,
            fixed = TRUE,
            value = TRUE
          )
          return(matches)
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
          ),
          matching_flag = TRUE
        )

        plot <- plot %<+% matching_dict

        plot <- plot +
          geom_tiplab(
            aes(
              subset = isTip,
              label = label,
              color = ifelse(matching_flag == TRUE, color, NA),
              family = font
            ),
            size = tip_label_size,
            align.tip.label = TRUE,
            nudge_x = nudge_x_label
          ) +
          scale_color_identity()

        if (length(pattern_id) < 10) {
          cat(
            "Printing phylogenetic tree plot with tip labels matching",
            pattern_id,
            "!\n"
          )
        } else {
          print("Printing phylogenetic tree plot with >10 tip label patterns!")
        }
      } else {
        message(cat(
          "Provided",
          pattern_id,
          "was not found among tip labels! Printing phylogenetic tree plot without tip labels!"
        ))
      }
    }
  } else {
    if (!is.null(pattern_id)) {
      stop("Please enable tip label printing with tiplabels = TRUE.")
    } else {
      print("Printing phylogenetic tree plot without tip labels!")
    }
  }

  # Visualize phylogenetic tree with color and shape mappings
  if (!is.null(color) || !is.null(shape)) {
    # Use the taxon group separator to generate mappings
    if (!is.null(taxon_group_separator)) {
      # First field in tip label corresponds to the taxon name (Default: TAXON_NAME|taxon_group_separator|OTHER_INFO)
      if (taxon_group_field == 1) {
        taxa_dict <- lapply(tree_obj@phylo$tip.label, function(label) {
          list(
            tip_label = label,
            group = sub(paste0(taxon_group_separator, ".*"), "", label)
          )
        })
      } else if (taxon_group_field == 2) {
        # Second field in tip label corresponds to the taxon name (OTHER_INFO|taxon_group_separator|TAXON_NAME)
        taxa_dict <- lapply(tree_obj@phylo$tip.label, function(label) {
          list(
            tip_label = label,
            group = sub(paste0(taxon_group_separator, ".*"), "", label)
          )
        })
      } else {
        taxa_dict <- lapply(tree_obj@phylo$tip.label, function(label) {
          list(tip_label = label, group = label)
        })
      }
    }

    if (length(taxa_dict) > 0) {
      tip_labels <- sapply(taxa_dict, function(entry) entry$tip_label)
      taxa_names <- sapply(taxa_dict, function(entry) entry$group)

      if (!is.null(color)) {
        missing_species_color <- setdiff(
          unique(sapply(taxa_dict, function(entry) entry$group)),
          names(color)
        )

        if (length(missing_species_color) > 0) {
          warning(
            "\033[34m\033[1m\033[30m",
            paste(
              "Color mapping not found for the following species: ",
              paste(missing_species_color, collapse = ", ")
            )
          )
        }

        tip_colors_df <- data.frame(
          label = tip_labels,
          taxa_colors = taxa_names,
          s_color = color[
            match(
              sapply(
                taxa_dict,
                function(entry) entry$group
              ),
              names(color)
            )
          ]
        )
      }

      if (!is.null(shape)) {
        # Associative vector of starshapes and numbers
        starshape_table <- c(
          "pentagram" = 1,
          "magen david" = 2,
          "seven pointed star" = 3,
          "anise star" = 4,
          "regular pentagon" = 5,
          "hexagon" = 6,
          "regular heptagon" = 7,
          "regular octagon" = 8,
          "anise star2" = 9,
          "anise star3" = 10,
          "regular triangle" = 11,
          "rhombus" = 12,
          "square" = 13,
          "four-pointed star" = 14,
          "circle" = 15,
          "heart" = 16,
          "left-triangle1" = 17,
          "right-triangle1" = 18,
          "left-triangle2" = 19,
          "right-triangle2" = 20,
          "rectangle" = 21,
          "triangle star" = 22,
          "regular triangle down" = 23,
          "hexagonal star" = 24,
          "ellipse" = 25,
          "thin triangle" = 26,
          "anise star4" = 27,
          "square diamond" = 28,
          "plus filled" = 29,
          "antiparallelogram" = 30,
          "semicircle" = 31
        )

        lapply(names(shape), function(name) {
          if (!(shape[name] %in% names(starshape_table))) {
            print(show_starshapes())
            stop(cat(shape[name], "is not found among ggstar's starshapes."))
          }
        })

        if (tiplabels == TRUE) {
          warning(
            "If both tip labels and shapes are enabled, shape mapping will override the text tip labels provided."
          )
        }

        missing_taxa_shape <- setdiff(
          unique(sapply(taxa_dict, function(entry) entry$group)),
          names(shape)
        )

        if (length(missing_taxa_shape) > 0) {
          warning(paste(
            "Shape mapping not found for the following species: ",
            paste(missing_taxa_shape, collapse = ", ")
          ))
        }

        tip_shapes_df <- data.frame(
          label = tip_labels,
          taxa_shapes = taxa_names,
          s_shape = shape[match(
            sapply(taxa_dict, function(entry) entry$group),
            names(shape)
          )]
        )
      }

      # Add ggtree aesthetics based on provided arguments
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

        # Merge color and shape mappings
        tip_mapping_df <- data.frame(
          label = tree_obj@phylo$tip.label,
          taxa_group = sapply(tree_obj@phylo$tip.label, extract_group),
          stringsAsFactors = FALSE
        )

        tip_mapping_df$s_color <- color[tip_mapping_df$taxa_group]
        tip_mapping_df$s_shape <- shape[tip_mapping_df$taxa_group]

        plot <- plot %<+% tip_mapping_df

        plot <- plot +
          geom_star(
            mapping = aes(
              subset = isTip,
              fill = s_color,
              starshape = s_shape
            ),
            size = tip_shape_size,
            color = "black"
          ) +
          scale_fill_identity(
            name = legend_title,
            breaks = unique(tip_mapping_df$s_color[
              !is.na(tip_mapping_df$s_color)
            ]),
            labels = unique(tip_mapping_df$taxa_group[
              !is.na(tip_mapping_df$s_color)
            ]),
            guide = guide_legend(
              override.aes = list(
                # In starshapes associated with tips and colors
                starshape = tip_mapping_df$s_shape[match(
                  unique(tip_mapping_df$s_color),
                  tip_mapping_df$s_color
                )]
              ),
              theme = theme(
                legend.position = legend_position,
                legend.direction = legend_orientation,
                legend.text = element_text(
                  face = legend_font_face,
                  family = font,
                  size = legend_fontsize,
                  color = "black"
                ),
                legend.key.size = unit(legend_key_size, "pt"),
                legend.key.spacing.x = unit(legend_spacing_x, "pt"),
                legend.key.spacing.y = unit(legend_spacing_y, "pt"),
                legend.key.width = unit(legend_key_width, "pt"),
                legend.title = element_text(
                  face = legend_title_font_face,
                  family = font,
                  size = legend_title_fontsize,
                  hjust = legend_title_hjust
                )
              )
            )
          ) +
          scale_starshape_identity(guide = "none")
      } else if (!is.null(color) && is.null(shape)) {
        plot <- plot %<+% tip_colors_df

        plot <- plot +
          geom_star(
            mapping = aes(
              subset = isTip,
              fill = ifelse(!is.na(taxa_colors), s_color, NA)
            ),
            color = "black",
            starshape = "circle",
            size = tip_shape_size
          ) +
          scale_fill_identity(
            breaks = tip_colors_df$s_color,
            labels = gsub("_", " ", tip_colors_df$taxa_colors),
            name = legend_title,
            guide = guide_legend(
              theme = theme(
                legend.position = 'none',
                legend.direction = legend_orientation,
                legend.text = element_text(
                  face = legend_font_face,
                  family = font,
                  size = legend_fontsize,
                  color = "black"
                ),
                legend.key.size = unit(legend_key_size, "cm"),
                legend.key.spacing.x = unit(legend_spacing_x, "cm"),
                legend.key.spacing.y = unit(legend_spacing_y, "cm"),
                legend.key.width = unit(legend_key_width, "cm"),
                legend.title = element_text(
                  face = legend_title_font_face,
                  family = font,
                  size = legend_title_fontsize,
                  hjust = legend_title_hjust
                )
              )
            )
          )
      } else if (is.null(color) && !is.null(shape)) {
        # Translate shapes from words to numbers
        tip_shapes_df$s_shape <- starshape_table[shape[match(
          sapply(taxa_dict, function(entry) entry$group),
          names(shape)
        )]]

        plot <- plot %<+% tip_shapes_df

        plot <- plot +
          geom_star(
            mapping = aes(
              subset = isTip & !is.na(s_shape),
              starshape = ifelse(!is.na(taxa_shapes), s_shape, NA)
            ),
            fill = "black",
            color = "black",
            size = tip_shape_size
          ) +
          scale_starshape_identity(
            labels = unique(names(tip_shapes_df$s_shape)),
            breaks = unique(tip_shapes_df$s_shape),
            name = legend_title,
            guide = guide_legend(
              theme = theme(
                legend.position = legend_position,
                legend.direction = legend_orientation,
                legend.text = element_text(
                  face = legend_font_face,
                  family = font,
                  size = legend_fontsize,
                  color = "black"
                ),
                legend.key.size = unit(legend_key_size, "cm"),
                legend.key.spacing.x = unit(legend_spacing_x, "cm"),
                legend.key.spacing.y = unit(legend_spacing_y, "cm"),
                legend.key.width = unit(legend_key_width, "cm"),
                legend.title = element_text(
                  face = legend_title_font_face,
                  family = font,
                  size = legend_title_fontsize,
                  hjust = legend_title_hjust
                )
              )
            )
          )
      }
    }
  } else {
    # When using highlighting features, handle basic tip visualization
    if (!tiplabels) {
      plot <- plot +
        geom_star(
          aes(subset = isTip, starshape = "circle"),
          fill = tip_fill_color,
          color = tip_border_color,
          size = tip_shape_size
        ) +
        scale_starshape_identity()
    } else {
      if (is.null(tip_label_size)) {
        tip_label_size <- 1
      }
      plot <- plot +
        geom_tiplab(
          size = tip_label_size,
          color = "black",
          align.tip.label = TRUE,
          nudge_x = nudge_x_label,
          family = font
        )
    }
  }

  # Add clade labels
  if (!is.null(clades)) {
    for (i in names(clades)) {
      plot <- plot +
        geom_cladelab(
          node = clades[i],
          label = i,
          align = TRUE,
          fill = 'black',
          offset.text = labeldist,
          barsize = 0.9,
          offset = bardist,
          clade_label_size = clade_label_size,
          family = font
        )
    }
  }

  # Modify legend style
  plot <- plot +
    theme(
      legend.position = legend_position,
      legend.direction = legend_orientation,
      legend.text = element_text(
        face = legend_font_face,
        family = font,
        size = legend_fontsize,
        color = "black"
      ),
      legend.key.size = unit(legend_key_size, "cm"),
      legend.spacing.x = unit(legend_spacing_x, "cm"),
      legend.spacing.y = unit(legend_spacing_y, "cm"),
      legend.key.width = unit(legend_key_width, "cm"),
      legend.title = element_text(
        face = legend_title_font_face,
        family = font,
        size = legend_title_fontsize,
        hjust = legend_title_hjust
      )
    )

  # Export plot with provided options
  export_plot(
    plot = plot,
    output = output,
    format = format,
    dpi = dpi,
    plot_width = plot_width,
    plot_height = plot_height
  )
}
