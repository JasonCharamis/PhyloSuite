#Library of functions for advanced tree manipulation and visualization using ggtree, ape, phytools, and other related tools.

# Function to check if a package is installed, and if not, install it.
# Function to install or load a package
load_packages <- function( tools ) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built, 1, 1)))
  tmp <- tmp[as.numeric(substr(tmp$Built, 1, 1)) == max_version, ]

  if ( is.vector(tools) == TRUE ) {
    for ( pkg in tools ) {
      if ( pkg %in% tmp$Package ) {
        library(tools, character.only = TRUE)
      } else {
          print(sprintf("%s %s", tools, "is not installed. Installing it!"))
          is_available <- BiocManager::available(tools)
          
          if (any(is_available == "TRUE")) {
            BiocManager::install(tools, dependencies = TRUE, update = TRUE)
          } else {
            install.packages(tools, dependencies = TRUE, ask = FALSE, reinstall = TRUE)
        }
      }
    }
  }
}

# Load required packages or install them if necessary
dependencies <- c("ape", "phytools", "treeio", "TreeTools", "ggstar", 
                  "ggtree", "ggplot2", "dplyr", "stringi", "stringr",
                  "tidytree")

for (d in dependencies ) {
  load_packages(d)
}


#==================================== TREE MANIPULATION FUNCTIONS ====================================#

# Function to read and preprocess a tree
read_tree <- function(nwk) {
  tree <- ape::read.tree(nwk)
  return(Preorder(tree))
}

# Function to print a tree with node IDs
node_ids <- function(tree, reference = NULL) {
  # Create the base tree plot with node labels
  tree_plot <- ggtree(tree, layout = "rectangular") + 
    geom_nodelab(aes(label = node), hjust = -0.1, color = "red")
  
  if (is.null(reference)) {
    return(tree_plot) 
  } else {
    # Filter the tips based on the reference string
    tips_to_label <- tree$tip.label[grep(reference, tree$tip.label)]
    
    # Create the plot with tip labels for the specified tips
    labeled_tree <- tree_plot + geom_tiplab2(
      aes(subset = label %in% tips_to_label, label = label, color = "red"), 
      geom = "label"
    )
    return(labeled_tree)
  }
}

# Function to collapse nodes based on bootstrap support, returns a "phylo" object
bootstrap_collapse <- function(tree, cutoff) { 
  return(as.polytomy(tree, feature = 'node.label', fun = function(x) as.numeric(x) < cutoff))
}

# Function to flip based on descendant nodes (internal nodes or leaves if node is terminal)
flip_node <- function(tree, node1, node2) {
  return(as.phylo(flip(ggtree(tree), node1, node2)))
}

# Function to group all descendant branches of node(s)
group_descendants <- function(tree, node1, node2 = "", node3 = "", node4 = "") {
  return(groupClade(tree, .node = c(node1, node2, node3, node4)))
}

# Function to extract a subtree by finding the MRCA of two anchor nodes while preserving branch lengths and bootstrap values 
extract_subtree <- function(t, tip1, tip2, bl = TRUE) {
  if (!(bl == "TRUE" || bl == "T")) {
    t$edge.length <- NULL
  } 
  return(ape::extract.clade(t, node = phytools::findMRCA(t, c(tip1, tip2))))
}

#==================================== TREE VISUALIZATION FUNCTIONS ====================================#

# Function to include bootstrap values as colored circles
bootstrap_circles <- function(x) {
  
  if (!is.null(node.label(x))) {
    bootstrap_colors <- c('(75,100]' = "black", '(50,75]' = "grey", '(0,50]' = "snow2")
    categories <- c('(75,100]', '(50,75]', '(0,50]')

    bs_tibble <- tibble(
      node = 1:Nnode(x@phylo) + Ntip(x@phylo),
      bootstrap = as.numeric(x@phylo$node.label),
      cut(bootstrap, c(0, 50, 75, 100))
    ) 
    return(data.frame(
      node = bs_tibble$node,
      bootstrap = as.numeric(bs_tibble$bootstrap),
      category = cut(bs_tibble$bootstrap, c(0, 50, 75, 100)),
      b_color = bootstrap_colors[match(cut(bs_tibble$bootstrap, c(0, 50, 75, 100)), categories)]
    ))
  } else {
      print ("Bootstrap values do not exist.")
      return (NULL)
    }
}

# Function to highlight nodes on a tree
highlight_tree <- function(tree, highlight_nodes, colors = NULL, layout = "circular", name = NULL, ...) {
  # Check if highlight_nodes is a list or a vector
  if (is.list(highlight_nodes)) {
    # It's a list, create a data frame with random colors
    if (is.null(colors)) {
      # Generate random colors if colors are not provided
      highlight <- data.frame(
        Groups = names(highlight_nodes),  
        Label = highlight_nodes,
        Color = sample(colors(), length(highlight_nodes))
      )
    } else {
      # Use provided colors
      highlight <- data.frame(
        Groups = names(highlight_nodes),  
        Label = highlight_nodes,
        Color = colors
      )
    }
  } else if (is.vector(highlight_nodes)) {
    # It's a vector, create a data frame with provided colors or generate random colors
    if (is.null(colors)) {
      # Generate random colors if colors are not provided
      highlight <- data.frame(
        Groups = NULL,  # You can specify the groups if needed
        Label = highlight_nodes,
        Color = sample(colors(), length(highlight_nodes))
      )
    } else {
      # Use provided colors
      highlight <- data.frame(
        Groups = NULL,  # You can specify the groups if needed
        Label = highlight_nodes,
        Color = colors
      )
    }
  } else {
    stop("highlight_tree must be either a list or a vector.")
  }
  
  # Create the ggtree plot
  if (any(grepl(".nwk|.tre", tree))) {
    x <- read_tree(tree)
    tree <- ggtree(x, layout = layout) + 
      geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) + 
      scale_starshape_identity() + scale_fill_identity()
  } else {
    tree <- ggtree(tree, layout = layout) + 
      geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) + 
      scale_starshape_identity() + scale_fill_identity()
  }
  
  # Highlight the nodes
  if (!is.null(highlight$Groups)) {
    treef <- tree + geom_hilight(
      data = highlight,
      aes(node = Label, fill = Groups),
      alpha = 0.3, extend = 0.10, linetype = 1, linewidth = 0.9, 
      colour = "black", show.legend = TRUE
    ) + scale_fill_manual(values = highlight$Color)
  } else {
    treef <- tree + geom_hilight(
      data = highlight,
      aes(node = Label, fill = highlight$Color),
      alpha = 0.3, extend = 0.10, linetype = 1, linewidth = 0.9, 
      colour = "black", show.legend = TRUE
    ) 
  }
  
  # Save the highlighted plot
  if (exists("treef")) {
    if (any(grepl(".nwk|.tre", tree))) {
      ggsave(plot = treef, sprintf("%s_highlighted.svg", sub(".nwk|.tre", "", tree)), dpi = 600)
      sprintf("Tree plotted and saved as %s_highlighted.svg", sub(".nwk|.tre", "", tree))
    } else {
      ggsave(plot = treef, "tree_plot_highlighted.svg", dpi = 600)
      print("Tree plotted and saved as tree_plot_highlighted.svg!")
    }
    return(treef)
  }
}

# Node coloring is deprecated from the visualize_tree function and will be implemented only in the highlight_tree argument.
# Now, I have to modify the code, to accept visualize tree objects in the highlight_tree function
# Function to visualize a tree with a wide variety of options and customizable features
# Reference taxons should be defined as reference = 'taxon_ID'.

visualize_tree <- function(tree, color = NULL, shape = NULL, bootstrap_circles = TRUE,
                           clades = NULL, labels = NULL, fontsize = 7, labeldist = 0.1, 
                           ...) {

  tree_df <- as_tibble(tree)

  # Initialize graphics object
  plot <- ggtree(tree) 
  #+ geom_tiplab()
  
  # Manipulate tip labels to generate species and discard tip labels which match reference
  references <- list(...)
  reference_args <- names(references)
  
  if (length(reference_args) == 0) {
    species_dict <- lapply(seq_along(tree$tip.label), function(i) {
      list(tip_label = tree$tip.label[i], species = sub("_.*", "", tree$tip.label[i]))
    })
  } else {  
    reference <- unlist(references)
    ref_species <- tree$tip.label[tree$tip.label %in% reference]
    non_ref_species <- tree$tip.label[!(tree$tip.label %in% reference)]
    
    species_dict <- lapply(non_ref_species, function(label) {
      list(tip_label = label, species = sub("_.*", "", label))
    })
  }
  
  # Include text in reference species tip labels 
  if ( !(length ( reference_args ) == 0) ) {
    plot <- plot + geom_text(data = data.frame(tip_label = ref_species), aes(label = tip_label), 
                             position = position_nudge(x = 0.08), vjust = 0.5, hjust = 0.9, size = 3.5, 
                             family = "Arial", fontface = "bold")
  }

  if ( is.null(color) && is.null(shape) ) {
    print ('No color or shape mappings provided. Will just print the phylogenetic tree.')
  }
  
  # Add color and/or shape mappings
  else {
    species_names <- sapply(species_dict, function(entry) entry$species)
    x <- data.frame()
  
    if ( !is.null(color) && !is.null(shape)) {
        tip_colors_df <- data.frame(
        label = sapply(species_dict, function(entry) entry$tip_label),
        species_colors = species_names,
        color = color[match(sapply(species_dict, function(entry) entry$species), unique(species_names))]
      )
      
      tip_shapes_df <- data.frame(
        label = sapply(species_dict, function(entry) entry$tip_label),
        species_shapes = species_names,
        shape = shape[match(sapply(species_dict, function(entry) entry$species), unique(species_names))]
      )
      
      # Join tree, color and shapes dataframe columns to generate a unique mapping tibble df,
      # which will be used as reference for producing the color and shape mappings in ggtree.
      # Join tip colors and shapes in tree object based on tip label
      
      x <- full_join(tree_df, tip_colors_df, by = 'label')
      x <- full_join(x, tip_shapes_df, by = 'label')
      plot <- plot + geom_star(aes(x=x, subset = isTip, fill = species_colors, starshape = shape), size = 3, show.legend = TRUE)
    }
    
    if ( !is.null(color) && is.null(shape)) {
      tip_colors_df <- data.frame(
        label = sapply(species_dict, function(entry) entry$tip_label),
        species_colors = species_names,
        color = color[match(sapply(species_dict, function(entry) entry$species), unique(species_names))]
      )

      x <- full_join(tree_df, tip_colors_df, by = 'label')
      plot <- plot + geom_star(aes(x=x, subset = isTip, fill = species_colors, starshape = '15'), size = 3, show.legend = TRUE)
    }
  
    if ( is.null(color) && !is.null(shape)) {
      tip_shapes_df <- data.frame(
        label = sapply(species_dict, function(entry) entry$tip_label),
        species_shapes = species_names,
        shape = shape[match(sapply(species_dict, function(entry) entry$species), unique(species_names))]
      )
      
      x <- full_join(tree_df, tip_shapes_df, by = 'label')
      plot <- plot + geom_star(aes(x=x, subset = isTip, fill = species_colors, starshape = shape), size = 3, show.legend = TRUE)
    }
  }

    if ( bootstrap_circles ) {
      bootstrap_legend = 'Legend'
      
      # Include bootstrap values as circles 
      bs_values <- bootstrap_circles(tree)

      # Represent bootstrap values with colored circles
      if ( !is.null(bs_values) ) {
        x <- full_join(x, bs_values, by = 'node')
        bootstrap_colors <- c('(75,100]' = "black", '(50,75]' = "grey", '(0,50]' = "snow2")
    
        root <- rootnode(tree)
        plot <- plot + geom_point2(aes(subset = !isTip & node != root, fill = category), shape = 21, size = 2) + theme_tree(legend.position = c(0.8, 0.2)) + 
                       scale_fill_manual(values = c(color, bootstrap_colors), guide = bootstrap_legend, name = 'Bootstrap Support (BS)',
                                         breaks = c('(75,100]', '(50,75]', '(0,50]'), labels = expression("BS >=75", "50 <= BS < 75", "BS < 50")) 
      }
    }
  
    # Add clades and labels for specific leaves
    if (!is.null(clades) && !is.null(labels)) {
        plot <- plot + geom_cladelab(node = clades, label = labels, align = TRUE, fill = 'white',
                                     offset.text = labeldist, barsize = 0.9, offset.bar = 0.5, fontsize = fontsize)
    }
  
    if (exists("plot")) {
      return(plot)
    }
 }


tree <- read_tree("Sandfly_1275_P450s_plus_Agambiae.nwk")
cyp6acj <- extract_subtree(tree, "arabicus_TRINITY_DN1067_c1_g4_i1_p1_501", "schwetzi_TRINITY_DN7802_c0_g1_i5_p1_511", bl=TRUE)
visualize_tree(cyp6acj)
