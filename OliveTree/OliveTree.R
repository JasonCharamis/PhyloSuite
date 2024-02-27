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
read_tree <- function(newick) {
  t <- phytools::read.newick(newick)
  return(Preorder(t))
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
bootstrap_collapse <- function(tree, cutoff = 0.5) { 
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
  t <- Preorder(tree)
  if (!(bl == "TRUE" || bl == "T")) {
    t$edge.length <- NULL
  } 
  return(ape::extract.clade(t, node = phytools::findMRCA(t, c(tip1, tip2))))
}

#==================================== TREE VISUALIZATION FUNCTIONS ====================================#

bootstrap_circles <- function(tree_obj) {
  if (!is.null(tree_obj$node.label)) {
    bootstrap_colors <- c('(75,100]' = "black", '(50,75]' = "grey", '(0,50]' = "snow2")
    categories <- c('(75,100]', '(50,75]', '(0,50]')
    
    internal_nodes <- phytools::drop.leaves(tree_obj)
    
    bs_tibble <- tibble(
      node = 1:tree_obj$Nnode,
      bootstrap = as.numeric(tree_obj$node.label),
      cut(bootstrap, c(0, 50, 75, 100)))
    
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


visualize_tree <- function(tree, form = "rectangular", tiplabels = TRUE, pattern=NULL,
                           bootstrap_circles = TRUE, bootstrap_numbers = FALSE,
                           color = NULL, shape = NULL,
                           clades = NULL, labels = NULL, mappings_legend = TRUE, 
                           save = TRUE, output = NULL, ...) {
  
  # Open phylogenetic tree file and/or object
  if (any(grepl(".newick|.nwk|.tre|.support", tree))) {
    tree_obj <- read_tree(tree)
  } else if (is(tree,"phylo")) {
    tree_obj <- tree
  } else {
    print ("Provided file is not a newick file or a phylo object.")
    return()
  }
  
  # Manipulate bootstrap values
  if (bootstrap_circles == TRUE && bootstrap_numbers != TRUE ) {
    bootstrap_legend <- 'Legend'
    root <- rootnode(tree_obj)
    bs_values <- bootstrap_circles(tree_obj)
    tree_with_bs <- full_join(as_tibble(tree_obj), bs_values, by = 'node')
    
    plot <- ggtree(tree_with_bs, layout = form) + geom_point(data = bs_values, aes( node != isTip & node != root & !is.na(b_color), fill = category), shape = 21, size = 2) + 
                                                  scale_fill_manual(values=c("black", "grey", "snow2"), guide='legend', 
                                                  name='Bootstrap Support(BP)', 
                                                  breaks=c('(75,100]', '(50,75]', '(0,50]'), 
                                                  labels=expression(BP>=75, 50 <= BP * " < 75", BP < 50)
                                                  )
    } else if (bootstrap_circles != TRUE && bootstrap_numbers == TRUE) {
        plot <- ggtree(tree_obj, layout = form)
        plot <- plot + geom_nodelab()
    } else if (bootstrap_circles == TRUE && bootstrap_numbers == TRUE) {
        stop ("bootstrap_circles = TRUE (default) and bootstrap_numbers = TRUE. Please select either bootstrap circles or bootstrap numbers option.")
    } else {
        print ("Bootstrap values will not be printed.")
        plot <- ggtree(tree_obj, layout = form)
    }
  
   # Visualize phylogenetic tree with color and shape mappings
   if (is.null(color) && is.null(shape) && is.null(clades) && is.null(labels) ) {
     if (tiplabels == TRUE) {
       if (is.null(pattern)) {
         return(plot + geom_tiplab())
       } else {
           return(plot + geom_tiplab(aes(subset=grepl(pattern, label))))
       }
     } else {
         return (plot)
     } 
   } else {
       # Manipulate tip labels to generate species_names and isolate tip labels which match reference
       references <- list(...)
       reference_args <- names(references)
       
       if (length(reference_args) == 0) {
         species_dict <- lapply(seq_along(tree_obj$tip.label), function(i) {
           list(tip_label = tree_obj$tip.label[i], species = sub("_.*", "", tree_obj$tip.label[i]))
         })
       } else { 
           reference <- as.character(unlist(references))
           ref_species <- tree_obj$tip.label[which(sapply(tree_obj$tip.label, function(x) any(grepl(reference, x))))]
           ref_dict <- data.frame(
             tip_label = ref_species,
             reference_flag = TRUE
           )
        
           plot_ref <- plot %<+% ref_dict 
           plot <- plot_ref + geom_tiplab(aes(label = ifelse(reference_flag, label, "")), color = "black")
           
           non_ref_species <- setdiff(tree_obj$tip.label, ref_species)
           species_dict <- lapply(non_ref_species, function(label) {
                                  list(tip_label = label, species = sub("_.*", "", label))
                                  })
           
           species_names <- sapply(species_dict, function(entry) entry$species)

           if (!is.null(color) && !is.null(shape)) {
             missing_species_color <- setdiff(unique(sapply(species_dict, function(entry) entry$species)), names(color))
             if (length(missing_species_color) > 0) {
               warning(paste("Color mapping not found for the following species: ", paste(missing_species_color, collapse = ", ")))
             }
             
             # Create data frames for tip colors and shapes
             tip_colors_df <- data.frame(
               label = sapply(species_dict, function(entry) entry$tip_label),
               species_colors = species_names,
               s_color = color[match(sapply(species_dict, function(entry) entry$species), names(color))]
             )
             
             # Check for missing species in color dataframe
             missing_species_shape <- setdiff(unique(sapply(species_dict, function(entry) entry$species)), names(shape))
             if (length(missing_species_shape) > 0) {
               warning(paste("Shape mapping not found for the following species: ", paste(missing_species_shape, collapse = ", ")))
             }
             
             tip_shapes_df <- data.frame(
               label = sapply(species_dict, function(entry) entry$tip_label),
               species_shapes = species_names,
               s_shape = shape[match(sapply(species_dict, function(entry) entry$species), names(shape))]
             )
             
             # Join data frames to create a unique mapping tibble
             plot <- plot %<+% tip_colors_df 
             plot <- plot %<+% tip_shapes_df  
             
             plot <- plot + geom_star(mapping = aes( subset = isTip & !is.na(s_color) & !is.na(s_shape),
                                                     fill = ifelse(!is.na(species_colors), s_color, "black"),
                                                     starshape = ifelse(!is.na(species_shapes), s_shape, "circle")),
                                     size = 3, show.legend = mappings_legend) +
                                     scale_starshape_identity() + 
                                     scale_fill_identity()

          } else if ( !is.null(color) && is.null(shape)) {
              tip_colors_df <- data.frame(
                label = sapply(species_dict, function(entry) entry$tip_label),
                species_colors = species_names,
                s_color = color[match(sapply(species_dict, function(entry) entry$species), names(color))]
                )
              
              plot <- plot %<+% tip_colors_df
              plot <- plot + geom_star(mapping = aes( subset = isTip & !is.na(s_shape),
                                                      fill = ifelse(!is.na(species_colors), s_color, "black"),
                                                      starshape = ifelse(!is.na(species_shapes), s_shape, "circle")),
                                       size = 3, show.legend = mappings_legend) +
                                       scale_fill_identity()
              
          } else if ( is.null(color) && !is.null(shape)) {
              tip_shapes_df <- data.frame(
                label = sapply(species_dict, function(entry) entry$tip_label),
                species_shapes = species_names,
                s_shape = shape[match(sapply(species_dict, function(entry) entry$species), names(shape))]
                )
              
              plot <- plot %<+% tip_shapes_df
              plot <- plot + geom_star(mapping = aes( subset = isTip & !is.na(s_shape),
                                                      fill = ifelse(!is.na(species_colors), s_color, "black"),
                                                      starshape = ifelse(!is.na(species_shapes), s_shape, "circle")),
                                       size = 3, show.legend = mappings_legend) +
                                       scale_starshape_identity()
          }
        }
        
        if (!is.null(clades) && !is.null(labels)) {
          plot <- plot + geom_cladelab(node = clades, label = labels, align = TRUE, fill = 'white',
                                         offset.text = labeldist, barsize = 0.9, offset.bar = 0.5, fontsize = fontsize)
        } else if (!is.null(clades) && is.null(labels)) {
          plot <- plot + geom_cladelab(node = clades, label = "", align = TRUE, fill = 'white',
                                       offset.text = labeldist, barsize = 0.9, offset.bar = 0.5, fontsize = fontsize)
        } else if (is.null(clades) && !is.null(labels)) {
          plot <- plot + geom_cladelab(node = "", label = labels, align = TRUE, fill = 'white',
                                       offset.text = labeldist, barsize = 0.9, offset.bar = 0.5, fontsize = fontsize)
        }
   }    
       if (exists("plot")) {
         if (save == TRUE) {
           if (is.null(output)) {
             if (any(grepl(".newick|.nwk|.tre|.support", tree))) {
               ggsave(plot = plot, sprintf("%s_visualized.svg", sub(".newick|.nwk|.tre|.support", "", tree)), dpi = 600)
               message <- sprintf("Tree plotted and saved as %s_visualized.svg", sub(".newick|.nwk|.tre|.support", "", tree))
               print(message)
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
   }
