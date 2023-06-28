
## Library of functions for advanced tree manipulation and visualization using ggtree.

package_install <- function ( package ) {
  if (requireNamespace(package, quietly = TRUE)) {
      library(package, character.only = TRUE)
  }
  
  else {
    print ( (sprintf("%s %s",package, "is not installed. Installing it!")))
    BiocManager::install(package)  
  }
}


## load dependencies - if not present will install them
dependencies <- c("ape","phytools","treeio","TreeTools","ggstar","ggtree","ggplot2","dplyr","stringi","stringr")

for ( i in dependencies ) { 
  package_install(i)
  }


## read and order tree 
read_tree <- function ( nwk ) {
  tree <- ape::read.tree(nwk)
  t <- Preorder(tree)
  return (t)
}


## print tree with node ids - very useful for multiple tasks
node_ids <- function ( tree ) {
  ggtree(tree, layout="circular") + geom_text(aes(label = node), size = 3)
}

## collapse nodes based on bootstrap support - returns a "phylo" object
bootstrap_collapse <- function ( tree, cutoff ) { 
  return ( as.polytomy(tree, feature='node.label', fun=function(x) as.numeric(x) < cutoff) )
  }


#flip based on descendant nodes (internal nodes or leaves if node is terminal)
flip_node <- function (tree, node1, node2) {
      flipped_tree <- flip(ggtree(tree), node1, node2)
      return (as.phylo(flipped_tree))
}


## group all descendant branches of node(s)
group_descendants <- function (tree, node1, node2="", node3="",node4="") {
    tree <- groupClade(tree, .node=c(node1, node2, node3, node4))
    return ( tree )
}


## extract subtree by finding the MRCA of two anchor nodes while preserving branch lengths and bootstrap values 
extract_subtree <- function (t, tip1, tip2) {
  subtree <- ape::extract.clade(t,node = phytools::findMRCA(t,c(tip1, tip2)))
  return ( subtree )
}


## include bootstrap values as colored circles
bootstrap_circles <- function ( x ) {
  
  bootstrap_colors <- c('(75,100]' = "black",
                        '(50,75]' = "grey",
                        '(0,50]' = "snow2")
  
  categories <- c('(75,100]' ,'(50,75]','(0,50]')
  
  bs_tibble <- tibble( node=1:Nnode(x@phylo) + Ntip(x@phylo), bootstrap = as.numeric(x@phylo$node.label), cut(bootstrap,c(0, 50, 75, 100))) ## create dataframe with bootstrap support for each node
  p <- data.frame(node=bs_tibble$node, bootstrap=as.numeric(bs_tibble$bootstrap), category = cut(bs_tibble$bootstrap, c(0, 50, 75, 100)), 
                  b_color = bootstrap_colors[match( cut(bs_tibble$bootstrap, c(0, 50, 75, 100)), categories)] ) ## assign each node to a support category with associated color
  return (p)
}

## visualize tree with a wide variety of options and customizable features
draw_single_tree <- function ( tree, node1=NULL, node2=NULL, node3=NULL, reference1=NULL, reference2=NULL, bootstrap_legend=TRUE, ... ) {
  
  ## if species matches reference discard from rest of tip labels - at least two references, adjust for more
  if (!is.null(reference2)) {
      reference <- c(reference1, reference2)
      ref_idx <- grep(paste(reference, collapse = "|"), tree$tip.label) 
  }
  
  if (is.null(reference2) & !is.null(reference1)) {
      reference <- reference1
      ref_idx <- grep(paste(reference), tree$tip.label) 
  }
  
  ref_labs <- tree$tip.label[ref_idx]
  ref_species <- data.frame(node=ref_idx, name= c(sub("_.*","",ref_labs)), label=ref_labs)
  species <- sub("_.*","",tree$tip.label[-ref_idx])
  
  if (is.null(reference1) & is.null(reference2)) {
      ref_labs <- NULL
      ref_species <- NULL
      ref_idx <- NULL
      species <- sub("_.*","",tree$tip.label)
  }
  
  ## get index of identified species
  species_idx <- grep(paste(species, collapse = "|"), tree$tip.label)
  
  ## else, keep name for tip color and shape mapping per species
  unique_species <- sort(unique(species))
  
  ## customize tip colors and shapes 
  ## ideally species names should be sorted to always make correct color and shape mappings
  group_colors <- c(arabicus="purple",
                    argentipes= "cyan",
                    duboscqi="salmon2",
                    longipalpis="darkgreen",
                    migonei= "green",
                    orientalis= "gold",
                    papatasi= "darkred",
                    perniciosus= "lightgreen",
                    schwetzi ="blue",
                    sergenti= "orange",
                    tobbi= "red")
  
  ## shape mappings
  tip_label_shape <- c(arabicus="15",
                       argentipes="15",
                       duboscqi="15",
                       longipalpis="13", 
                       migonei="13",
                       orientalis="15",
                       papatasi="15",
                       perniciosus="15",
                       schwetzi="10",
                       sergenti="15",
                       tobbi="15"
  )
  
  
  ## create associative dataframes of tip color and shape 
  tip_colors_df <- data.frame(label = tree$tip.label[species_idx], species=species, colour = group_colors[match(species, unique_species)])
  tip_shapes_df <- data.frame(label = tree$tip.label[species_idx], species=species, shape=tip_label_shape[match(species, unique_species)])

  ## join tip colors and shapes in tree object based on tip label
  x <- data.frame()
  x <- full_join(tree, tip_colors_df, by='label')
  x <- full_join(x, tip_shapes_df, by='label')
  
  ## include bootstrap values as circles 
  bs_values <- bootstrap_circles (x)
  x <- full_join(x, bs_values, by='node')
  
  bootstrap_colors <- c('(75,100]' = "black",
                        '(50,75]' = "grey",
                        '(0,50]' = "snow2")
  
  node_ids(tree) ## print the node ids to select which nodes to highlight
  
  nodes <- as.character(x@data$node)
  
  if (bootstrap_legend == TRUE) {
      bootstrap_legend = 'legend'
  }
  
  else {
      bootstrap_legend = "none"
  }
  
  if ( !is.null(node3) ) {
    highlights <- c("olivedrab2","rosybrown1","wheat1")
  }
  
  else {
    highlights <- c("rosybrown1","wheat1")
  }
  
  root <- rootnode(tree)
  ## draw tree with bootstrap nodes, color and shape mappings, and highlighted nodes
  ggtree(x) + geom_star(aes(x,subset=isTip,starshape=shape, fill=species.x), size=3, show.legend = F) +
    geom_point2(aes(subset=!isTip & node != root, fill=category), shape=21, size=2) + theme_tree(legend.position=c(0.8, 0.2)) +
    scale_fill_manual(values=c(group_colors,bootstrap_colors), guide=bootstrap_legend, name='Bootstrap Support (BS)',
                      breaks=c('(75,100]', '(50,75]', '(0,50]'), labels=expression("BS >=75", "50 <= BS < 75", "BS < 50")) +
    geom_highlight(node=c( node1, node2, node3), fill=highlights, alpha=0.3, extend=0.10, linetype=23, colour="black") +
    geom_text(aes(label = ifelse(nodes %in% ref_species$node, ref_species$name, "")), 
              position = position_nudge(x = 0.08), vjust = 0.5, hjust=0.9, size = 3.5, family="Arial", fontface="bold")
}


draw_pairs <- function (fig1, fig2, horizontal = F, ...) {
  package_install("gridExtra")
  
  if ( horizontal == T) {
    multi_panel_figure <- grid.arrange( fig1, fig2, nrow = 2)
  }
  
  else {
    multi_panel_figure <- grid.arrange(fig1, fig2, nrow = 1)
  }
}


