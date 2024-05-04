source("~/bin/EvoSuite/OliveTree/R/OliveTree2.R")
setwd("~/Documents/OliveTrees")

    
# Squamate species tree with 6,885 species from Title & Singhal et al. (2024). The macroevolutionary singularity of snakes. Science. 383(6685), 918–923.

# First read and visualize the whole tree with 6,885 species
squamate_tree <- read_tree(input_file = "squamates_Title_Science2024_ultrametric_constrained.tre")

visualize_tree (
  tree = "squamates_Title_Science2024_ultrametric_constrained.tre",
  save = FALSE
)

# Load file with full species taxonomy
squamates <- read.csv("squamateTreetaxa.csv")

# Highlight snakes and other major taxonomic categories from huge tree
snake_families <- c("Acrochordidae", "Aniliidae", "Anomalepididae", "Anomochilidae",
      "Bipedidae", "Boidae", "Bolyeridae", "Colubridae", "Cylindrophiidae",
      "Elapidae", "Gerrhopilidae", "Helodermatidae", "Homalopsidae",
      "Lamprophiidae", "Leptotyphlopidae", "Loxocemidae", "Pythonidae", "Typhlopidae", 
      "Uropeltidae", "Viperidae", "Xenopeltidae", "Xenophidiidae", "Xenotyphlopidae")

lizard_families <- c("Agamidae", "Alopoglossidae", "Amphisbaenidae", "Anguidae", "Anniellidae", "Anomalepididae", "Anomochilidae",
       "Bipedidae", "Blanidae", "Bolyeridae", "Cadeidae", "Carphodactylidae", "Chamaeleonidae", "Cordylidae", "Corytophanidae",
       "Crotaphytidae", "Cyclocoridae", "Dactyloidae", "Dibamidae", "Diplodactylidae", "Diploglossidae", "Eublepharidae", "Gekkonidae",
       "Gerrhopilidae", "Gerrhosauridae", "Gymnophthalmidae", "Helodermatidae", "Hoplocercidae", "Iguanidae", "Lacertidae", "Lanthanotidae",
       "Leiocephalidae", "Leiosauridae", "Liolaemidae", "Opluridae", "Pareidae", "Phrynosomatidae", "Phyllodactylidae", "Polychrotidae", "Pygopodidae",
       "Rhineuridae", "Scincidae", "Shinisauridae", "Sphaerodactylidae", "Teiidae", "Trogonophidae", "Tropiduridae", "Varanidae", "Xantusiidae", "Xenosauridae")


snake_species <- unlist(sapply(
  snake_families,
  function(family) {
    squamates[grepl(
      family, squamates[, 'family']
    ), ][, 1]
  }
))

lizard_species <- unlist(sapply(
  lizard_families,
  function(family) {
    squamates[grepl(
      family, squamates[, 'family']
    ), ][, 1]
  }
))


# Identified nodes corresponding to each category
print_internal_nodes(squamate_tree, form = "circular", references = crotalus_species)

highlight_nodes = c(
  "Snakes" = 6890, 
  "Crotalus" = 8340,
  "Lizards" = 8983
)

# Generate the highlighted tree with colors corresponding to the different clans
highlight_tree(
  squamate_tree,
  form = "circular",
  highlight_nodes = highlight_nodes,
  colors = c(
    "Crotalus" = "blue",
    "Lizards" = "green4"
  ),
  save = FALSE
)

# Extract Crotalus genus subtree
crotalus_species <- c(squamates[grepl('Crotalus', squamates[, 'genus']), ][, 1])
crotalus_subtree <- extract_subtree(tree = squamate_tree, taxon_list = crotalus_species)

# Generate random colors for each species
group_colors <- sapply ( crotalus_subtree$tip.label, function (x) { c( x = sample(colors(), 1) ) } )
names(group_colors) <- sapply( names(group_colors), function (x) { gsub(".x", "", x) } )

visualize_tree( 
  crotalus_subtree, 
  color = group_colors, 
  form = "circular", 
  taxonID_separator = NULL,
  tiplabels = FALSE, 
  tip_shape_size = 4, 
  branch_length = TRUE, 
  save = FALSE
)
