source("~/bin/EvolSuite/OliveTree/R/OliveTree.R")
#source("~/Documents/test/EvolSuite/OliveTree/R/OliveTree.R")
setwd("~/Documents/OliveTrees")

# Squamate species tree with 6,885 species from Title & Singhal et al. (2024). The macroevolutionary singularity of snakes. Science. 383(6685), 918–923.

# First read and visualize the whole tree with 6,885 species
squamate_tree <- read_tree(input_file = "squamates_Title_Science2024_ultrametric_constrained.tre")

'''
visualize_tree (
  tree = "squamates_Title_Science2024_ultrametric_constrained.tre",
  save = FALSE
)
'''

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

snake_species_tips <- unlist(sapply(
  squamate_tree@phylo$tip.label,
  
  function (species) {
    grep(species, snake_species, value = TRUE, fixed = TRUE)
  }
))

lizard_species_mappings <- setNames(rep("red", length(lizard_species_tips)), lizard_species_tips)
snake_species_mappings <- setNames(rep("red", length(snake_species_tips)), snake_species_tips)

visualize_tree(
  tree = squamate_tree,
  form = "circular",
  tiplabels = TRUE,
  pattern_id = names(snake_species_mappings),
  tip_label_color = snake_species_mappings,
  save = FALSE
)


# Identified nodes corresponding to each category
#print_internal_nodes(squamate_tree, form = "circular", references = crotalus_species)

highlight_nodes = c(
  "Snakes" = 6890, 
  "Lizards" = 8983
)

# Generate the highlighted tree with colors corresponding to the different clans
highlight_tree(
  squamate_tree,
  form = "circular",
  highlight_nodes = highlight_nodes,
  colors = c(
    "Lizards" = "green4",
    "Snakes" = "red"
  ),
  legend_direction = "vertical",
  legend_key_size = 0.8,
  legend_font_bold = "plain",
  legend_spacing_y = 1,
  legend_key_width = 0.8,
  save = FALSE
)

# Extract Crotalus genus subtree
crotalus_species <- c(squamates[grepl('Crotalus', squamates[, 'genus']), ][, 1])
crotalus_subtree <- extract_subtree(tree = squamate_tree, taxon_list = crotalus_species)

# Generate random colors for each species
group_colors <- setNames(sample(colors(), length(crotalus_subtree$tip.label)), crotalus_subtree$tip.label)
group_colors

visualize_tree( 
  crotalus_subtree, 
  color = group_colors, 
  form = "circular", 
  legend_position = "right",
  legend_text_position = "right",
  legend_orientation = "vertical",
  legend_key_size = 0.000001,
  legend_font = "Arial",
  legend_fontsize = 8.5,
  legend_font_bold = "italic",
  legend_spacing_x = 0.00001,
  legend_spacing_y = 0.00001,
  legend_key_width = 0.00001,
  legend_title_hjust = 0.0001,
  legend_name = "Species",
  tiplabels = FALSE, 
  taxon_group_separator = NULL,
  tip_label_color = NULL,
  tip_shape_size = 5, 
  branch_length = TRUE, 
  save = FALSE
)


