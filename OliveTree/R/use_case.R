
source ("~/bin/OliveTrees/OliveTree/R/OliveTree.R")
setwd("~/Documents/OliveTrees")

# Squamate species tree with 6,885 species from Title & Singhal et al. (2024). The macroevolutionary singularity of snakes. Science. 383(6685), 918–923. 

# First read and visualize the whole tree with 6,885 species
squamate_tree <- read_tree(input_file = "squamates_Title_Science2024_ultrametric_constrained.tre")

visualize_tree( 
  tree = squamate_tree,
  form = "circular",
  tiplabels = TRUE,
  tip_label_size = 1,
  bootstrap_numbers = TRUE,
  bootstrap_circles = FALSE,
  bootstrap_number_nudge_x = 0,
  bootstrap_number_nudge_y = 0.4,
  save = FALSE
)


# Show tip labels with Crotalus
visualize_tree(
  tree = squamate_tree,
  form = "circular",
  tiplabels = TRUE,
  pattern_id = "Crotalus",
  save = TRUE,
  output = "whole_squamate_tree_with_Crotalus.svg"
)


# Load file with full species taxonomy
squamates <- read.csv("squamateTreetaxa.csv")

# Highlight snakes and other major taxonomic categories from huge tree
snake_families <- c("Acrochordidae","Aniliidae","Anomalepididae","Anomochilidae",
                    "Bipedidae","Boidae","Bolyeridae","Colubridae","Cylindrophiidae",
                    "Elapidae","Gerrhopilidae","Helodermatidae","Homalopsidae",
                    "Lamprophiidae","Leptotyphlopidae","Loxocemidae","Pythonidae",
                    "Typhlopidae","Uropeltidae","Viperidae","Xenopeltidae","Xenophidiidae","Xenotyphlopidae")

lizard_families <- c("Agamidae","Alopoglossidae","Amphisbaenidae","Anguidae","Anniellidae","Anomalepididae","Anomochilidae",
                     "Bipedidae","Blanidae","Bolyeridae","Cadeidae","Carphodactylidae","Chamaeleonidae","Cordylidae","Corytophanidae",
                     "Crotaphytidae","Cyclocoridae","Dactyloidae","Dibamidae","Diplodactylidae","Diploglossidae","Eublepharidae","Gekkonidae",
                     "Gerrhopilidae","Gerrhosauridae","Gymnophthalmidae","Helodermatidae","Hoplocercidae","Iguanidae","Lacertidae","Lanthanotidae",
                     "Leiocephalidae","Leiosauridae","Liolaemidae","Opluridae","Pareidae","Phrynosomatidae","Phyllodactylidae","Polychrotidae","Pygopodidae",
                     "Rhineuridae","Scincidae","Shinisauridae","Sphaerodactylidae","Teiidae","Trogonophidae","Tropiduridae","Varanidae","Xantusiidae","Xenosauridae")


snakes <- unlist(sapply(
    snake_families,
    
    function (family) {
      squamates[grepl(
        family, squamates[,'family']),][,1]
    }
))

lizards <- unlist(sapply(
  lizard_families,
  
  function (family) {
    squamates[grepl(
      family, squamates[,'family']),][,1]
  }
))


print_internal_nodes(squamate_tree, form = "circular", references = c(snakes, lizards))

highlight_nodes = c("Snakes" = 6890,
                    "Lizards" = 8983)

# Generate the highlighted tree with colors corresponding to the different clans
highlight_tree(
  squamate_tree,
  form = "circular",
  highlight_nodes = highlight_nodes, 
  colors = c("Snakes"= "red3",
             "Lizards" = "green4"
  ),
  save = TRUE,
  output = "Squamate_tree_highlighted.svg"
)


# Extract Crotalus genus subtree
crotalus_species <- c(squamates[grepl('Crotalus', squamates[,'genus']),][,1])
crotalus_subtree <- extract_subtree(tree = squamate_tree, taxon_list = crotalus_species)

crotalus_species

visualize_tree(
  crotalus_subtree,
  form = "circular",
  tiplabels = TRUE,
  tip_label_size = 4,
  save = TRUE,
  output = "Crotalus_tree.svg"
)

