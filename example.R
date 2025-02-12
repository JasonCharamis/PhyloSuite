# Load required libraries
library(ape)
library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)

# Create example tree
set.seed(123)
tree <- rtree(10)

# Create example data for pie segments
tip_data <- data.frame(
  tip_label = tree$tip.label,
  condition_A = runif(10),
  condition_B = runif(10),
  condition_C = runif(10)
)

# Normalize the proportions
tip_data[,2:4] <- t(apply(tip_data[,2:4], 1, function(x) x/sum(x)))

# Convert to long format for mapping
tip_data_long <- tip_data %>%
  pivot_longer(
    cols = starts_with("condition"),
    names_to = "condition",
    values_to = "proportion"
  ) %>%
  group_by(tip_label) %>%
  # Calculate cumulative proportions for star segments
  mutate(
    start_angle = lag(cumsum(proportion), default = 0) * 2 * pi,
    end_angle = cumsum(proportion) * 2 * pi
  )

# Create the plot
p <- ggtree(tree) +
  geom_star(
    data = tip_data_long,
    aes(
      x = x,  # You'll need to join with tip coordinates
      y = y,
      starshape = condition,
      fill = condition,
      angle = start_angle,
      angle2 = end_angle
    ),
    size = 3,
    starstroke = 0.1
  ) +
  scale_fill_manual(
    values = c(
      "condition_A" = "#E41A1C",
      "condition_B" = "#377EB8",
      "condition_C" = "#4DAF4A"
    )
  ) +
  theme(legend.position = "right")

# Function to create tree with star pie charts
create_star_pie_tree <- function(tree, data, conditions, 
                               colors = NULL, star_size = 3) {
  # Validate input
  if(!all(conditions %in% colnames(data))) {
    stop("All conditions must be present in data")
  }
  
  # Get base tree plot
  p <- ggtree(tree)
  
  # Get tip coordinates
  tip_coords <- p$data %>%
    filter(isTip) %>%
    select(label = label, x, y)
  
  # Prepare data for star segments
  data_long <- data %>%
    pivot_longer(
      cols = all_of(conditions),
      names_to = "condition",
      values_to = "proportion"
    ) %>%
    group_by(label) %>%
    mutate(
      start_angle = lag(cumsum(proportion), default = 0) * 2 * pi,
      end_angle = cumsum(proportion) * 2 * pi
    )
  
  # Join with coordinates
  plot_data <- left_join(data_long, tip_coords, by = "label")
  
  # Set default colors if none provided
  if(is.null(colors)) {
    colors <- setNames(
      rainbow(length(conditions)),
      conditions
    )
  }
  
  # Create plot with star pie charts
  p + 
    geom_star(
      data = plot_data,
      aes(
        x = x,
        y = y,
        starshape = condition,
        fill = condition,
        angle = start_angle,
        angle2 = end_angle
      ),
      size = star_size,
      starstroke = 0.1
    ) +
    scale_fill_manual(values = colors) +
    theme(legend.position = "right")
}

# Example usage:
# my_tree <- read.tree("tree.nwk")
# my_data <- read.csv("proportions.csv")
# conditions <- c("typeA", "typeB", "typeC")
# colors <- c(
#   "typeA" = "red",
#   "typeB" = "blue",
#   "typeC" = "green"
# )
# tree_plot <- create_star_pie_tree(
#   my_tree,
#   my_data,
#   conditions,
#   colors,
#   star_size = 3
# )