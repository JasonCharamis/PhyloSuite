#' Library of plotting functions.
#' @import ggplot2
#' @import ggplotify
#' @import gridExtra
#'
#' @title load_packages: Function to check if a package is installed, and if not, install it.
#' @description This function checks if a package is installed. If not, it installs the package using BiocManager if available, otherwise using install.packages.
#' @param tools A character vector of package names to be checked and installed.
#' @return NULL
#' @export

# Function to check if a package is installed, and if not, install it.
# Function to install or load a package
load_packages <- function( tools ) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built, 1, 1)))
  tmp <- tmp[as.numeric(substr(tmp$Built, 1, 1)) == max_version, ] 
  
  for ( pkg in tools ) {
    if ( pkg %in% tmp$Package ) {
      library (pkg, character.only = TRUE)
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

dependencies <- c("ggplot2", "ggplotify", "gridExtra")

load_packages(dependencies)


multipanel <- function ( horizontal = "T", ... ) {
  plots <- list(...)
  
  if ( horizontal == "T") { # Arrange horizontally
    multi_panel_figure <- grid.arrange(grobs = plots, ncol = length(plots), nrow = 1)
  } else { # Arrange vertically
    multi_panel_figure <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)
  }
  return (multi_panel_figure)
}
