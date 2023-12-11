## Library with functions related to plotting.

package_install <- function ( package_name ) {
  if (requireNamespace(package_name, quietly = TRUE)) {
    library(package_name, character.only = TRUE)
  }
  
  else {
    print ( (sprintf("%s %s",package_name, "is not installed. Installing it!")))

    if ( package_name %in% BiocManager::available() ) {
      BiocManager::install(package_name)
    }
    
    else {
      install.packages(package_name)
    }
    
  }
}

dependencies <- c("ggplot2","ggplotify")

for ( i in dependencies ) { 
  package_install(i)
}


multipanel <- function ( horizontal = "T", ... ) {
  package_install("gridExtra")

  plots <- list(...)
  
  ## arrange horizontally
  if ( horizontal == "T") {
    multi_panel_figure <- grid.arrange(grobs = plots, ncol = length(plots), nrow = 1)
  }
  
  ## arrange vectically
  else {
    multi_panel_figure <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)
  }
  return (multi_panel_figure)
}
