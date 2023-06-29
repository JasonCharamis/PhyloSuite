## Library with functions related to plotting.

package_install <- function ( package_name ) {
  if (requireNamespace(package_name, quietly = TRUE)) {
    library(package_name, character.only = TRUE)
  }
  
  else {
    print ( (sprintf("%s %s",package_name, "is not installed. Installing it!")))

    if ( BiocManager::available(package_name) ) {
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

draw_pairs <- function (fig1, fig2, horizontal = "F", ...) {
  package_install("gridExtra")
  
  if ( horizontal == "T") {
    multi_panel_figure <- grid.arrange( fig1, fig2, nrow = 2, ncol = 1)
  }
  
  else {
    multi_panel_figure <- grid.arrange(fig1, fig2, nrow = 1, ncol = 2)
  }
  return (multi_panel_figure)
}
