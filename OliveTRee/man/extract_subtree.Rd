% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OliveTRee.R
\name{extract_subtree}
\alias{extract_subtree}
\title{extract_subtree
Function to extract a subtree by finding the MRCA (Most Recent Common Ancestor) of two anchor leaves or
the subtree that contains a user-provided list of tip labels,  while preserving branch lengths and bootstrap values.}
\usage{
extract_subtree(tree, tip1 = NULL, tip2 = NULL, taxon_list = NULL)
}
\arguments{
\item{tree}{A phylogenetic tree object of class 'phylo' or 'treedata'.}

\item{tip1, tip2}{The labels or pattern of the two anchor tips.}

\item{taxon_list}{User-provided list of tip labels, the subtree which contains them will be extracted.}
}
\value{
A treedata object representing the extracted subtree with preserved branch lengths and bootstrap values.
}
\description{
extract_subtree
Function to extract a subtree by finding the MRCA (Most Recent Common Ancestor) of two anchor leaves or
the subtree that contains a user-provided list of tip labels,  while preserving branch lengths and bootstrap values.
}
\examples{
# Extract a subtree preserving branch lengths and bootstrap values
tree <- read.tree("path/to/your/treefile.newick")
subtree <- extract_subtree(tree, "species_A", "species_B", branch_length = TRUE)
visualize_tree(subtree)

# Extract a cladogram with equal branch lengths
tree <- read.tree("path/to/your/treefile.newick")
subtree <- extract_subtree(tree, "tip1", "tip")
visualize_tree(subtree)

}
\seealso{
\code{\link{visualize_tree}} for visualizing phylogenetic trees.
\code{\link{findMRCA}} for finding the Most Recent Common Ancestor (MRCA) of multiple nodes.
}
