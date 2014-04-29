\name{GetAncestors}
\alias{GetAncestors}
\title{Get Ancestors of a single tip taxon}
\description{This function will return a vector of nodes throughout a tree for a single taxon}
\usage{
GetAncestors(treeEdgeMatrix, tip)
}
\arguments{
  \item{treeEdgeMatrix}{An edge matrix from a tree, such as tree$edge}
  \item{tip}{An single tip node number (not taxon name)}
}
\details{}
\value{Returns a vector of nodes in a tree.
}
\author{B.Banbury}
