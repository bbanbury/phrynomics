\name{RemoveInvariantSites}
\alias{RemoveInvariantSites}
\title{Remove Invariant SNP Sites}
\description{This function will remove sites that have no variation.}
\usage{
RemoveInvariantSites(SNPdataset, chatty=FALSE)
}
\arguments{
  \item{SNPdataset}{A SNP dataset}
  \item{chatty}{Optional print to screen update}
}
\details{Ambiguity codes are taken into account. A SNP with a single base plus an ambiguity code that includes that base (for example, AAM) will be dropped, because M can be either A or C.}
\value{Returns a subset of the original data frame
}
\examples{
  data(fakeData)
   RemoveInvariantSites(fakeData)
}
\author{B.Banbury}
