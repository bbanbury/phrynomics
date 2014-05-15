\name{removeOutgroups}
\alias{removeOutgroups}
\title{Remove Outgroups}
\description{This function will remove outgroups from the dataset}
\usage{
removeOutgroups(SNPdataset, taxon="PH")
}
\arguments{
  \item{SNPdataset}{A SNP dataset}
  \item{taxon}{A tag for the ingroup}
}
\details{This function will keep data in which the taxon is included in the rowname, otherwise it is assumed to be an outgroup.}
\value{Returns a subset data frame.}
\examples{
  data(fakeData)
  removeOutgroups(fakeData, taxon="in")
}
\author{B.Banbury}
