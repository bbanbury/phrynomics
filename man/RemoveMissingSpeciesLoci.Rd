\name{RemoveMissingSpeciesLoci}
\alias{RemoveMissingSpeciesLoci}
\title{RemoveMissingSpeciesLoci}
\description{This function will remove entire loci (including multiple SNPs) when a species is missing data}
\usage{
RemoveMissingSpeciesLoci(data)
}
\arguments{
  \item{data}{A SNP dataset}
}
\details{For SNAPP analyses (maybe more) an entire species can not be missing loci. If all of the individuals within a species are missing data, then the entire locus is dropped from the dataset.}
\value{Returns a subset data frame}
\examples{
  data(fakeData)
  RemoveMissingSpeciesLoci(fakeData)
}
\author{B.Banbury}
