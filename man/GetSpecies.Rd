% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetSpecies.R
\name{GetSpecies}
\alias{GetSpecies}
\title{Get Species Names in a Dataset}
\usage{
GetSpecies(taxa)
}
\arguments{
\item{taxa}{A character vector of taxon names. Numbers in species names represent individuals of the same species, for example "taxon1" and "taxon2" are both the same species "taxon".}
}
\value{
Returns a vector of unique species names
}
\description{
This function will detect the names of species in a dataset by removing the number identifyer.
}
\examples{
taxa <- c("taxon1", "taxon2", "taxon3")
GetSpecies(taxa)

data(fakeData)
GetSpecies(rownames(fakeData))
}
\seealso{
\link{ReadSNP} \link{IsMissingSpecies} \link{RemoveMissingSpeciesLoci}
}
