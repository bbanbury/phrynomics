% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReturnAmbyCode.R
\name{ReturnAmbyCode}
\alias{ReturnAmbyCode}
\title{Return Ambiguity Codes}
\usage{
ReturnAmbyCode(bases)
}
\arguments{
\item{bases}{A vector of bases}
}
\value{
Returns a single character with the ambiguity code.
}
\description{
This function will take a set of bases and return the IUPAC ambiguity code.
}
\examples{
bases <- c("A", "G")
ReturnAmbyCode(bases)
ReturnAmbyCode(c("A","G","C"))
}
\seealso{
\link{ReadSNP} \link{WriteSNP} \link{ReturnNucs}
}
