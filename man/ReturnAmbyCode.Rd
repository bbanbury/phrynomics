\name{ReturnAmbyCode}
\alias{ReturnAmbyCode}
\title{Return Ambiguity Code}
\description{This function will take a set of bases and return the IUPAC ambiguity code}
\usage{
ReturnAmbyCode(bases)
}
\arguments{
  \item{bases}{A vector of bases}
}
\details{See \link{http://www.bioinformatics.org/sms/iupac.html}}
\value{Returns a single character with IUPAC code.
}
\author{B.Banbury}
\seealso{
  \code{\link{ReturnNucs}}  
  \code{\link{ReturnMrBayesAmbyCode}}
}
\examples{
ReturnAmbyCode(c("A", "T"))
ReturnAmbyCode(c("G","C","T"))
}