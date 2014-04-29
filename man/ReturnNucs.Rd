\name{ReturnNucs}
\alias{ReturnNucs}
\title{Return Nucleotide Bases}
\description{This function will return possible nucleotides given an IUPAC ambiguity code}
\usage{
ReturnNucs(NucCode, forSNAPP=FALSE)
}
\arguments{
  \item{NucCode}{A single SNP Base or ambiguity code}
  \item{forSNAPP}{Whether to trat missing characters as possible bases (FALSE), or as missing data (TRUE), see example below}
}
\details{See \link{http://www.bioinformatics.org/sms/iupac.html}}
\value{Returns possible nucleotides given an ambiguity code. 
}
\author{B.Banbury}
\seealso{
  \code{\link{ReturnAmbyCode}}
  \code{\link{ReturnMrBayesAmbyCode}}
}
\examples{
ReturnNucs("R")
ReturnNucs("N", forSNAPP=FALSE)
ReturnNucs("N", forSNAPP=TRUE)
}