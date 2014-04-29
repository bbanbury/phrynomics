\name{ReturnUniqueBases}
\alias{ReturnUniqueBases}
\title{Return Unique Bases}
\description{This function will return the unique nucleotide bases of an SNP accounting for an ambiguity code}
\usage{
ReturnUniqueBases(SNP)
}
\arguments{
  \item{SNP}{A single SNP}
}
\details{See http://www.bioinformatics.org/sms/iupac.html}
\value{Returns vector of bases
}
\author{B.Banbury}
\seealso{
  \code{\link{ReturnAmbyCode}}
  \code{\link{ReturnMrBayesAmbyCode}}
}
\examples{
ReturnUniqueBases(c("R", "Y", "W", "D", "S", "M", "K"))
ReturnUniqueBases(c("A"))
}