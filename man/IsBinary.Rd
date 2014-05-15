\name{IsBinary}
\alias{IsBinary}
\title{Is Binary SNP}
\description{This function will detect whether a SNP is binary}
\usage{
IsBinary(SNP)
}
\arguments{
  \item{SNP}{A single SNP}
}
\details{This function will determine if a SNP site is binary (TRUE) or not (FALSE). For example, c("A", "A", "T") will return TRUE; c("A", "A", "A") or c("A", "G", "T") will return FALSE 

Ambiguity codes are taken into account to mean either heterozygotes or uncertainty. Such that c("A", "A", "R") will return TRUE; c("A", "A", "Y") will return FALSE}
\value{Boolean return for whether SNP is binary
}
\examples{
  data(fakeData)
  splitData <- SplitSNP(fakeData)  #the IsVariable function will only work on single SNPs
  IsBinary(splitData[,1])
  IsBinary(splitData[,2])
  apply(splitData, 2, IsBinary)
}
\author{B.Banbury}
