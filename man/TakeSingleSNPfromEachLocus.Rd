\name{TakeSingleSNPfromEachLocus}
\alias{TakeSingleSNPfromEachLocus}
\title{Take a Single SNP from Each Locus}
\description{This function will semi-randomly prune SNPs from a single locus}
\usage{
TakeSingleSNPfromEachLocus(data)
}
\arguments{
  \item{data}{A SNP dataset}
}
\details{For some analyses, SNPs should be independent of one another and so should not come from the same locus. This function will semi-randomly prune a single SNP from one locus by first determining which SNP has the most data first and then from equally good choices choose a random site. Note that if data come in as a concatenated sequence, it will look like a single locus and only choose one SNP.}
\value{Returns a list of three: [[1]] is the final reduced data matrix, [[2]] is the SNP sites per locus that were kept, and [[3]] is a vector of whether it was taken from using max amount of data ("M") or random ("R").
}
\examples{
  data(fakeData)
  TakeSingleSNPfromEachLocus(fakeData)
}
\author{B.Banbury}
