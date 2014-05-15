\name{GetBaseFrequencies}
\alias{GetBaseFrequencies}
\title{Get Base Frequencies}
\description{This function will get base frequencies for a single SNP}
\usage{
GetBaseFrequencies(SNP)
}
\arguments{
  \item{SNP}{A single SNP site}
}
\value{Returns a named vector with bases and frequencies}
\examples{
  data(fakeData)
  splitData <- SplitSNP(fakeData)  #the GetBaseFrequencies function will only work on single SNPs
  GetBaseFrequencies(splitData[,1])
  GetBaseFrequencies(splitData[,2])
  apply(splitData, 2, GetBaseFrequencies)
}
\author{B.Banbury}
