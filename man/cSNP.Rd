\name{cSNP}
\alias{cSNP}
\alias{SplitSNP}
\title{Concatenate or Split SNP Datasets}
\description{This function will concatenate or split datasets, and reduce them if necessary by the keepVector and/or remove spaces that indicate loci}
\usage{
  cSNP(splitSNP, KeepVector=NULL, maintainLoci=TRUE)
  SplitSNP(SNPdataset)
}
\arguments{
  \item{splitSNP}{A split SNP data matrix (using SplitSNP function)}
  \item{KeepVector}{Boolean vector to keep or remove individual SNPs}
  \item{maintainLoci}{If TRUE, it will keep any spaces in the dataset which indicate locus breaks}
  \item{SNPdataset}{An SNP dataset}
}
\details{This function will return a data frame of concatenated or split SNPs}
\value{Returns a data frame}
\examples{
  data(fakeData)
  splitData <- SplitSNP(fakeData)
  cSNP(splitData, maintainLoci=TRUE)

  cSNP(splitData, maintainLoci=FALSE)
}
\author{B.Banbury}

