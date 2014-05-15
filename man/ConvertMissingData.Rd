\name{ConvertMissingData}
\alias{ConvertMissingData}
\title{Convert Missing Data}
\description{This function will convert missing data characters}
\usage{
ConvertMissingData(SNPdataset, oldmissing="N", newmissing="?")
}
\arguments{
  \item{SNPdataset}{The SNP data frame}
  \item{oldmissing}{The old character to replace}
  \item{newmissing}{The new character}
}
\details{This function will return the same data frame as input, but with replacement characters}
\value{Returns a data frame}
\examples{
  data(fakeData)
  ConvertMissingData(fakeData, "-", "N")
}
\author{B.Banbury}
