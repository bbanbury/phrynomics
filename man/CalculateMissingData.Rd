\name{CalculateMissingData}
\alias{CalculateMissingData}
\title{Calculate Missing Data}
\description{This function will calculate the amount of missing data in an SNP matrix}
\usage{
  CalculateMissingData(data, missingChar="N")
}

\arguments{
  \item{data}{SNP dataset}
  \item{missingChar}{Character that represents missing data ("N", "?", "-")}
}
\examples{
  data(fakeData)
  CalculateMissingData(fakeData, missingChar="N")
}
\author{B.Banbury}