\name{RemoveNonBinary}
\alias{RemoveNonBinary}
\title{RemoveNonBinary}
\description{This function will remove any non-binary SNPs from the dataset}
\usage{
RemoveNonBinary(data)
}
\arguments{
  \item{data}{A SNP dataset}
}
\value{Returns a subset data frame.}
\examples{
  data(fakeData)
  RemoveNonBinary(fakeData)
}
\author{B.Banbury}
