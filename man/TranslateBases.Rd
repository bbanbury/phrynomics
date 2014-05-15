\name{TranslateBases}
\alias{TranslateBases}
\title{Translate Bases}
\description{This function will translate the bases in various ways}
\usage{
TranslateBases(data, translateMissing=TRUE, translateMissingChar="-", ordered=FALSE)
}
\arguments{
  \item{data}{SNP dataset}
  \item{translateMissing}{Whether to translate missing characters}
  \item{translateMissingChar}{If \code{translateMissing} is TRUE, then it will convert all "N", "?", and "-" to this character}
  \item{ordered}{Whether character states should be ordered by frequency (used for SNAPP)}
}
\details{Ordered transformations will only transform binary states; the most common allele gets transformed to a 0, the less common allele to a 2, and the heterozygote to a 1. This is used for SNAPP analyses.  For MrBayes Mk/Mkv model, the states do not need to be ordered nor binary}
\value{Returns a data frame}
\examples{
  data(fakeData)
  
  #for MrBayes Mkv model
  TranslateBases(fakeData)
}
\author{B.Banbury}
