\name{PreProcess}
\alias{PreProcess}
\title{PreProcess SNP data}
\description{This function will "preprocess" a SNP dataset. }
\usage{
PreProcess(SNPdataset)
}
\arguments{
  \item{SNPdataset}{A SNP dataset}
}
\details{SNP data matrices that come out of other software programs sometimes have off characters that can not be read into R. THis function will strip these out of the matrix.}
\value{Returns a data frame without odd characters.  If you find your matrix still includes wonky characters, please contact me or fill out a bug report, so that I can add it to the list. 
}
\author{B.Banbury}
