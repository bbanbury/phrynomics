% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WriteSNP.R
\name{WriteSNP}
\alias{WriteSNP}
\title{Write SNP}
\usage{
WriteSNP(SNPdataset, file = "", format = "phylip", missing = "N")
}
\arguments{
\item{SNPdataset}{A data matrix of the class "matrix", "data.frame", or "snp"}

\item{file}{A file name specified by either a variable or a double-quoted string. If not specified (ie, file=""), then  data is printed to the console.}

\item{format}{Format for file, either "nexus" or "phylip"}

\item{missing}{A character denoting missing data, usually "N", "?", or "-"}
}
\description{
Write a SNP File
}
\details{
This function writes SNP datasets into nexus or phylip formatted files.
}
\examples{
a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
write.table(a, file="ex.data", col.names=FALSE)
b <- ReadSNP("ex.data")
WriteSNP(b, file="refakeData.phy")

data(fakeData)
WriteSNP(fakeData, file="refakeData.nex", format="nexus")
}
\seealso{
\link{ReadSNP}
}
