\name{IsVariable}
\alias{IsVariable}
\title{Is Variable SNP}
\description{This function will detect whether a SNP is variable}
\usage{
IsVariable(SNP)
}
\arguments{
  \item{SNP}{A single SNP}
}
\details{This function will determine if a SNP site is variable (TRUE) or not (FALSE). Ambiguity codes are taken into account to mean either heterozygotes or uncertainty. For example, c("A", "A", "S") will return TRUE; c("A", "A", "M") will return FALSE, because M can be either A or C. 
}
\value{Boolean return for whether SNP is variable
}
\author{B.Banbury}
