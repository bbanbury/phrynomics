\name{ReturnMrBayesAmbyCode}
\alias{ReturnMrBayesAmbyCode}
\alias{ReturnMrBayesCode}
\title{Return MrBayes Ambiguity Code}
\description{This function will return MrBayes Mk/Mkv model ambiguity codes}
\usage{
ReturnMrBayesCode(base)
ReturnMrBayesAmbyCode(NucCode)
}
\arguments{
  \item{NucCode}{A single SNP Base or ambiguity code}
  \item{base}{A single SNP Base, A, T, G, or C}
}
\details{MrBayes Mk/Mkv models use numeric character assignments, so this function will transform a data A, T, G, C into 1, 2, 3, 4's respectively. Ambiguous assignemnts are coded between parenthesis. See http://mrbayes.sourceforge.net/wiki/index.php/Evolutionary_Models_Implemented_in_MrBayes_3#Standard_Discrete_.28Morphology.29_Model}
\value{Returns a single character
}
\author{B.Banbury}
\examples{
ReturnMrBayesAmbyCode("A")
ReturnMrBayesAmbyCode("D")
}