#' Return MrBayes Ambiguity Codes
#' 
#' This function will take an IUPAC ambiguity code and return a set of bases in MrBayes parenthetical format
#' @param NucCode An ambiguity code
#' @export
#' @return Returns a character vector with numerical base possibilities surrounded by parentheses. SNPs will be transformed where A=1, T=2, G=3, and C=4.
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{ReturnAmbyCode}
#' @examples
#' ReturnMrBayesAmbyCode("A")
#' ReturnMrBayesAmbyCode("R")

ReturnMrBayesAmbyCode <- function(NucCode) {
  possibilities <- NULL
  if(NucCode == "A" || NucCode == "G" || NucCode == "C" || NucCode == "T" || NucCode == "U")
    return(ReturnMrBayesCode(NucCode))
  if(NucCode == "N" || NucCode == "-" || NucCode == "?") {
    return(NucCode)
  }
  if(NucCode == "R")  possibilities <- c("A", "G")
  if(NucCode == "Y")  possibilities <- c("C", "T")
  if(NucCode == "W")  possibilities <- c("A", "T")
  if(NucCode == "S")  possibilities <- c("G", "C")
  if(NucCode == "M")  possibilities <- c("A", "C")
  if(NucCode == "K")  possibilities <- c("G", "T")
  if(NucCode == "B")  possibilities <- c("G", "C", "T")
  if(NucCode == "H")  possibilities <- c("A", "C", "T")
  if(NucCode == "D")  possibilities <- c("A", "G", "T")
  if(NucCode == "V")  possibilities <- c("A", "G", "C")
  return(paste("(", paste(sapply(possibilities, ReturnMrBayesCode), collapse=""), ")", sep=""))
}
