#' Return Ambiguity Codes
#' 
#' This function will take an IUPAC ambiguity code and return a set of bases
#' @param NucCode An ambiguity code
#' @param forSNAPP Logical. If FALSE (default), then missing data characters will be returned as all possibilities. If TRUE, the return for missing data will be returned "-".
#' @export
#' @return Returns a character vector with base possibilities.
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{ReturnAmbyCode}
#' @examples
#' ReturnNucs("N", forSNAPP=FALSE)
#' ReturnNucs("N", forSNAPP=TRUE)
#' ReturnNucs("K")

ReturnNucs <- function(NucCode, forSNAPP=FALSE) {
  possibilities <- NULL
  if(NucCode == "A" || NucCode == "G" || NucCode == "C" || NucCode == "T" || NucCode == "U")
    possibilities <- NucCode
  if(NucCode == "N" || NucCode == "-" || NucCode == "?") {
    if(forSNAPP)  possibilities <- "-"
    else  possibilities <- c("A", "G", "C", "T", "U")
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
  return(possibilities)
}
