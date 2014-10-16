#' Return Ambiguity Codes
#' 
#' This function will take a set of bases and return the IUPAC ambiguity code. 
#' @param bases A vector of bases
#' @export
#' @return Returns a single character with the ambiguity code.
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{ReturnNucs}
#' @examples
#' bases <- c("A", "G")
#' ReturnAmbyCode(bases)
#' ReturnAmbyCode(c("A","G","C"))

ReturnAmbyCode <- function(bases){
  if(all(bases %in% c("A","G")))  return("R")
  if(all(bases %in% c("C","T")))  return("Y")
  if(all(bases %in% c("A","T")))  return("W")
  if(all(bases %in% c("G","C")))  return("S")
  if(all(bases %in% c("A","C")))  return("M")
  if(all(bases %in% c("G","T")))  return("K")
  if(all(bases %in% c("G","C","T")))  return("B")
  if(all(bases %in% c("A","C","T")))  return("H")
  if(all(bases %in% c("A","G","T")))  return("D")
  if(all(bases %in% c("A","G","C")))  return("V")
  if(all(bases %in% c("A","G","C","T")))  return("N")
  if(all(bases %in% c("A","G","C","T","U")))  return("N")
}
