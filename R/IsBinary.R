#' Is Binary
#' 
#' This function will determine if a SNP site is binary (TRUE) or not (FALSE). For example, c("A", "A", "T") will return TRUE; c("A", "A", "A") or c("A", "G", "T") will return FALSE. Ambiguity codes are taken into account to mean either heterozygotes or uncertainty. Such that c("A", "A", "R") will return TRUE; c("A", "A", "Y") will return FALSE.
#' @param SNP A single SNP, as a vector of bases
#' @export
#' @return Boolean response whether binary
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{IsVariable}
#' @examples
#' IsBinary(c("A", "A", "R"))
#' IsBinary(c("A", "A", "Y"))

IsBinary <- function(SNP){
  if(all(SNP == " "))
    return (TRUE)
  bases <- unique(c(sapply(SNP, ReturnNucs, forSNAPP=TRUE), recursive=T)) #forSNAPP arg
  if(any(bases == "-"))  #remove any missing data
    bases <- bases[-which(bases == "-")]
  if(length(bases) == 2)  return(TRUE)
  else return(FALSE)
}
