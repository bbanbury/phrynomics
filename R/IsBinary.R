#' Is Binary
#' 
#' This function will determine if a site is binary (TRUE) or not (FALSE). For example, c("A", "A", "T") will return TRUE; c("A", "A", "A") or c("A", "G", "T") will return FALSE. Ambiguity codes are taken into account to mean either heterozygotes or uncertainty. Such that c("A", "A", "R") will return TRUE; c("A", "A", "Y") will return FALSE.
#' @param site A single SNP site, as a vector of bases
#' @export
#' @return Boolean response whether binary
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{IsVariable}
#' @examples
#' IsBinary(c("A", "A", "R"))
#' IsBinary(c("A", "A", "Y"))

IsBinary <- function(site){
  if(all(site == " "))
    return (TRUE)
  bases <- unique(c(sapply(site, ReturnNucs, forSNAPP=TRUE), recursive=T)) #forSNAPP arg
  if(any(bases == "-"))  #remove any missing data
    bases <- bases[-which(bases == "-")]
  if(length(bases) == 2)  return(TRUE)
  else return(FALSE)
}
