#' Is Variable
#' 
#' This function will determine if a site is variable (TRUE) or not (FALSE). Ambiguity codes are taken into account to mean either heterozygotes or uncertainty. For example, c("A", "A", "S") will return TRUE, because "S" reduces to "G" and "C" and so either of those is variable with "A"; c("A", "A", "M") will return FALSE, because M can be either A or C. If the "M" is uncertain and is an "A" then it is not variable.  
#' @param SNP A single SNP, as a vector of bases
#' @export
#' @return Boolean response of variability
#' @seealso \link{ReadSNP} \link{WriteSNP}
#' @examples
#' IsVariable(c("A", "A", "A", "S"))
#' IsVariable(c("A", "A", "A", "M"))

IsVariable <- function(SNP){
  var <- FALSE
  bases <- c("A", "C", "G", "T", "U")
  basesInSNP <- bases[which(bases %in% SNP)]
  if(all(SNP == " "))
    return(TRUE)
  if(length(basesInSNP) == 0)  #all "N"
    return(FALSE)
  if(length(basesInSNP) > 1)  #more than one base anyway, so skip checking ambigs
    return(TRUE)
  if(length(basesInSNP) == 1) { #if only one base plus ambigs
    for(i in 2:length(SNP)) {
      var <- c(var, !any(ReturnNucs(SNP[i]) %in% basesInSNP))
    }
  }
  return(any(var))
}
