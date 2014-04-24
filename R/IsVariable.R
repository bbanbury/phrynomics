IsVariable <- function(SNP){
# Checks ambiguous sites for certain unique bases
# AAAM will return FALSE (because M can be either A or C--not unique)
# But AAAS will return T
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
