ReturnUniqueBases <- function(SNP){
  bases <- unique(c(sapply(SNP, ReturnNucs, forSNAPP=TRUE), recursive=T)) #forSNAPP arg
  if(any(bases == "-"))
    bases <- bases[-which(bases == "-")]
  return(bases)
}
