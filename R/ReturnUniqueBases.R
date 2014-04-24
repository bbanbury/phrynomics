ReturnUniqueBases <- function(SNP){
#function only works with binary data
  bases <- unique(c(sapply(SNP, ReturnNucs, forSNAPP=TRUE), recursive=T)) #forSNAPP arg
  if(any(bases == "-"))
    bases <- bases[-which(bases == "-")]
  return(bases)
}
