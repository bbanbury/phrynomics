IsBinary <- function(SNP){
  if(all(SNP == " "))
    return (TRUE)
  bases <- unique(c(sapply(SNP, ReturnNucs, forSNAPP=TRUE), recursive=T)) #forSNAPP arg
  if(any(bases == "-"))  #remove any missing data
    bases <- bases[-which(bases == "-")]
  if(length(bases) == 2)  return(TRUE)
  else return(FALSE)
}
