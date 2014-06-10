RemoveMissingSpeciesLoci <- function(data, chatty=FALSE) {
  speciesWithN <- apply(data, 2, MissingSpeciesVector, SpeciesNames=rownames(data))
  if(chatty)
    print(paste("removed", length(which(speciesWithN == "FALSE")), "of", length(speciesWithN), "loci"))
  data2 <- as.data.frame(data[,speciesWithN])
  rownames(data2) <- rownames(data)  
  return(data2) 
}
