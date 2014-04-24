RemoveMissingSpeciesLoci <- function(data) {
  speciesWithN <- apply(data, 2, MissingSpeciesVector, SpeciesNames=rownames(data))
  print(paste("removed", length(which(speciesWithN == "FALSE")), "of", length(speciesWithN), "loci"))
  data2 <- as.data.frame(data[,speciesWithN])
  rownames(data2) <- rownames(data)  
  return(data2) 
}
