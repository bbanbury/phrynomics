convertPhyDataToNexData <- function(phyData) {
  phyData <- SplitSNP(phyData)
  if(any(phyData[1,] == " "))
    phyData <- SplitSNP(phyData)[,-which(SplitSNP(phyData)[1,] == " ")]
  nexData <- vector("list", dim(phyData)[1])
  for(i in 1:length(nexData)){
    nexData[[i]] <- phyData[i,] #nexus parser doesn't like spaces between loci
    names(nexData)[i] <- rownames(phyData)[i]
  }
  return(nexData)
}
