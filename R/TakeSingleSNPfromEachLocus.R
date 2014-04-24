TakeSingleSNPfromEachLocus <- function(data) { 
#this function will take a SNP from each locus by 1)solving for the SNP with the least amount of missing data and 2)taking a random SNP from any ties.
  nuSites <- GetNumberOfSitesForLocus(data)
  singleSNPfromLocus <- matrix(nrow=dim(data)[1], ncol=length(nuSites))
  rownames(singleSNPfromLocus) <- rownames(data)
  whichRandomSites <- NULL
  randomOrMax <- NULL
  for(locus in sequence(length(nuSites))){
    singleLocus <- as.matrix(data[,locus])
    keepSNPs <- which(CalculateMissingData(singleLocus) == min(CalculateMissingData(singleLocus)))
    if(length(keepSNPs) > 1) {  #if there are ties, then take a random SNP
      randomSite <- floor(runif(1, min=1, max=1+length(keepSNPs)))
      whichRandomSites <- c(whichRandomSites, randomSite)
      randomOrMax <- c(randomOrMax, "R")
      singleSNPfromLocus[,locus] <- SplitSNP(singleLocus)[, randomSite]    
    }
    else {  #this is if there is only one SNP with the most data
      whichRandomSites <- c(whichRandomSites, keepSNPs)
      randomOrMax <- c(randomOrMax, "M")
      singleSNPfromLocus[,locus] <- SplitSNP(singleLocus)[, keepSNPs]
    }
  }
  singleSNPfromLocus <- cSNP(singleSNPfromLocus)
  return(list(singleSNPfromLocus, whichRandomSites, randomOrMax))
}
