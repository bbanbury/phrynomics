##' Take Single SNP From Each Locus
#' 
#' This function will keep a single site from each locus in the dataset.  For example, if your dataset contains 450 sites in 3 loci, your final dataset will only include 3 sites. This is to prevent linkage issues across sites within the same locus. This function will initially caluclate the percentage of species that have data for each site within a locus and it will choose to keep the site with the most data.  If there is equally good coverage across sites, then one is chosen by random. 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @export
#' @return Returns a list of three items.  The first is the subset dataset with only a single site per locus. The second is a vector of the position of the chosen site within the locus. The third is a vector that details if the site was chosen based on maximum coverage ("M") or randomly ("R"). 
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{CalculateMissingData}
#' @examples
#' data(fakeData)
#' TakeSingleSNPfromEachLocus(fakeData)
#' TakeSingleSNPfromEachLocus(fakeData)[[1]]

TakeSingleSNPfromEachLocus <- function(SNPdataset) { 
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  nuSites <- GetNumberOfSitesForLocus(SNPdataset)
  singleSNPfromLocus <- matrix(nrow=dim(SNPdataset)[1], ncol=length(nuSites))
  rownames(singleSNPfromLocus) <- rownames(SNPdataset)
  whichRandomSites <- NULL
  randomOrMax <- NULL
  for(locus in sequence(length(nuSites))){
    singleLocus <- as.matrix(SNPdataset[,locus])
    keepSNPs <- which(CalculateMissingData(singleLocus, "sites") == min(CalculateMissingData(singleLocus, "sites")))
    if(length(keepSNPs) > 1) {
      randomSite <- floor(runif(1, min=1, max=1+length(keepSNPs)))
      whichRandomSites <- c(whichRandomSites, randomSite)
      randomOrMax <- c(randomOrMax, "R")
      singleSNPfromLocus[,locus] <- SplitSNP(singleLocus)[, randomSite]    
    }
    else {
      whichRandomSites <- c(whichRandomSites, keepSNPs)
      randomOrMax <- c(randomOrMax, "M")
      singleSNPfromLocus[,locus] <- SplitSNP(singleLocus)[, keepSNPs]
    }
  }
  singleSNPfromLocus <- cSNP(singleSNPfromLocus)
  if(snpclass == "snp")
    return(list(snpdata=ReadSNP(singleSNPfromLocus), position=whichRandomSites, randomOrMax=randomOrMax))
  return(list(snpdata=singleSNPfromLocus, position=whichRandomSites, randomOrMax=randomOrMax))
}
