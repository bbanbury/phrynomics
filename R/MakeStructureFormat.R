##' Make Structure Formatted Files
#' 
#' This function will create tables that are needed as Structure input. Alleles are coded
#' numerically (A=1, C=2, G=3, T=4, NA=-9). Each sample's unlinked file gets reduced to 
#' its individual alleles (i.e. ambiguity codes are assumed to be heterozygous snps).
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @export
#' @return Returns a table with alleles. 
#' @seealso \link{ReadSNP} \link{AddAFlag} \link{ReturnAlleleCounts}
#' @examples
#' data(fakeData)
#' fakeData <- RemoveNonBinary(RemoveInvariantSites(fakeData))
#' MakeStructureFormat(fakeData)
MakeStructureFormat <- function(SNPdataset){
  snpclass <- "table"
  if (inherits(SNPdataset, "snp")) {
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  splitdata <- SplitSNP(SNPdataset)
  newSNPdataset <- NULL
  for(i in sequence(dim(splitdata)[1])){
    uniques <- sapply(splitdata[i,], ReturnNucs)
    toreduce <- which(lapply(uniques, length) == 2)
    res <- rbind(splitdata[i,], splitdata[i,])
    if(length(toreduce > 0)){
      for(j in toreduce){
        res[,j] <- uniques[[j]]
      }
    }
    newSNPdataset <- rbind(newSNPdataset, res)
  }
  newSNPdataset <- TranslateBases(newSNPdataset)
  newSNPdataset <- SplitSNP(newSNPdataset)
  newSNPdataset  <- apply(newSNPdataset, 2, function(x) return(gsub("-", "-9", x)))
  return(newSNPdataset)
}
