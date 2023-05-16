#' Reduce Data by Minimum Individuals
#' 
#' This function will reduce a dataset to include sites or loci that meet the minimum number of individuals.  For example, if you want only data for which at least 75% of the individuals in your dataset are represented, then it will remove those loci/sites that have less.     
#'
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param calc Calculate missing data either by individual "sites" or "loci". If you choose to calculate missing data by loci, any site with information will count (for example, "NNNNA" will count as non-missing)
#' @param threshold The minimum number of individuals that must be present for data to be retained.  This can be either a percent of individuals or a whole integer. When threshold=1, it means at least one individual, not 100\% of the data.  
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @param chatty Optional print to screen messages
#' @export
#' @return Returns a subset dataset.  
#' @seealso \link{ReadSNP} \link{MakeMinIndTable} \link{CalculateMissingData} \link{plotMissing}
#' @examples
#' data(fakeData)
#' MakePresentAbsent(fakeData, "sites")
#' MakePresentAbsent(fakeData, "loci")

ReduceMinInd <- function(SNPdataset, calc="loci", threshold=0.9, missingChar="any", chatty=FALSE){
  snpclass <- "table"
  if(inherits(SNPdataset, "snp")){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  initialsites <- sum(nchar(SNPdataset[1,]))
  if(threshold <= 0.99)
    threshold <- ceiling((threshold * dim(SNPdataset)[1]))
  if(chatty)
    message(paste(calc, "must have", threshold, "individuals to be retained."))
  if(calc == "sites")
    SNPdataset <- SplitSNP(SNPdataset)
  presabsentdata <- apply(SNPdataset, c(1,2), IsMissing, missingChar=missingChar)
  tots <- apply(presabsentdata, 2, sum)
  keepVector <- tots >= threshold
  if(calc == "sites"){
    newSNPdataset <- cSNP(SNPdataset, keepVector, maintainLoci=TRUE)
    if(chatty)
      message(paste("removed", initialsites - sum(nchar(newSNPdataset[1,])), "of", initialsites, calc))
  }
  if(calc == "loci"){
    newSNPdataset <- SNPdataset[which(keepVector)]
    if(chatty)
      message(paste("removed", dim(SNPdataset)[2] - dim(newSNPdataset)[2], "of", dim(SNPdataset)[2], calc))
  }
  if(snpclass == "snp")
    return(ReadSNP(newSNPdataset))
  else
    return(newSNPdataset)
}