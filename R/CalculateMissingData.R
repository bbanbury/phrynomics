#' Calculate Amount of Missing Data
#' 
#' This function will calculate what percentage of the dataset is missing.  You must give it the missing character value (ex: "N", "-"), and you can not have more than one value in the dataset. 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param calc Calculate missing data either by individual "sites" or "loci". If you choose to calculate missing data by loci, any site with information will count (for example, "NNNNA" will count as non-missing)
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @export
#' @return Returns a vector of missing data for each site or locus
#' @seealso \link{ReadSNP} \link{SplitSNP}
#' @examples
#' data(fakeData)
#' CalculateMissingData(fakeData, calc="sites")
#' CalculateMissingData(fakeData, calc="loci")

CalculateMissingData <- function(SNPdataset, calc="sites", missingChar="any"){
  if (missingChar == "any")
    missingChar <- c("N", "-", "?")
  char <- match.arg(arg=missingChar, choices=c("N", "-", "?"), several.ok=TRUE)
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data 
  missingSites <- MakePresentAbsent(SNPdataset, calc="loci", returnData="presabsentdata")
  PercentMissing <- apply(missingSites, 2, function(x) (length(x) - sum(x))/length(x))
  return(PercentMissing)
}
