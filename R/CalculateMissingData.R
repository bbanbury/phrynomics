##' Calculate Amount of Missing Data
#' 
#' This function will calculate what percentage of the dataset is missing.  You must give it the missing character value (ex: "N", "-"), and you can not have more than one value in the dataset. 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @export
#' @return Returns a vector of missing data for each SNP
#' @seealso \link{ReadSNP} \link{ConvertMissingData}
#' @examples
#' data(fakeData)
#' CalculateMissingData(fakeData)

CalculateMissingData <- function(SNPdataset, missingChar="any"){
  if (missingChar == "any")
    missingChar <- c("N", "-", "?")
  char <- match.arg(arg=missingChar, choices=c("N", "-", "?"), several.ok=TRUE)
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data 
  splitdata <- SplitSNP(SNPdataset)
  MissingSites <- rep(0, dim(splitdata)[2])
  for(i in sequence(length(char))){
    MissingSites <- MissingSites + apply(splitdata, 2, function(x) length(which(x == char[i])))
  }
  PercentMissing <- MissingSites/dim(splitdata)[1]
  return(PercentMissing)
}
